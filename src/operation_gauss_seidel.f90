!=====================================================================!
! The gauss-seidel iteration on the tower, swept by colour.
!
! Jacobi corrects every cell from the old state; gauss-seidel lets
! each correction see the ones already made. On a graph the safe
! order is the colouring the sweep_order delegation already answers:
! all cells of one colour share no face, so a whole colour updates
! at once, each colour seeing every colour before it,
!
!      for each colour:  r = rhs - A x       (x already partly new)
!                        x <- x + omega * r / diag   on that colour
!
! The residual measured at the top of an iteration is the first
! colour's, so a sweep of c colours costs c products with the
! operator. Which colouring is swept is a question the type answers
! for itself - jacobi is this iteration over one colour class, and
! says so by extension.
!
! SOR is not another solver. It is this one at omega away from one -
! a parameter, absorbed, exactly as the admission law orders.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_gauss_seidel

  use util_precision  , only : dp
  use operation_minimization, only : minimizer
  use util_factorisation    , only : dense_factorisation
  use util_tally, only : tally_record, linear_solves

  implicit none

  private
  public :: gauss_seidel

  type, extends(minimizer) :: gauss_seidel

     real(dp) :: omega = 1.0_dp

     ! the diagonal, its factorisations and the colouring belong to the
     ! attached operator, not to one solve: probing them costs
     ! maxval(colours) * block_width operator applications, so they are
     ! held until the operator is attached again
     real(dp)                 , allocatable, private :: held_diagonal(:,:,:)
     type(dense_factorisation), allocatable, private :: held_block(:)
     integer                  , allocatable, private :: held_colours(:)

   contains

     procedure :: name => gauss_seidel_name
     procedure :: colouring
     procedure :: solve

  end type gauss_seidel

contains

  pure function gauss_seidel_name(this) result(name)

    class(gauss_seidel), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'gauss-seidel'

  end function gauss_seidel_name

  !===================================================================!
  ! The colour classes the sweep runs over: one entry per unknown,
  ! from the coupling's own colouring.
  !===================================================================!

  subroutine colouring(this, n, colours)

    class(gauss_seidel), intent(in)  :: this
    integer            , intent(in)  :: n
    integer, allocatable, intent(out) :: colours(:)

    call this % sweep_order(colours)

    if (size(colours) /= n) then
       error stop 'gauss_seidel: one colour per unknown'
    end if

  end subroutine colouring


  subroutine solve(this, rhs, x, achieved)

    class(gauss_seidel), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: d(:,:,:), r(:), piece(:)
    type(dense_factorisation), allocatable :: block(:)
    integer , allocatable :: colours(:)
    integer :: it, col, b, w, nb, i

    call tally_record(linear_solves)

    ! the diagonal, a block at a time: a block of width one is the
    ! number itself, and a wider one is factorised once for the operator
    w  = this % block_width
    nb = size(x) / w

    if (.not. this % diagonal_valid) then

       call this % block_diagonal(d)

       if (w == 1) then
          do b = 1, nb
             if (abs(d(1, 1, b)) < tiny(1.0_dp)) d(1, 1, b) = huge(1.0_dp)
          end do
       else
          allocate(block(nb))
          do b = 1, nb
             call block(b) % factorise(d(:, :, b), tiny(1.0_dp))
             if (block(b) % singular()) then
                error stop 'gauss_seidel: a block on the diagonal is singular'
             end if
          end do
       end if

       call this % colouring(nb, colours)

       this % held_diagonal = d
       if (allocated(block)) then
          if (allocated(this % held_block)) deallocate(this % held_block)
          allocate(this % held_block, source=block)
       end if
       this % held_colours   = colours
       this % diagonal_valid = .true.

    else

       d       = this % held_diagonal
       colours = this % held_colours
       if (allocated(this % held_block)) allocate(block, source=this % held_block)

    end if
    call this % begin_imbalance()

    do it = 1, this % max_iterations

       call this % imbalance(rhs, x, r)
       achieved = this % norm(r)
       if (this % halted(achieved, it)) return

       do col = 1, maxval(colours)
          if (col > 1) then
             call this % imbalance(rhs, x, r)
          end if
          do b = 1, nb
             if (colours(b) /= col) cycle
             if (w == 1) then
                x(b) = x(b) + this % omega * r(b) / d(1, 1, b)
             else
                call block(b) % substitute(r((b - 1) * w + 1:b * w), piece, transposed=.false.)
                do i = 1, w
                   x((b - 1) * w + i) = x((b - 1) * w + i) + this % omega * piece(i)
                end do
             end if
          end do
       end do

    end do

    call this % imbalance(rhs, x, r)
    achieved = this % norm(r)

  end subroutine solve

end module operation_gauss_seidel
