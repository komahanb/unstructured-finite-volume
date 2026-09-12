!=====================================================================!
! The gauss-seidel iteration on the tower, swept by colour.
!
! Jacobi corrects every cell from the old state; gauss-seidel lets
! each correction read the corrections already made. On a graph the
! valid order is the colouring the sweep_order delegation already
! returns: all cells of one colour share no face, so a whole colour
! updates at once, each colour reading every colour before it,
!
!      for each colour:  r = rhs - A x       (x already partly new)
!                        x <- x + omega * r / diag   on that colour
!
! The residual measured at the top of an iteration is the first
! colour's, so a sweep of c colours costs c products with the
! operator. Which colouring is swept is determined by the type
! itself - jacobi is this iteration over one colour class, and
! states so by extension.
!
! SOR is not another solver. SOR is this one at omega not equal to
! one - a parameter, absorbed, as the admission law requires.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_gauss_seidel

  use iso_fortran_env , only : int64
  use util_precision  , only : dp
  use operation_minimization, only : minimizer
  use util_factorisation    , only : dense_factorisation
  use util_tally, only : linear_solves, factorisations

  implicit none

  private
  public :: gauss_seidel

  type, extends(minimizer) :: gauss_seidel

     real(dp) :: omega = 1.0_dp

     ! the diagonal, its factorisations and the colouring belong to the
     ! attached operator, not to one solve: computing them costs
     ! maxval(colours) * block_width operator applications, so they are
     ! stored until the operator is attached again
     real(dp)                 , allocatable, private :: stored_diagonal(:,:,:)
     type(dense_factorisation), allocatable, private :: stored_block(:)
     integer                  , allocatable, private :: stored_colours(:)

   contains

     procedure :: name => gauss_seidel_name
     procedure :: colouring
     procedure :: storage_entries => gauss_seidel_storage_entries
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

  !===================================================================!
  ! The block diagonal, one square block of the stated width per
  ! block of unknowns, and its factorised copy where wider than one.
  !===================================================================!

  pure integer(int64) function gauss_seidel_storage_entries(this, num_unknowns) result(entries)

    class(gauss_seidel), intent(in) :: this
    integer            , intent(in) :: num_unknowns

    entries = int(num_unknowns, int64) * int(max(this % block_width, 1), int64)
    if (this % block_width > 1) entries = 2_int64 * entries

  end function gauss_seidel_storage_entries


  subroutine solve(this, rhs, x, achieved)

    class(gauss_seidel), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    real(dp), allocatable :: d(:,:,:), r(:), block_solution(:)
    type(dense_factorisation), allocatable :: block(:)
    integer , allocatable :: colours(:)
    integer :: it, col, b, w, nb, i

    call this % record_event(linear_solves)

    call this % initialize_residual_history()
    call this % imbalance(rhs, x, r)
    achieved = this % norm(r)
    if (this % terminated(achieved, 0)) return

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
             call this % record_event(factorisations)
             call block(b) % factorise(d(:, :, b), tiny(1.0_dp))
             if (block(b) % singular()) then
                error stop 'gauss_seidel: a block on the diagonal is singular'
             end if
          end do
       end if

       call this % colouring(nb, colours)

       this % stored_diagonal = d
       if (allocated(block)) then
          if (allocated(this % stored_block)) deallocate(this % stored_block)
          allocate(this % stored_block, source=block)
       end if
       this % stored_colours   = colours
       this % diagonal_valid = .true.

    else

       d       = this % stored_diagonal
       colours = this % stored_colours
       if (allocated(this % stored_block)) allocate(block, source=this % stored_block)

    end if
    do it = 1, this % max_iterations

       do col = 1, maxval(colours)
          if (col > 1) then
             call this % imbalance(rhs, x, r)
          end if
          do b = 1, nb
             if (colours(b) /= col) cycle
             if (w == 1) then
                x(b) = x(b) + this % omega * r(b) / d(1, 1, b)
             else
                call block(b) % substitute(r((b - 1) * w + 1:b * w), block_solution, transposed=.false.)
                do i = 1, w
                   x((b - 1) * w + i) = x((b - 1) * w + i) + this % omega * block_solution(i)
                end do
             end if
          end do
       end do

       call this % imbalance(rhs, x, r)
       achieved = this % norm(r)
       if (this % terminated(achieved, it)) return

    end do

  end subroutine solve

end module operation_gauss_seidel
