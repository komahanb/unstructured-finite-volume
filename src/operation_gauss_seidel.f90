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
  use operation_stencil     , only : stencil
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

    character(len=250) :: message

    call this % sweep_order(colours)

    if (size(colours) /= n) then
       write(message,'(a,i0,a,i0)') 'gauss_seidel: one colour is required per unknown; &
            &size(colours) = ', size(colours), ', n = ', n
       error stop trim(message)
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

    real(dp), allocatable :: d(:,:,:), r(:), block_solution(:), dx(:)
    type(dense_factorisation), allocatable :: block(:)
    integer , allocatable :: colours(:), ordered(:), first(:)
    logical , allocatable :: active(:)
    integer :: it, col, b, w, nb, i, num_colours
    logical :: explicit, distributed
    character(len=250) :: message

    call this % record_event(linear_solves)

    ! an explicit operator keeps the residual current by the product
    ! of the columns each colour class changes; any other recomputes
    ! it by a full product before each class
    explicit = .false.
    select type (a => this % action)
    type is (stencil)
       explicit = .true.
    end select

    call this % initialize_residual_history()
    ! the residual at the zero state is the right-hand side itself
    if (all(x == 0.0_dp)) then
       r = rhs
    else
       call this % imbalance(rhs, x, r)
    end if
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
                write(message,'(a,i0,a,i0)') 'gauss_seidel: the diagonal block at position ', b, &
                     & ' of ', nb, ' is singular'
                error stop trim(message)
             end if
          end do
       end if

       call this % colouring(nb, colours)

       ! the diagonal, its factorisation and the colouring are stored
       ! for the operator and read in place by every solve on it
       call move_alloc(d, this % stored_diagonal)
       if (allocated(this % stored_block)) deallocate(this % stored_block)
       if (allocated(block)) call move_alloc(block, this % stored_block)
       call move_alloc(colours, this % stored_colours)
       this % diagonal_valid = .true.

    end if
    num_colours = maxval(this % stored_colours)
    allocate(block_solution(w))
    distributed = allocated(this % distribution)
    if (explicit) then
       allocate(dx(size(x)), source=0.0_dp)
       ! the columns whose change reaches an owned row, by colour:
       ! every column on one image, the owned and the halo over
       ! several
       allocate(active(size(x)), source=.true.)
       if (distributed) then
          active = this % distribution % is_owned
          active(this % distribution % halo) = .true.
       end if
       call columns_by_colour(this % stored_colours, w, num_colours, active, ordered, first)
    end if

    do it = 1, this % max_iterations

       do col = 1, num_colours
          if (col > 1 .and. .not. explicit) then
             call this % imbalance(rhs, x, r)
          end if
          do b = 1, nb
             if (this % stored_colours(b) /= col) cycle
             if (distributed) then
                if (.not. this % distribution % is_owned((b - 1) * w + 1)) cycle
             end if
             if (w == 1) then
                block_solution(1) = r(b) / this % stored_diagonal(1, 1, b)
             else
                call this % stored_block(b) % substituted(r((b - 1) * w + 1:b * w), block_solution, &
                     & transposed=.false.)
             end if
             do i = 1, w
                x((b - 1) * w + i) = x((b - 1) * w + i) + this % omega * block_solution(i)
                if (explicit) dx((b - 1) * w + i) = this % omega * block_solution(i)
             end do
          end do
          if (explicit) then
             if (distributed) call this % distribution % update(dx)
             associate (columns => ordered(first(col):first(col + 1) - 1))
               call subtracted_columns(this, columns, dx, r)
               dx(columns) = 0.0_dp
             end associate
          end if
       end do

       if (.not. explicit) call this % imbalance(rhs, x, r)
       achieved = this % norm(r)
       if (this % terminated(achieved, it)) return

    end do

  end subroutine solve

  ! the active unknowns ordered by the colour of their block: those of
  ! colour col are ordered(first(col):first(col+1)-1)
  pure subroutine columns_by_colour(colours, w, num_colours, active, ordered, first)

    integer, intent(in)  :: colours(:), w, num_colours
    logical, intent(in)  :: active(:)
    integer, allocatable, intent(out) :: ordered(:), first(:)

    integer, allocatable :: at(:)
    integer :: b, i, col

    allocate(first(num_colours + 1), source=0)
    do b = 1, size(colours)
       do i = 1, w
          if (active((b - 1) * w + i)) first(colours(b) + 1) = first(colours(b) + 1) + 1
       end do
    end do
    first(1) = 1
    do col = 1, num_colours
       first(col + 1) = first(col + 1) + first(col)
    end do
    at = first(1:num_colours)
    allocate(ordered(first(num_colours + 1) - 1))
    do b = 1, size(colours)
       col = colours(b)
       do i = 1, w
          if (.not. active((b - 1) * w + i)) cycle
          ordered(at(col)) = (b - 1) * w + i
          at(col) = at(col) + 1
       end do
    end do

  end subroutine columns_by_colour

  ! r = r - A(:, columns) dx(columns), on the explicit operator
  subroutine subtracted_columns(this, columns, dx, r)

    class(gauss_seidel), intent(in)    :: this
    integer            , intent(in)    :: columns(:)
    real(dp)           , intent(in)    :: dx(:)
    real(dp)           , intent(inout) :: r(:)

    select type (a => this % action)
    type is (stencil)
       if (allocated(this % distribution)) then
          call a % column_product(columns, dx, r, factor=-1.0_dp, rows=this % distribution % is_owned)
       else
          call a % column_product(columns, dx, r, factor=-1.0_dp)
       end if
    end select

  end subroutine subtracted_columns

end module operation_gauss_seidel
