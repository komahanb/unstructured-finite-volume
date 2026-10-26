!=====================================================================!
! Module that states the Gauss-Seidel iteration on the tower, swept by
! colour.
!
! Jacobi corrects every cell from the old state; Gauss-Seidel lets
! each correction read the corrections already made. On a graph the
! valid order is the colouring the sweep_order delegation already
! returns: all cells of one colour share no face, so a whole colour
! updates at once, each colour reading every colour before it,
!
!      for each colour:  r = rhs - A x       (x already partly new)
!                        x <- x + omega * r / diag   on that colour
!
!   one iteration, three colours
!
!   colour 1:   x(1) += omega diag^-1 r            r from the old x
!   colour 2:   r updated, x(2) += ...             reads colour 1's new x
!   colour 3:   r updated, x(3) += ...             reads colours 1 and 2
!   then:       the residual measured, the stop tested
!
! The residual measured at the top of an iteration is the first
! colour's, so a sweep of c colours costs c products with the
! operator. Which colouring is swept is determined by the type itself:
! jacobi is this iteration over one colour class, and states so by
! extension.
!
! SOR is not another solver. SOR is this one at omega not equal to
! one: a parameter, absorbed, as the admission law requires.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_gauss_seidel

  use iso_fortran_env    , only : int64
  use util_precision     , only : dp
  use operation_minimizer, only : minimizer
  use operation_stencil  , only : stencil
  use relation_binary    , only : group_by_key
  use util_factorisation , only : dense_factorisation
  use util_tally         , only : linear_solves, factorisations

  implicit none

  private

  public :: gauss_seidel

  !-------------------------------------------------------------------!
  ! Type that models the Gauss-Seidel iteration, a member of the
  ! minimizer hierarchy
  !-------------------------------------------------------------------!

  type, extends(minimizer) :: gauss_seidel

     real(dp) :: omega = 1.0_dp

     ! The diagonal, its factorisations and the colouring belong to the
     ! attached operator, not to one solve: computing them costs
     ! maxval(colours) * block_width operator applications, so they are
     ! stored until the operator is attached again
     real(dp)                 , allocatable, private :: stored_diagonal(:,:,:)
     type(dense_factorisation), allocatable, private :: stored_block(:)
     integer                  , allocatable, private :: stored_colours(:)

   contains

     procedure :: name            => gauss_seidel_name
     procedure :: colouring
     procedure :: storage_entries => gauss_seidel_storage_entries
     procedure :: solve

  end type gauss_seidel

  ! A Gauss-Seidel iteration stated by its sweeps
  interface gauss_seidel
     module procedure gauss_seidel_of
  end interface gauss_seidel

contains

  !===================================================================!
  ! Returns the name of the iteration
  !===================================================================!

  pure function gauss_seidel_name(this) result(name)

    class(gauss_seidel), intent(in)  :: this
    character(len=:)   , allocatable :: name

    ! Names this, which the name does not read, to keep the compiler
    ! from warning that the argument is unused
    associate (u1 => this); end associate

    name = 'gauss-seidel'

  end function gauss_seidel_name

  !===================================================================!
  ! Returns the colour classes the sweep runs over: one entry per
  ! unknown, from the coupling's own colouring
  !===================================================================!

  subroutine colouring(this, n, colours)

    class(gauss_seidel) , intent(in)  :: this
    integer             , intent(in)  :: n
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
  ! Returns the entries stored: the block diagonal, one square block of
  ! the stated width per block of unknowns, and its factorised copy
  ! where the width is above one
  !===================================================================!

  pure integer(int64) function gauss_seidel_storage_entries(this, num_unknowns) result(entries)

    class(gauss_seidel), intent(in) :: this
    integer            , intent(in) :: num_unknowns

    entries = int(num_unknowns, int64) * int(max(this % block_width, 1), int64)
    if (this % block_width > 1) entries = 2_int64 * entries

  end function gauss_seidel_storage_entries

  !===================================================================!
  ! Solves A x = rhs by sweeps over the colours, and returns the norm
  ! of the residual achieved.
  !
  ! An explicit operator, a stencil, keeps the residual current by
  ! subtracting the product of the columns each colour class changes;
  ! any other operator recomputes it by a full product before each
  ! class after the first:
  !
  !   stencil:          r <- r - A(:, changed columns) dx(changed)
  !   other operators:  r <- rhs - A x
  !
  ! The diagonal is taken a block at a time: a block of width one is
  ! the number itself, and a wider one is factorised once for the
  ! operator. Over several images every image takes the same branches,
  ! and exchanges after each colour of the whole.
  !===================================================================!

  subroutine solve(this, rhs, x, achieved)

    class(gauss_seidel), intent(inout) :: this
    real(dp)           , intent(in)    :: rhs(:)
    real(dp)           , intent(inout) :: x(:)
    real(dp)           , intent(out)   :: achieved

    real(dp), allocatable                  :: d(:,:,:), r(:), block_solution(:), dx(:)
    type(dense_factorisation), allocatable :: block(:)
    integer , allocatable                  :: colours(:), ordered(:), first(:)
    logical , allocatable                  :: active(:)
    integer                                :: it, col, b, w, nb, i, num_colours
    logical                                :: explicit, distributed, zero_start
    character(len=250)                     :: message

    call this % record_event(linear_solves)

    ! An explicit operator keeps the residual current by the product of
    ! the columns each colour class changes; any other recomputes it by
    ! a full product before each class
    explicit = .false.
    select type (a => this % action)
    type is (stencil)
       explicit = .true.
    end select

    call this % initialize_residual_history()
    ! The residual at the zero state is the right-hand side itself
    zero_start = all(x == 0.0_dp)
    ! Over several images every image takes the same branch, the one
    ! that exchanges or the one that does not
    if (allocated(this % distribution)) zero_start = this % distribution % maximum(merge(0.0_dp, 1.0_dp, zero_start)) == 0.0_dp
    if (zero_start) then
       r = rhs
    else
       call this % imbalance(rhs, x, r)
    end if
    achieved = this % norm(r)
    if (this % terminated(achieved, 0)) return

    ! The diagonal, a block at a time: a block of width one is the
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
             ! A block of another image's is never relaxed here, and its
             ! rows are not in this image's part
             if (allocated(this % distribution)) then
                if (.not. this % distribution % is_owned((b - 1) * w + 1)) cycle
             end if
             call this % record_event(factorisations)
             call block(b) % factorise(d(:, :, b), tiny(1.0_dp))
             if (block(b) % singular()) then
                write(message,'(a,i0,a,i0,a)') 'gauss_seidel: the diagonal block at position ', b, &
                     & ' of ', nb, ' is singular'
                error stop trim(message)
             end if
          end do
       end if

       ! Over several images, the colours of the whole blocks, the same on
       ! every image, restricted to this image's part
       if (allocated(this % distribution)) then
          colours = this % distribution % block_colour
       else
          call this % colouring(nb, colours)
       end if

       ! The diagonal, its factorisation and the colouring are stored for
       ! the operator and read in place by every solve on it
       call move_alloc(d, this % stored_diagonal)
       if (allocated(this % stored_block)) deallocate(this % stored_block)
       if (allocated(block)) call move_alloc(block, this % stored_block)
       call move_alloc(colours, this % stored_colours)
       this % diagonal_valid = .true.

    end if
    num_colours = maxval(this % stored_colours)
    allocate(block_solution(w))
    distributed = allocated(this % distribution)
    ! Every image exchanges after each colour of the whole
    if (distributed) num_colours = this % distribution % num_colours
    if (explicit) then
       allocate(dx(size(x)), source=0.0_dp)
       ! The columns whose change reaches an owned row, by colour: every
       ! column of the part, the owned and the halo
       allocate(active(size(x)), source=.true.)
       call columns_by_colour(this % stored_colours, w, num_colours, active, ordered, first)
    end if

    do it = 1, this % max_iterations

       do col = 1, num_colours
          if (col > 1 .and. .not. explicit) then
             call this % imbalance(rhs, x, r)
          end if
          ! Every block of this colour, corrected by its diagonal
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
             ! The owned increments of this colour sent and the halo's
             ! received, then subtracted from the owned rows, and the
             ! other parts' increments in the summed rows this part owns,
             ! in the same sum
             if (distributed) call this % distribution % gather(dx, colour=col, summed_into=r, factor=-1.0_dp)
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

  !===================================================================!
  ! Returns the active unknowns ordered by the colour of their block:
  ! those of colour col are ordered(first(col):first(col+1)-1),
  !
  !   colour       1               2               3
  !   ordered      [ . . . . ]     [ . . . ]       [ . . ]
  !                ^first(1)       ^first(2)       ^first(3)
  !===================================================================!

  pure subroutine columns_by_colour(colours, w, num_colours, active, ordered, first)

    integer             , intent(in)  :: colours(:), w, num_colours
    logical             , intent(in)  :: active(:)
    integer, allocatable, intent(out) :: ordered(:), first(:)

    integer :: j

    ! A column that is not active has key zero and belongs to no colour
    call group_by_key(num_colours, merge([(colours((j - 1) / w + 1), j = 1, size(active))], 0, active), &
         & [(j, j = 1, size(active))], first, ordered)

  end subroutine columns_by_colour

  !===================================================================!
  ! Subtracts from r the product of the explicit operator's columns
  ! with the increments: r = r - A(:, columns) dx(columns)
  !===================================================================!

  subroutine subtracted_columns(this, columns, dx, r)

    class(gauss_seidel), intent(in)    :: this
    integer            , intent(in)    :: columns(:)
    real(dp)           , intent(in)    :: dx(:)
    real(dp)           , intent(inout) :: r(:)

    select type (a => this % action)
    type is (stencil)
       call a % column_product(columns, dx, r, factor=-1.0_dp)
    end select

  end subroutine subtracted_columns

  !===================================================================!
  ! Creates a Gauss-Seidel iteration stated by the sweeps it takes; its
  ! blocks are those of the system it is paired with
  ! (operation_minimization)
  !===================================================================!

  function gauss_seidel_of(sweeps) result(this)
    integer, intent(in) :: sweeps
    type(gauss_seidel)  :: this
    this % max_iterations = sweeps
  end function gauss_seidel_of

end module operation_gauss_seidel
