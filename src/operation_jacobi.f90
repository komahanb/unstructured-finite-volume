!=====================================================================!
! Module that states the Jacobi iteration on the tower.
!
! Every cell is corrected by its own residual divided by its own
! diagonal, all cells at once:
!
!      x  <-  x + omega * (rhs - A x) / diag
!
! This is the Gauss-Seidel sweep over one colour class. With every
! cell the same colour, every correction reads the old state and none
! reads another correction:
!
!   gauss-seidel, two colours              jacobi, one colour
!
!   o x o x o    o corrected first;        o o o o o    every cell
!   x o x o x    x then reads the          o o o o o    corrected at
!   o x o x o    new o                     o o o o o    once from the
!                                                        old state
!
! The module states this by extension. The sweep, the diagonal and
! the evaluation of each residual are Gauss-Seidel's; this module
! specifies only which colouring is swept.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_jacobi

  use operation_gauss_seidel, only : gauss_seidel

  implicit none

  private

  public :: jacobi

  !-------------------------------------------------------------------!
  ! Type that models the Jacobi iteration: Gauss-Seidel over one
  ! colour class
  !-------------------------------------------------------------------!

  type, extends(gauss_seidel) :: jacobi

   contains

     procedure :: name      => jacobi_name
     procedure :: colouring => one_colour

  end type jacobi

contains

  !===================================================================!
  ! Returns the name of the iteration
  !===================================================================!

  pure function jacobi_name(this) result(name)

    class(jacobi)   , intent(in)  :: this
    character(len=:), allocatable :: name

    ! Names this, which the name does not read, to keep the compiler
    ! from warning that the argument is unused
    associate (u1 => this); end associate
    name = 'jacobi'

  end function jacobi_name

  !===================================================================!
  ! Returns one colour class for n unknowns: every unknown corrected at
  ! once from the state the iteration began at,
  !
  !   colours = [1, 1, 1, ..., 1]
  !
  ! The sweep does not read the coupling; the block diagonal of an
  ! operator other than a stencil still does.
  !===================================================================!

  subroutine one_colour(this, n, colours)

    class(jacobi)       , intent(in)  :: this
    integer             , intent(in)  :: n
    integer, allocatable, intent(out) :: colours(:)

    associate (u1 => this); end associate
    allocate(colours(n))
    colours = 1

  end subroutine one_colour

end module operation_jacobi
