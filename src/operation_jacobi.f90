!=====================================================================!
! The jacobi iteration on the tower.
!
! The oldest idea in iterative solving: correct every cell by its
! own residual over its own diagonal, all cells at once,
!
!      x  <-  x + omega * (rhs - A x) / diag
!
! which is the gauss-seidel sweep over one colour class - every cell
! the same colour, so every correction sees the old state and none
! sees another's. That is stated by extension: the sweep, the
! diagonal, the judgement of each residual are gauss-seidel's, and
! this file says only which colouring is swept.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_jacobi

  use operation_gauss_seidel, only : gauss_seidel

  implicit none

  private
  public :: jacobi

  type, extends(gauss_seidel) :: jacobi

   contains

     procedure :: name      => jacobi_name
     procedure :: colouring => one_colour

  end type jacobi

contains

  pure function jacobi_name(this) result(name)

    class(jacobi), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'jacobi'

  end function jacobi_name

  !===================================================================!
  ! One colour class: every unknown corrected at once from the state
  ! the iteration began at. The coupling is never asked.
  !===================================================================!

  subroutine one_colour(this, n, colours)

    class(jacobi), intent(in)  :: this
    integer      , intent(in)  :: n
    integer, allocatable, intent(out) :: colours(:)

    associate (u1 => this); end associate
    allocate(colours(n))
    colours = 1

  end subroutine one_colour

end module operation_jacobi
