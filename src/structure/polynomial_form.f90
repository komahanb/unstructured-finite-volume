!=====================================================================!
! The polynomial graph_form: the constant and the three coordinates,
! reckoned about the point of interest - the Taylor shape at degree
! one, whose span is every linear field.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module structure_polynomial_form

  use iso_fortran_env    , only : dp => REAL64
  use structure_graph, only : GRAPH_SIDE_VERTEX
  use structure_support, only : support
  use structure_graph_form        , only : graph_form

  implicit none

  private
  public :: polynomial_form

  type, extends(graph_form) :: polynomial_form

   contains

     procedure :: size_of => polynomial_size
     procedure :: values  => polynomial_values
     procedure :: slopes  => polynomial_slopes

  end type polynomial_form

  interface polynomial_form
     module procedure create
  end interface polynomial_form

contains

  ! Born with every table entry standing: the members are the four.
  pure type(polynomial_form) function create() result(this)

    integer :: m

    this % support = support(GRAPH_SIDE_VERTEX, [(m, m = 1, 4)])

  end function create

  pure integer function polynomial_size(this)

    class(polynomial_form), intent(in) :: this

    associate (u1 => this); end associate

    polynomial_size = 4

  end function polynomial_size

  pure subroutine polynomial_values(this, x, at, phi)

    class(polynomial_form), intent(in) :: this
    real(dp), intent(in)  :: x(3), at(3)
    real(dp), intent(out) :: phi(:)

    associate (u1 => this); end associate

    phi(1)   = 1.0_dp
    phi(2:4) = x - at

  end subroutine polynomial_values

  pure subroutine polynomial_slopes(this, x, at, direction, dphi)

    class(polynomial_form), intent(in) :: this
    real(dp), intent(in)  :: x(3), at(3), direction(3)
    real(dp), intent(out) :: dphi(:)

    associate (u1 => this, u2 => x, u3 => at); end associate

    dphi(1)   = 0.0_dp
    dphi(2:4) = direction

  end subroutine polynomial_slopes

end module structure_polynomial_form
