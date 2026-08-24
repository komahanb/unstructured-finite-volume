!=====================================================================!
! The step scaling of an edge: the power of the step that carries a
! source's units into its constraint's.
!
!      scaling  =  dt_head ** (source_degree - determines)
!
! A constraint that determines a derivative of degree d and reads a
! source of degree d' is a statement about d' - d differentiations,
! and each of those carries one power of the step. The exponent is
! therefore fixed by the two degrees the edge joins and by nothing
! else - not by the family, not by the order, not by the position in
! the history. A backward difference for the velocity reads values
! into a first derivative and scales by dt^-1; the same family's
! acceleration row scales by dt^-2; an Adams quadrature reads an
! acceleration into a velocity and scales by dt^+1; a source of the
! constraint's own degree scales by dt^0 and is carried unchanged.
!
! This is the tau side of the product that forms a scheme's weights.
! The alpha side is the family's dimensionless coefficient, and
! operation_weight multiplies the two.
!
!             WHAT IS REFUSED
!
! A zero or negative step at an edge's head, where the exponent is
! negative: the reciprocal is not defined and the quotient stops the
! program. A zero exponent never reads the step, so a vertex with no
! step ending at it is refused only by an edge that actually divides
! by it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_step_scaling

  use iso_fortran_env         , only : dp => REAL64
  use operation_edge_function , only : edge_function
  use util_derivative_terms   , only : derivative_terms, integer_power

  implicit none

  private
  public :: step_scaling

  type, extends(edge_function) :: step_scaling

   contains

     procedure :: name             => scaling_name
     procedure :: edge_coefficient => scaling_edge_coefficient

  end type step_scaling

  interface step_scaling
     module procedure create
  end interface step_scaling

contains

  function create() result(this)

    type(step_scaling) :: this

    call this % declare_arguments(3)

  end function create

  pure function scaling_name(this) result(name)

    class(step_scaling), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate
    name = 'step scaling'

  end function scaling_name

  pure function scaling_edge_coefficient(this, dt, tail, head, &
       & source_degree, determines) result(c)

    class(step_scaling)   , intent(in) :: this
    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: tail, head, source_degree, determines
    type(derivative_terms) :: c

    associate (u1 => this, u2 => tail); end associate

    c = integer_power(dt(head), source_degree - determines)

  end function scaling_edge_coefficient

end module operation_step_scaling
