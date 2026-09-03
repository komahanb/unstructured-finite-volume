!=====================================================================!
! The weight of an edge: the product of a family's dimensionless
! coefficient with the step scaling that converts the source's units
! into its constraint's.
!
!      weight  =  alpha  *  dt_head ** (tail_degree - head_degree)
!
! This is the product of the two fields the coupling stores - the
! tau field of step powers and the alpha field of coefficients - and
! it is the entry the constraint row stores for that source. The row
! itself reads
!
!      (the component this constraint head_degree)
!          -  sum over the edges into it of weight * (its source)   =  0
!
! so the determined component enters with one and every source enters
! with its weight. A backward difference for the velocity gives
! alpha_j / dt_k on q_(k-j); an Adams velocity row gives one on
! q'_(k-1) and dt_k alpha_i on q"_(k-i); a Runge-Kutta stage gives
! one on the incoming component and dt a_ij on the stage above it.
!
! A weight is an edge function over the same coupling as the
! coefficients it multiplies, so it is substitutable wherever they
! are, and the caller chooses whether the step powers are already
! included.
!
!             THE STEP AT A STAGE
!
! The step is read at the edge's head, so a coupling whose vertices
! are stages must store the step size of the step being taken at
! every one of those stage vertices, not only at the two instants. A
! stage vertex left at zero stops the program in any edge whose
! exponent is negative, and produces the wrong scale without an
! error in any edge whose exponent is positive.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_weight

  use util_precision  , only : dp
  use operation_edge_function, only : edge_function
  use util_derivative_terms  , only : derivative_terms, integer_power, operator(*)

  implicit none

  private
  public :: scheme_weight

  type, extends(edge_function) :: scheme_weight

     class(edge_function), private, allocatable :: coefficients

   contains

     procedure :: edge_coefficient => weight_edge_coefficient

  end type scheme_weight

  interface scheme_weight
     module procedure create
  end interface scheme_weight

contains

  !===================================================================!
  ! The weights of a family, named after it. The family is copied,
  ! so the weights outlive the argument passed here.
  !===================================================================!

  function create(coefficients) result(this)

    class(edge_function), intent(in) :: coefficients
    type(scheme_weight) :: this

    allocate(this % coefficients, source=coefficients)
    call this % declare_edge_arguments(coefficients % name() // ' weight')

  end function create

  !===================================================================!
  ! The product. Both factors store their own partials in the steps,
  ! so the product rule is applied by the arithmetic and the result
  ! is exact at every degree.
  !===================================================================!

  pure function weight_edge_coefficient(this, dt, tail, head, &
       & tail_degree, head_degree) result(c)

    class(scheme_weight)  , intent(in) :: this
    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: tail, head, tail_degree, head_degree
    type(derivative_terms) :: c

    c = this % coefficients % edge_coefficient(dt, tail, head, &
         & tail_degree, head_degree) &
         & * integer_power(dt(head), tail_degree - head_degree)

  end function weight_edge_coefficient

end module operation_weight
