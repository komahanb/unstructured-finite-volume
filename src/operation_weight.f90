!=====================================================================!
! The weight of an edge: the product of a family's dimensionless
! coefficient with the step scaling that carries the source's units
! into its constraint's.
!
!      weight  =  alpha  *  dt_head ** (source_degree - determines)
!
! This is the product of the two fields the coupling carries - the
! tau field of step powers and the alpha field of coefficients - and
! it is the entry the constraint row holds for that source. The row
! itself reads
!
!      (the component this constraint determines)
!          -  sum over the edges into it of weight * (its source)   =  0
!
! so the determined component enters with one and every source enters
! with its weight. A backward difference for the velocity gives
! alpha_j / dt_k on q_(k-j); an Adams velocity row gives one on
! q'_(k-1) and dt_k alpha_i on q"_(k-i); a Runge-Kutta stage gives
! one on the incoming component and dt a_ij on the stage above it.
!
! A weight is itself a family: it declares the same march as the
! coefficients it multiplies, so it plugs in wherever they do and the
! caller chooses whether the step powers are already carried.
!
!             THE STEP AT A STAGE
!
! The step is read at the edge's head, so a coupling whose vertices
! are stages must carry the step of the step being taken at every one
! of those stage vertices, not only at the two instants. A stage
! vertex left at zero stops the program in any edge whose exponent is
! negative, and silently carries the wrong scale in any edge whose
! exponent is positive.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_weight

  use util_precision  , only : dp
  use operation_family      , only : family
  use util_derivative_terms , only : derivative_terms, integer_power, operator(*)

  implicit none

  private
  public :: scheme_weight

  type, extends(family) :: scheme_weight

     class(family), private, allocatable :: coefficients

   contains

     procedure :: name             => weight_name
     procedure :: history_depth    => weight_history_depth
     procedure :: num_stages       => weight_num_stages
     procedure :: primary_degree   => weight_primary_degree
     procedure :: row_pattern      => weight_row_pattern
     procedure :: edge_coefficient => weight_edge_coefficient

  end type scheme_weight

  interface scheme_weight
     module procedure create
  end interface scheme_weight

contains

  !===================================================================!
  ! The weights of a family. The family is copied in, so the weights
  ! outlive whatever was handed here.
  !===================================================================!

  function create(coefficients) result(this)

    class(family), intent(in) :: coefficients
    type(scheme_weight) :: this

    allocate(this % coefficients, source=coefficients)
    call this % declare_arguments(3)

  end function create

  pure function weight_name(this) result(name)

    class(scheme_weight), intent(in) :: this
    character(len=:), allocatable :: name

    name = this % coefficients % name() // ' weight'

  end function weight_name

  !===================================================================!
  ! The march is the one the coefficients declare: multiplying by a
  ! power of the step changes no reach, no stage count, no unknown
  ! and no sparsity.
  !===================================================================!

  pure integer function weight_history_depth(this, equation_degree)

    class(scheme_weight), intent(in) :: this
    integer             , intent(in) :: equation_degree

    weight_history_depth = this % coefficients % history_depth(equation_degree)

  end function weight_history_depth

  pure integer function weight_num_stages(this)

    class(scheme_weight), intent(in) :: this

    weight_num_stages = this % coefficients % num_stages()

  end function weight_num_stages

  pure integer function weight_primary_degree(this, equation_degree)

    class(scheme_weight), intent(in) :: this
    integer             , intent(in) :: equation_degree

    weight_primary_degree = this % coefficients % primary_degree(equation_degree)

  end function weight_primary_degree

  pure subroutine weight_row_pattern(this, determines, equation_degree, &
       & offset, source_degree)

    class(scheme_weight), intent(in) :: this
    integer             , intent(in) :: determines, equation_degree
    integer, allocatable, intent(out) :: offset(:), source_degree(:)

    call this % coefficients % row_pattern(determines, equation_degree, &
         & offset, source_degree)

  end subroutine weight_row_pattern

  !===================================================================!
  ! The product. Both factors carry their own partials in the steps,
  ! so the product rule is applied by the arithmetic and the result
  ! is exact at every degree.
  !===================================================================!

  pure function weight_edge_coefficient(this, dt, tail, head, &
       & source_degree, determines) result(c)

    class(scheme_weight)  , intent(in) :: this
    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: tail, head, source_degree, determines
    type(derivative_terms) :: c

    c = this % coefficients % edge_coefficient(dt, tail, head, &
         & source_degree, determines) &
         & * integer_power(dt(head), source_degree - determines)

  end function weight_edge_coefficient

end module operation_weight
