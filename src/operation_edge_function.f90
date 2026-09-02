!=====================================================================!
! A differentiable function on the edges of a coupling, determined by
! the step field.
!
!      input graph    the coupling: one vertex per instant, or per
!                     stage, and one edge from each vertex a vertex
!                     reads into that vertex
!      input 1        dt, a real field on the vertices: the step that
!                     ends at each vertex, read at an edge's head
!      input 2        the derivative degree of each edge's source
!      input 3        the degree each edge's constraint determines
!      output         one real per edge
!
! A concretion supplies the rule for one edge and nothing else. The
! edges are traversed, the inputs are read, and the result is placed
! on the coupling's edge set here once, for every concretion.
!
!             THE PARTIALS IN THE STEPS, TO ANY DEGREE
!
! The rule is evaluated over derivative_terms, so the value and every
! mixed partial in the steps are computed together and the result is
! the coefficient of the full subset: with no directions that is the
! value, with n directions it is the n-th mixed partial. apply and
! partial_action are therefore one code path, no perturbation is
! made, and the partials are exact.
!
!             WHAT IS REFUSED
!
! Missing inputs, and a variation on anything but the steps, stop the
! program. A step that a rule divides by must be positive, and is
! refused where it is read: the first vertex of a chain has no step
! ending at it and has step zero, which is only an error if a rule
! reads it. An edge a concretion does not define is refused by
! that concretion.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_edge_function

  use util_precision  , only : dp
  use operation_action      , only : operation, variation, contract
  use operation_action      , only : binding, bound_integer_vector, seeded_argument
  use operation_action, only : emit_real
  use view_directed         , only : directed_graph
  use field_calculus        , only : field, FIELD_REAL, FIELD_INTEGER
  use graph_fractal         , only : graph
  use util_derivative_terms , only : derivative_terms, mixed_partial, &
       & max_subset_width

  implicit none

  private
  public :: edge_function

  type, abstract, extends(operation) :: edge_function

   contains

     procedure(edge_coefficient_interface), deferred :: edge_coefficient

     procedure :: domain         => edge_domain
     procedure :: apply          => edge_apply
     procedure :: partial_action => edge_partial_action
     procedure :: declare_edge_arguments

  end type edge_function

  abstract interface

     !----------------------------------------------------------------!
     ! The value on one edge, from the step terms, the two vertices
     ! the edge joins, the derivative degree of its source and the
     ! degree its constraint determines.
     !----------------------------------------------------------------!

     pure function edge_coefficient_interface(this, dt, tail, head, &
          & source_degree, determines) result(c)
       import :: edge_function, derivative_terms
       class(edge_function)  , intent(in) :: this
       type(derivative_terms), intent(in) :: dt(:)
       integer               , intent(in) :: tail, head, source_degree, determines
       type(derivative_terms) :: c
     end function edge_coefficient_interface

  end interface

contains

  !
  ! The coefficients are defined on the coupling's edges.
  !===================================================================!

  subroutine edge_domain(this, input_graph, domain, num_entries)

    class(edge_function), intent(in)  :: this
    class(directed_graph), intent(in)  :: input_graph
    type(graph)          , intent(out) :: domain
    integer              , intent(out) :: num_entries

    associate (u1 => this); end associate
    domain      = input_graph % edge_set()
    num_entries = input_graph % num_edges()

  end subroutine edge_domain

  !===================================================================!
  ! The three arguments every edge function reads, the label it
  ! reports, and the highest exact degree: the width the subset
  ! masks can index.
  !===================================================================!

  subroutine declare_edge_arguments(this, label)

    class(edge_function), intent(in out) :: this
    character(len=*)    , intent(in)     :: label

    call this % declare_arguments(3, [contract(FIELD_REAL, 1), &
         & contract(FIELD_INTEGER, 1), contract(FIELD_INTEGER, 1)], &
         & label=label, max_degree=max_subset_width())

  end subroutine declare_edge_arguments

  !===================================================================!
  ! The coefficient of the full subset on every edge: the value when
  ! there are no directions, the mixed partial along all of them
  ! otherwise.
  !===================================================================!

  subroutine full_terms(this, input_graph, dt, source_degree, determines, output)

    class(edge_function), intent(in)  :: this
    class(directed_graph) , intent(in) :: input_graph
    type(derivative_terms), intent(in) :: dt(:)
    integer               , intent(in) :: source_degree(:), determines(:)
    class(field), allocatable, intent(inout) :: output

    type(derivative_terms) :: c
    real(dp), allocatable :: values(:)
    integer :: e

    allocate(values(input_graph % num_edges()))

    do e = 1, input_graph % num_edges()
       c = this % edge_coefficient(dt, input_graph % edge_tail(e), &
            & input_graph % edge_head(e), source_degree(e), determines(e))
       values(e) = mixed_partial(c)
    end do

    call emit_real(this % name() // ' coefficients', input_graph % edge_set(), &
         & input_graph % num_edges(), values, output)

  end subroutine full_terms

  subroutine edge_apply(this, input_graph, inputs, output)

    class(edge_function), intent(in)  :: this
    class(directed_graph), intent(in)        :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    call this % value_by_partial_action(input_graph, inputs, output)

  end subroutine edge_apply

  !===================================================================!
  ! The steps as derivative terms seeded by the variations, which
  ! must all be on the steps, the first argument; one naming another
  ! argument stops the program.
  !===================================================================!

  subroutine edge_partial_action(this, input_graph, inputs, variations, output)

    class(edge_function), intent(in)  :: this
    class(directed_graph), intent(in)        :: input_graph
    type(binding)         , intent(in)        :: inputs(:)
    type(variation)      , intent(in)        :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(derivative_terms), allocatable :: dt(:)
    integer , allocatable :: source_degree(:), determines(:)
    integer :: consumed

    call this % require_variations(variations)
    call seeded_argument(this, inputs, variations, 1, dt, consumed)
    if (consumed < size(variations)) then
       error stop 'operation_edge_function: the coefficients vary with the steps alone'
    end if
    call bound_integer_vector(inputs, this % argument(2), source_degree)
    call bound_integer_vector(inputs, this % argument(3), determines)
    call full_terms(this, input_graph, dt, source_degree, determines, output)

  end subroutine edge_partial_action

end module operation_edge_function
