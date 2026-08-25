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
! mixed partial in the steps are carried together and the result is
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
! ending at it and carries zero, which is only an error if a rule
! asks for it. An edge a concretion does not define is refused by
! that concretion.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_edge_function

  use util_precision  , only : dp
  use operation_action      , only : operation, variation
  use view_directed         , only : directed_graph
  use field_calculus        , only : field
  use graph_fractal         , only : graph
  use field_stored          , only : stored_field
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
     procedure :: max_degree     => edge_max_degree
     procedure :: partial_action => edge_partial_action

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
  ! The coefficients live on the coupling's edges.
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
  ! The degree the subset masks can index: the bit width less the
  ! sign bit and the bit the full mask would carry past it.
  !===================================================================!

  pure integer function edge_max_degree(this)

    class(edge_function), intent(in) :: this

    associate (u1 => this); end associate
    edge_max_degree = max_subset_width()

  end function edge_max_degree

  !===================================================================!
  ! The three inputs, which must all be given.
  !===================================================================!

  subroutine read_inputs(input_data, dt, source_degree, determines)

    class(field), intent(in) :: input_data(:)
    real(dp), allocatable, intent(out) :: dt(:)
    integer , allocatable, intent(out) :: source_degree(:), determines(:)

    if (size(input_data) < 3) then
       error stop 'operation_edge_function: the steps, source degrees and conditions are given'
    end if

    call input_data(1) % real_vector(dt)
    call input_data(2) % integer_vector(source_degree)
    call input_data(3) % integer_vector(determines)

  end subroutine read_inputs

  !===================================================================!
  ! The steps as derivative terms: the value at each vertex, and along
  ! direction i the i-th variation's entry there. Every variation
  ! must be on the steps, the operation's first argument.
  !===================================================================!

  subroutine step_terms(this, dt, variations, terms)

    class(edge_function), intent(in)  :: this
    real(dp)       , intent(in) :: dt(:)
    type(variation), intent(in) :: variations(:)
    type(derivative_terms), allocatable, intent(out) :: terms(:)

    real(dp), allocatable :: v(:)
    integer :: n, i, k

    n = size(variations)
    allocate(terms(size(dt)))

    do k = 1, size(dt)
       terms(k) = derivative_terms(dt(k), n)
    end do

    do i = 1, n
       if (.not. variations(i) % argument_is(this % argument(1))) then
          error stop 'operation_edge_function: the coefficients vary with the steps alone'
       end if
       call variations(i) % direction(v)
       do k = 1, size(dt)
          call terms(k) % set_direction(i, v(k))
       end do
    end do

  end subroutine step_terms

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

    type(stored_field) :: out
    type(derivative_terms) :: c
    real(dp), allocatable :: values(:)
    integer :: e

    allocate(values(input_graph % num_edges()))

    do e = 1, input_graph % num_edges()
       c = this % edge_coefficient(dt, input_graph % edge_tail(e), &
            & input_graph % edge_head(e), source_degree(e), determines(e))
       values(e) = mixed_partial(c)
    end do

    out = stored_field(this % name() // ' coefficients', input_graph % edge_set(), &
         & input_graph % num_edges())
    call out % set_real_vector(values)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine full_terms

  subroutine edge_apply(this, input_graph, input_data, output)

    class(edge_function), intent(in)  :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(variation), allocatable :: none(:)
    type(derivative_terms), allocatable :: dt(:)
    real(dp), allocatable :: steps(:)
    integer , allocatable :: source_degree(:), determines(:)

    if (.not. present(input_data)) then
       error stop 'operation_edge_function: the steps, source degrees and conditions are given'
    end if

    call read_inputs(input_data, steps, source_degree, determines)
    allocate(none(0))
    call step_terms(this, steps, none, dt)
    call full_terms(this, input_graph, dt, source_degree, determines, output)

  end subroutine edge_apply

  subroutine edge_partial_action(this, input_graph, input_data, variations, output)

    class(edge_function), intent(in)  :: this
    class(directed_graph), intent(in)        :: input_graph
    class(field)         , intent(in)        :: input_data(:)
    type(variation)      , intent(in)        :: variations(:)
    class(field), allocatable, intent(inout) :: output

    type(derivative_terms), allocatable :: dt(:)
    real(dp), allocatable :: steps(:)
    integer , allocatable :: source_degree(:), determines(:)

    call this % require_owned(variations)

    if (size(variations) > this % max_degree()) then
       error stop 'operation_edge_function: the requested order is within max_degree'
    end if

    call read_inputs(input_data, steps, source_degree, determines)
    call step_terms(this, steps, variations, dt)
    call full_terms(this, input_graph, dt, source_degree, determines, output)

  end subroutine edge_partial_action
end module operation_edge_function
