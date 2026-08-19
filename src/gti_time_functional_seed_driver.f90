!=====================================================================!
! GTI TIME FUNCTIONAL SEED DRIVER
!
! Functional-to-seed construction for the time graph. A time
! functional is a list of scalar terms on time vertices,
!
!      F_time = sum_k F_k(q_{v_k}, xi, t_k),
!
! one functional form serving every term in this phase, and the
! seed driver converts the terms into the three things the
! reverse pass consumes or reports:
!
!      value                sum_k F_k, the scalar itself
!      vertex q-seeds       vertex_seed(v) += F_{k,q}^T for every
!                           term on v - duplicate terms add
!      explicit design      sum_k F_{k,xi}[eta]
!
! The full design gradient is assembled OUTSIDE this seat as
!
!      explicit F_xi[eta]  +  reverse residual contribution,
!
! and this driver never calls the reverse pass - it builds seeds;
! whoever owns them walks the graph.
!
! Each term is evaluated on a point holding exactly the vertex's
! state, the design, and the term's evaluation time - a DIRECT
! vertex functional, no motif relation anywhere. The q-gradient
! is one order-1 state partial action per basis direction, each
! answer a scalar because F is one - never a finite difference.
!
! Seeds live outside the graph, as ever: the caller owns the
! vertex_seed array, empty meaning zero, and the graph's primal
! content is never touched.
!
! The functional representation carries term addresses and times,
! and nothing else: no form, no graph pointer, no solver, no
! seeds, no design. The driver carries nothing at all.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_functional_seed_drivers

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, &
       & GTI_ARG_STATE, GTI_ARG_DESIGN
  use gti_form_evaluators  , only : gti_form_evaluator
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_time_local_schemes, only : gti_evaluation_point
  use gti_time_graphs      , only : gti_time_graph

  implicit none

  private
  public :: gti_time_functional_term
  public :: gti_time_functional
  public :: gti_time_functional_seed_result
  public :: gti_time_functional_seed_driver

  !===================================================================!
  ! One scalar term's address: which vertex, at which time.
  !===================================================================!

  type :: gti_time_functional_term

     integer  :: vertex_index = 0
     real(dp) :: evaluation_time = 0.0_dp

  end type gti_time_functional_term

  !===================================================================!
  ! The time functional: a list of terms, and nothing else.
  !===================================================================!

  type :: gti_time_functional

     type(gti_time_functional_term), allocatable :: term(:)

   contains

     procedure :: num_terms
     procedure :: validate => functional_validate

  end type gti_time_functional

  !===================================================================!
  ! What one seeding pass reports: the accumulated scalar value
  ! and the explicit design action - the seeds themselves land in
  ! the caller's array.
  !===================================================================!

  type :: gti_time_functional_seed_result

     logical  :: completed = .false.
     real(dp) :: value = 0.0_dp
     type(gti_value_buffer) :: explicit_design_action

  end type gti_time_functional_seed_result

  !===================================================================!
  ! The stateless verb-set. The types keep their public singular
  ! names; Fortran denies a type its host module's name, so the
  ! module speaks in the plural.
  !===================================================================!

  type :: gti_time_functional_seed_driver

   contains

     procedure :: build_point
     procedure :: value         => term_value
     procedure :: q_gradient
     procedure :: design_action
     procedure :: seed_all

  end type gti_time_functional_seed_driver

contains

  !===================================================================!
  ! How many scalar terms the functional carries.
  !===================================================================!

  pure function num_terms(this) result(n)

    class(gti_time_functional), intent(in) :: this
    integer :: n

    if (allocated(this % term)) then
       n = size(this % term)
    else
       n = 0
    end if

  end function num_terms

  !===================================================================!
  ! The representation laws: terms exist, every term's vertex is
  ! in the graph, solved, and providing q.
  !===================================================================!

  subroutine functional_validate(this, graph)

    class(gti_time_functional), intent(in) :: this
    type(gti_time_graph)      , intent(in) :: graph

    integer :: k, vertex_index

    if (this % num_terms() < 1) then
       error stop 'gti_time_functional: at least one term is required'
    end if

    do k = 1, this % num_terms()

       vertex_index = this % term(k) % vertex_index
       if (vertex_index < 1 .or. vertex_index > graph % num_vertices()) then
          error stop 'gti_time_functional: term vertex is in range'
       end if

       if (.not. graph % vertex(vertex_index) % has_solution) then
          error stop 'gti_time_functional: term vertex has solution'
       end if

       if (.not. graph % vertex(vertex_index) % sample % state % &
            & has_component(GTI_STATE_Q)) then
          error stop 'gti_time_functional: term vertex provides q'
       end if

    end do

  end subroutine functional_validate

  !===================================================================!
  ! The point one term sees: exactly the vertex's state, the
  ! design, and the term's time. A direct vertex functional - no
  ! motif relation is constructed, and the graph is not touched.
  !===================================================================!

  subroutine build_point(this, graph, term, design, point)

    class(gti_time_functional_seed_driver), intent(in)    :: this
    type(gti_time_graph)                  , intent(in)    :: graph
    type(gti_time_functional_term)        , intent(in)    :: term
    type(gti_design_bundle)               , intent(in)    :: design
    type(gti_evaluation_point)            , intent(inout) :: point

    if (term % vertex_index < 1 .or. &
         & term % vertex_index > graph % num_vertices()) then
       error stop 'gti_time_functional: term vertex is in range'
    end if

    if (.not. graph % vertex(term % vertex_index) % sample % state % &
         & has_component(GTI_STATE_Q)) then
       error stop 'gti_time_functional: term vertex provides q'
    end if

    point % state     = graph % vertex(term % vertex_index) % sample % state
    point % design    = design
    point % time      = term % evaluation_time
    point % window_id = 0
    point % step_id   = 0
    point % stage_id  = 0

  end subroutine build_point

  !===================================================================!
  ! One term's scalar value, held to scalarness.
  !===================================================================!

  subroutine term_value(this, functional_form, graph, term, design, output)

    class(gti_time_functional_seed_driver), intent(in)    :: this
    class(gti_differentiable_form)        , intent(in)    :: functional_form
    type(gti_time_graph)                  , intent(in)    :: graph
    type(gti_time_functional_term)        , intent(in)    :: term
    type(gti_design_bundle)               , intent(in)    :: design
    type(gti_value_buffer)                , intent(inout) :: output

    type(gti_form_evaluator)   :: evaluator
    type(gti_evaluation_point) :: point

    call this % build_point(graph, term, design, point)

    call evaluator % value(functional_form, point, output)

    call require_scalar(output, &
         & 'gti_time_functional_seed_driver: functional value is scalar')

  end subroutine term_value

  !===================================================================!
  ! One term's q-gradient: the dense g with g^T dq = dF_k/dq_v,
  ! one order-1 state partial action per basis direction, each a
  ! scalar because F is one. Never a finite difference.
  !===================================================================!

  subroutine q_gradient(this, functional_form, graph, term, design, gradient)

    class(gti_time_functional_seed_driver), intent(in)    :: this
    class(gti_differentiable_form)        , intent(in)    :: functional_form
    type(gti_time_graph)                  , intent(in)    :: graph
    type(gti_time_functional_term)        , intent(in)    :: term
    type(gti_design_bundle)               , intent(in)    :: design
    type(gti_value_buffer)                , intent(inout) :: gradient

    type(gti_form_evaluator)   :: evaluator
    type(gti_evaluation_point) :: point
    type(gti_partial_request)  :: request
    type(gti_direction_bundle) :: direction
    type(gti_value_buffer)     :: term_buffer

    real(dp), allocatable :: basis(:), term_values(:), gradient_values(:)
    integer :: nvalues, ncomp, j

    call this % build_point(graph, term, design, point)

    call extract_vertex_q_shape(graph, term % vertex_index, nvalues, ncomp)

    allocate(basis(nvalues), gradient_values(nvalues))

    do j = 1, nvalues

       basis    = 0.0_dp
       basis(j) = 1.0_dp

       direction % argument_kind   = GTI_ARG_STATE
       direction % state_component = GTI_STATE_Q
       call direction % values % set_real(basis, ncomp=ncomp)

       request = gti_partial_request(order=1, &
            & argument_kind  =[GTI_ARG_STATE], &
            & state_component=[GTI_STATE_Q])

       call evaluator % partial_action(functional_form, point, request, &
            & [direction], term_buffer)

       call require_scalar(term_buffer, &
            & 'gti_time_functional_seed_driver: functional action is scalar')

       call term_buffer % get_real(term_values)
       gradient_values(j) = term_values(1)

    end do

    call gradient % set_real(gradient_values)

  end subroutine q_gradient

  !===================================================================!
  ! One term's explicit design action, F_xi[eta], held to
  ! scalarness.
  !===================================================================!

  subroutine design_action(this, functional_form, graph, term, design, &
       & design_direction, action)

    class(gti_time_functional_seed_driver), intent(in)    :: this
    class(gti_differentiable_form)        , intent(in)    :: functional_form
    type(gti_time_graph)                  , intent(in)    :: graph
    type(gti_time_functional_term)        , intent(in)    :: term
    type(gti_design_bundle)               , intent(in)    :: design
    type(gti_value_buffer)                , intent(in)    :: design_direction
    type(gti_value_buffer)                , intent(inout) :: action

    type(gti_form_evaluator)   :: evaluator
    type(gti_evaluation_point) :: point
    type(gti_partial_request)  :: request
    type(gti_direction_bundle) :: direction

    integer :: ndirection

    ndirection = 0
    if (allocated(design_direction % rvals)) ndirection = size(design_direction % rvals)
    if (ndirection == 0) then
       error stop 'gti_time_functional_seed_driver: design direction has values'
    end if

    call this % build_point(graph, term, design, point)

    direction % argument_kind = GTI_ARG_DESIGN
    direction % values        = design_direction

    request = gti_partial_request(order=1, argument_kind=[GTI_ARG_DESIGN])

    call evaluator % partial_action(functional_form, point, request, &
         & [direction], action)

    call require_scalar(action, &
         & 'gti_time_functional_seed_driver: functional action is scalar')

  end subroutine design_action

  !===================================================================!
  ! The whole seeding pass: every term's value gathered, every
  ! term's q-gradient landed on its vertex seed - duplicates add -
  ! and the explicit design actions summed into one scalar.
  !===================================================================!

  subroutine seed_all(this, functional_form, graph, functional, design, &
       & design_direction, vertex_seed, result)

    class(gti_time_functional_seed_driver), intent(in)    :: this
    class(gti_differentiable_form)        , intent(in)    :: functional_form
    type(gti_time_graph)                  , intent(in)    :: graph
    type(gti_time_functional)             , intent(in)    :: functional
    type(gti_design_bundle)               , intent(in)    :: design
    type(gti_value_buffer)                , intent(in)    :: design_direction
    type(gti_value_buffer)                , intent(inout) :: vertex_seed(:)
    type(gti_time_functional_seed_result) , intent(inout) :: result

    type(gti_value_buffer) :: value_buffer, gradient, action

    real(dp), allocatable :: scalar_values(:), gradient_values(:)
    real(dp) :: explicit_total
    integer :: k, ndirection

    call graph % validate()
    call functional % validate(graph)

    if (size(vertex_seed) /= graph % num_vertices()) then
       error stop 'gti_time_functional_seed_driver: vertex seed size matches graph vertices'
    end if

    ndirection = 0
    if (allocated(design_direction % rvals)) ndirection = size(design_direction % rvals)
    if (ndirection == 0) then
       error stop 'gti_time_functional_seed_driver: design direction has values'
    end if

    result % completed = .false.
    result % value     = 0.0_dp
    explicit_total     = 0.0_dp

    do k = 1, functional % num_terms()

       call this % value(functional_form, graph, functional % term(k), &
            & design, value_buffer)
       call value_buffer % get_real(scalar_values)
       result % value = result % value + scalar_values(1)

       call this % q_gradient(functional_form, graph, functional % term(k), &
            & design, gradient)
       call gradient % get_real(gradient_values)
       call add_to_seed(vertex_seed(functional % term(k) % vertex_index), &
            & gradient_values)

       call this % design_action(functional_form, graph, functional % term(k), &
            & design, design_direction, action)
       call action % get_real(scalar_values)
       explicit_total = explicit_total + scalar_values(1)

    end do

    result % completed = .true.
    call result % explicit_design_action % set_real([explicit_total])

  end subroutine seed_all

  !===================================================================!
  ! Scalarness: one entry, one component, one stored real.
  !===================================================================!

  subroutine require_scalar(buffer, law)

    type(gti_value_buffer), intent(in) :: buffer
    character(len=*)      , intent(in) :: law

    integer :: nstored

    nstored = 0
    if (allocated(buffer % rvals)) nstored = size(buffer % rvals)

    if (buffer % nentries /= 1 .or. buffer % ncomp /= 1 .or. nstored /= 1) then
       error stop law
    end if

  end subroutine require_scalar

  !===================================================================!
  ! The referenced vertex q's shape, with its presence proven:
  ! the seat holds a field, and the field holds reals.
  !===================================================================!

  subroutine extract_vertex_q_shape(graph, vertex_index, nvalues, ncomp)

    type(gti_time_graph), intent(in)  :: graph
    integer             , intent(in)  :: vertex_index
    integer             , intent(out) :: nvalues, ncomp

    real(dp), allocatable :: q_values(:)

    if (.not. graph % vertex(vertex_index) % sample % state % &
         & has_component(GTI_STATE_Q)) then
       error stop 'gti_time_functional: term vertex provides q'
    end if

    call graph % vertex(vertex_index) % sample % state % &
         & component(1 + GTI_STATE_Q) % value % get_real_vector(q_values)
    if (size(q_values) == 0) then
       error stop 'gti_time_functional_seed_driver: vertex q has values'
    end if

    nvalues = size(q_values)
    ncomp   = graph % vertex(vertex_index) % sample % state % &
         & component(1 + GTI_STATE_Q) % value % num_components()

  end subroutine extract_vertex_q_shape

  !===================================================================!
  ! Accumulate into one vertex seed: an empty seed is zero and
  ! takes the incoming shape; an occupied seed must agree in
  ! shape, and adds entrywise.
  !===================================================================!

  subroutine add_to_seed(target, values)

    type(gti_value_buffer), intent(inout) :: target
    real(dp)              , intent(in)    :: values(:)

    integer :: nstored

    nstored = 0
    if (allocated(target % rvals)) nstored = size(target % rvals)

    if (nstored == 0) then
       call target % set_real(values)
       return
    end if

    if (target % nentries /= size(values) .or. target % ncomp /= 1 .or. &
         & nstored /= size(values)) then
       error stop 'gti_time_functional_seed_driver: propagated seed shape matches target'
    end if

    target % rvals = target % rvals + values

  end subroutine add_to_seed

end module gti_time_functional_seed_drivers
