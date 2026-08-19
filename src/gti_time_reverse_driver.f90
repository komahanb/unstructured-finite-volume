!=====================================================================!
! GTI TIME REVERSE SEED DRIVER
!
! The reverse traversal over the explicit time graph, from
! supplied vertex adjoint seeds:
!
!      for r = num_relations, ..., 1:
!         solve   J_u^T lambda_r = seed_u
!         gather  g_xi[eta] += - lambda_r^T R_xi[eta]
!         push    seed_h    += - R_h^T lambda_r   for each history h
!
! One relation's constraint is R_r(q_u, q_h, xi, t) = 0; forward
! solved q_u from q_h, and reverse runs the same coupling
! backwards: the seed on the unknown - the covector dF/dq_u the
! caller supplies or an earlier reverse step accumulated - defines
! the relation adjoint through the TRANSPOSE of the very Jacobian
! Newton used, the design gradient collects one contraction per
! relation, and the history vertices inherit their seeds through
! the motif-row actions at THEIR local weights:
!
!      R_h[dq] = sum_m D_m R [ w_m(h) dq ],
!
! the same chain rule as J_u, at a different column of the motif.
! No derivative is ever finite-differenced, and no scheme is named.
!
! Seeds live OUTSIDE the graph: the caller hands a vertex_seed
! array, one buffer per vertex, empty meaning zero - the graph
! stays primal-only, storing q and time structure, while reverse
! accumulation is owned externally by whoever asked for it.
!
! This phase starts from supplied seeds; the functional-to-seed
! construction is a later seat. Reverse stored order is the
! execution order: no topological scheduler, no adaptivity, no
! traversal abstraction, and the forward driver is never named.
!
! The driver carries nothing: no form, no graph, no seeds, no
! solver state, no design, no map.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_reverse_drivers

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, GTI_ARG_STATE
  use gti_form_evaluators  , only : gti_form_evaluator
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif, &
       & gti_time_local_residual_evaluator, gti_evaluation_point
  use gti_time_local_newton_drivers , only : gti_time_local_newton_driver
  use gti_time_local_adjoint_drivers, only : gti_time_local_adjoint_driver
  use gti_time_graphs      , only : gti_time_graph
  ! the dense Jacobian rides class_graph_stencil as a weighted-edge
  ! graph inside this adapter, and dense_direct minimizes it there
  use class_graph_dense_direct, only : solve_dense_matrix_with_dense_direct

  implicit none

  private
  public :: gti_time_reverse_driver
  public :: gti_time_reverse_options
  public :: gti_time_reverse_step_result
  public :: gti_time_reverse_result

  !===================================================================!
  ! The one knob of the reverse pass.
  !===================================================================!

  type :: gti_time_reverse_options

     real(dp) :: singular_tolerance = 1.0e-14_dp

   contains

     procedure :: validate => options_validate

  end type gti_time_reverse_options

  !===================================================================!
  ! What one relation's adjoint reports: which relation, which
  ! vertex seeded it, the adjoint and the seed it answered, the
  ! design sensitivity it contracted against, and its scalar
  ! contribution.
  !===================================================================!

  type :: gti_time_reverse_step_result

     integer  :: relation_index = 0
     integer  :: unknown_vertex = 0
     real(dp) :: linear_residual_norm = huge(1.0_dp)
     type(gti_value_buffer) :: lambda
     type(gti_value_buffer) :: unknown_seed
     type(gti_value_buffer) :: residual_design_action
     type(gti_value_buffer) :: design_gradient_contribution

  end type gti_time_reverse_step_result

  !===================================================================!
  ! What the whole reverse pass reports: completion, every step's
  ! account, and the accumulated design-gradient action.
  !===================================================================!

  type :: gti_time_reverse_result

     logical :: completed = .false.
     integer :: completed_relations = 0
     type(gti_time_reverse_step_result), allocatable :: step(:)
     type(gti_value_buffer) :: design_gradient_action

  end type gti_time_reverse_result

  !===================================================================!
  ! The stateless verb-set. The types keep their public singular
  ! names; Fortran denies a type its host module's name, so the
  ! module speaks in the plural.
  !===================================================================!

  type :: gti_time_reverse_driver

   contains

     procedure :: solve_relation_adjoint
     procedure :: propagate_relation
     procedure :: reverse_all

  end type gti_time_reverse_driver

contains

  pure subroutine options_validate(this)

    class(gti_time_reverse_options), intent(in) :: this

    if (this % singular_tolerance <= 0.0_dp) then
       error stop 'gti_time_reverse_driver: singular tolerance is positive'
    end if

  end subroutine options_validate

  !===================================================================!
  ! One relation's adjoint: J_u^T lambda = seed_u, the design
  ! sensitivity R_xi[eta], and the scalar contribution
  ! -lambda^T R_xi[eta]. The design action is taken FIRST, so a
  ! misshapen residual sensitivity dies with its own law before
  ! any linear algebra runs. Nothing here mutates the seeds.
  !===================================================================!

  subroutine solve_relation_adjoint(this, residual_form, graph, relation_index, &
       & design, design_direction, vertex_seed, options, result)

    class(gti_time_reverse_driver)    , intent(in)    :: this
    class(gti_differentiable_form)    , intent(in)    :: residual_form
    type(gti_time_graph)              , intent(in)    :: graph
    integer                           , intent(in)    :: relation_index
    type(gti_design_bundle)           , intent(in)    :: design
    type(gti_value_buffer)            , intent(in)    :: design_direction
    type(gti_value_buffer)            , intent(in)    :: vertex_seed(:)
    type(gti_time_reverse_options)    , intent(in)    :: options
    type(gti_time_reverse_step_result), intent(inout) :: result

    type(gti_time_local_newton_driver)  :: newton
    type(gti_time_local_adjoint_driver) :: adjoint
    type(gti_time_sample), allocatable  :: samples(:)
    type(gti_value_buffer)              :: q_star, r_action
    real(dp) :: achieved

    real(dp), allocatable :: q_values(:), seed_values(:), r_values(:)
    real(dp), allocatable :: jacobian(:,:), transposed(:,:), lambda_values(:)
    integer :: n, unknown, unknown_vertex, ndirection

    call options % validate()
    call graph % validate()

    if (relation_index < 1 .or. relation_index > graph % num_relations()) then
       error stop 'gti_time_reverse_driver: relation index is in range'
    end if

    if (size(vertex_seed) /= graph % num_vertices()) then
       error stop 'gti_time_reverse_driver: vertex seed size matches graph vertices'
    end if

    ndirection = 0
    if (allocated(design_direction % rvals)) ndirection = size(design_direction % rvals)
    if (ndirection == 0) then
       error stop 'gti_time_reverse_driver: design direction has values'
    end if

    unknown        = graph % relation(relation_index) % unknown_sample
    unknown_vertex = graph % relation(relation_index) % unknown_vertex()

    if (.not. graph % vertex(unknown_vertex) % has_solution) then
       error stop 'gti_time_reverse_driver: unknown vertex has solution'
    end if

    call vertex_seed(unknown_vertex) % get_real(seed_values)
    if (size(seed_values) == 0) then
       error stop 'gti_time_reverse_driver: unknown seed has values'
    end if

    !----------------------------------------------------------------!
    ! Materialize the solved local state.
    !----------------------------------------------------------------!

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the materialization replaces it wholesale
    allocate(samples(0))

    call graph % build_samples(relation_index, samples)

    call samples(unknown) % state % component(1 + GTI_STATE_Q) % value % &
         & get_real_vector(q_values)
    n = size(q_values)
    call q_star % set_real(q_values, ncomp=samples(unknown) % state % &
         & component(1 + GTI_STATE_Q) % value % num_components())

    if (size(seed_values) /= n) then
       error stop 'gti_time_reverse_driver: unknown seed size matches unknown size'
    end if

    !----------------------------------------------------------------!
    ! The design sensitivity first: a misshapen action dies with
    ! its own law before any linear algebra runs.
    !----------------------------------------------------------------!

    call adjoint % residual_design_action(residual_form, &
         & graph % relation(relation_index) % motif, samples, unknown, &
         & q_star, design, design_direction, &
         & graph % relation(relation_index) % evaluation_time, r_action)
    call r_action % get_real(r_values)
    if (size(r_values) /= n) then
       error stop 'gti_time_reverse_driver: residual design action size matches unknown size'
    end if

    !----------------------------------------------------------------!
    ! The dense J_u from Newton's exact action, transposed and
    ! eliminated against the seed.
    !----------------------------------------------------------------!

    call newton % dense_jacobian(residual_form, &
         & graph % relation(relation_index) % motif, samples, unknown, &
         & q_star, design, &
         & graph % relation(relation_index) % evaluation_time, jacobian)

    transposed = transpose(jacobian)

    call solve_dense_matrix_with_dense_direct(transposed, seed_values, &
         & options % singular_tolerance, lambda_values, achieved)

    result % relation_index       = relation_index
    result % unknown_vertex       = unknown_vertex
    result % linear_residual_norm = &
         & norm2(matmul(transposed, lambda_values) - seed_values)

    call result % lambda % set_real(lambda_values, ncomp=q_star % ncomp)
    result % unknown_seed           = vertex_seed(unknown_vertex)
    result % residual_design_action = r_action
    call result % design_gradient_contribution % set_real( &
         & [-dot_product(lambda_values, r_values)])

  end subroutine solve_relation_adjoint

  !===================================================================!
  ! Push one relation's adjoint to its history: for each
  ! non-unknown sample position h, the covector h_i with
  !
  !      h_i(j) = lambda^T R_h[e_j],
  !      R_h[dq] = sum_m D_m R [ w_m(h) dq ],
  !
  ! lands on the history vertex as  seed_h += - h_i. The unknown's
  ! own seed is not touched here.
  !===================================================================!

  subroutine propagate_relation(this, residual_form, graph, relation_index, &
       & design, lambda, vertex_seed)

    class(gti_time_reverse_driver), intent(in)    :: this
    class(gti_differentiable_form), intent(in)    :: residual_form
    type(gti_time_graph)          , intent(in)    :: graph
    integer                       , intent(in)    :: relation_index
    type(gti_design_bundle)       , intent(in)    :: design
    type(gti_value_buffer)        , intent(in)    :: lambda
    type(gti_value_buffer)        , intent(inout) :: vertex_seed(:)

    type(gti_time_local_residual_evaluator) :: point_builder
    type(gti_time_sample), allocatable      :: samples(:)
    type(gti_evaluation_point)              :: point

    real(dp), allocatable :: lambda_values(:), action_values(:), basis(:), h(:)
    integer :: i, j, ndq, unknown, vertex_index

    call graph % validate()

    if (relation_index < 1 .or. relation_index > graph % num_relations()) then
       error stop 'gti_time_reverse_driver: relation index is in range'
    end if

    if (size(vertex_seed) /= graph % num_vertices()) then
       error stop 'gti_time_reverse_driver: vertex seed size matches graph vertices'
    end if

    call lambda % get_real(lambda_values)
    if (size(lambda_values) == 0) then
       error stop 'gti_time_reverse_driver: lambda has values'
    end if

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the materialization replaces it wholesale
    allocate(samples(0))

    call graph % build_samples(relation_index, samples)
    call point_builder % build_point(graph % relation(relation_index) % motif, &
         & samples, design, graph % relation(relation_index) % evaluation_time, &
         & point)

    unknown = graph % relation(relation_index) % unknown_sample

    do i = 1, graph % relation(relation_index) % arity()

       if (i == unknown) cycle
       vertex_index = graph % relation(relation_index) % sample_vertex(i)

       ndq = samples(i) % state % component(1 + GTI_STATE_Q) % value % &
            & num_entries() * samples(i) % state % component(1 + GTI_STATE_Q) % &
            & value % num_components()

       allocate(basis(ndq), h(ndq))

       do j = 1, ndq
          basis    = 0.0_dp
          basis(j) = 1.0_dp
          call residual_sample_action(residual_form, point, &
               & graph % relation(relation_index) % motif, i, basis, &
               & action_values)
          h(j) = dot_product(lambda_values, action_values)
       end do

       call add_to_seed(vertex_seed(vertex_index), -h)

       deallocate(basis, h)

    end do

  end subroutine propagate_relation

  !===================================================================!
  ! The whole reverse pass: relations in reverse stored order,
  ! each solved against its seed, its design contribution
  ! gathered, its history seeded - and the total design-gradient
  ! action handed back as one scalar.
  !===================================================================!

  subroutine reverse_all(this, residual_form, graph, design, design_direction, &
       & vertex_seed, options, result)

    class(gti_time_reverse_driver), intent(in)    :: this
    class(gti_differentiable_form), intent(in)    :: residual_form
    type(gti_time_graph)          , intent(in)    :: graph
    type(gti_design_bundle)       , intent(in)    :: design
    type(gti_value_buffer)        , intent(in)    :: design_direction
    type(gti_value_buffer)        , intent(inout) :: vertex_seed(:)
    type(gti_time_reverse_options), intent(in)    :: options
    type(gti_time_reverse_result) , intent(inout) :: result

    real(dp), allocatable :: contribution(:)
    real(dp) :: total
    integer :: r

    call options % validate()
    call graph % validate()

    if (graph % num_relations() < 1) then
       error stop 'gti_time_reverse_driver: at least one relation is required'
    end if

    if (size(vertex_seed) /= graph % num_vertices()) then
       error stop 'gti_time_reverse_driver: vertex seed size matches graph vertices'
    end if

    result % completed = .false.
    result % completed_relations = 0

    if (allocated(result % step)) deallocate(result % step)
    allocate(result % step(graph % num_relations()))

    total = 0.0_dp

    do r = graph % num_relations(), 1, -1

       call this % solve_relation_adjoint(residual_form, graph, r, design, &
            & design_direction, vertex_seed, options, result % step(r))

       call result % step(r) % design_gradient_contribution % get_real(contribution)
       total = total + contribution(1)

       call this % propagate_relation(residual_form, graph, r, design, &
            & result % step(r) % lambda, vertex_seed)

       result % completed_relations = result % completed_relations + 1

    end do

    result % completed = .true.
    call result % design_gradient_action % set_real([total])

  end subroutine reverse_all

  !===================================================================!
  ! R_i[dq] for one relation-local sample position: the motif-row
  ! chain rule at THAT position's weights,
  !
  !      R_i[dq] = sum_m D_m R [ w_m(i) dq ],
  !
  ! one order-1 state partial action per live row, summed. The
  ! same law as the unknown Jacobian, at a different column of the
  ! motif; never a finite difference.
  !===================================================================!

  subroutine residual_sample_action(residual_form, point, motif, &
       & local_position, dq_values, action_values)

    class(gti_differentiable_form), intent(in)  :: residual_form
    type(gti_evaluation_point)    , intent(in)  :: point
    type(gti_time_motif)          , intent(in)  :: motif
    integer                       , intent(in)  :: local_position
    real(dp)                      , intent(in)  :: dq_values(:)
    real(dp), allocatable         , intent(out) :: action_values(:)

    type(gti_form_evaluator)   :: evaluator
    type(gti_partial_request)  :: request
    type(gti_direction_bundle) :: direction
    type(gti_value_buffer)     :: term

    real(dp), allocatable :: term_values(:)
    real(dp) :: weight
    integer :: m
    logical :: started

    started = .false.

    do m = 1, motif % size()

       weight = motif % rule(m) % weight(local_position)
       if (weight == 0.0_dp) cycle

       direction % argument_kind   = GTI_ARG_STATE
       direction % state_component = motif % rule(m) % state_component
       call direction % values % set_real(weight * dq_values)

       request = gti_partial_request(order=1, &
            & argument_kind  =[GTI_ARG_STATE], &
            & state_component=[motif % rule(m) % state_component])

       call evaluator % partial_action(residual_form, point, request, &
            & [direction], term)

       call term % get_real(term_values)

       if (started) then
          action_values = action_values + term_values
       else
          action_values = term_values
          started       = .true.
       end if

    end do

    if (.not. started) then
       allocate(action_values(0))
    end if

  end subroutine residual_sample_action

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
       error stop 'gti_time_reverse_driver: propagated seed shape matches target'
    end if

    target % rvals = target % rvals + values

  end subroutine add_to_seed

end module gti_time_reverse_drivers
