!=====================================================================!
! GTI TIME DEGREE-R DIRECTIONAL TRAVERSAL
!
! The general forward derivative traversal over the time graph:
! given a solved G_time, one design direction eta, and a positive
! maximum degree r - tested through degree 8 - walk relations in
! stored forward order and compute every unknown vertex's
! directional derivatives
!
!      J_u q_u^(s) = -B^(s),        s = 1, ..., r,
!
! with ONE J_u per relation - built once from the local Jacobian
! action, before the degree loop, and eliminated r times. A higher
! degree is not a new driver; it is one more turn of the same
! loop, which is the architecture's central claim made executable:
! the Jacobian never changes with the degree, only the chain-rule
! right-hand side does, and the assembler generates every pattern
! the degree demands.
!
! The default path is affine: xi(eps) = xi + eps eta,
!
!      xi^(1) = eta,        xi^(s) = 0   for s >= 2,
!
! taken whenever the caller supplies no design path. A caller may
! instead occupy the design derivative seats directly, xi^(1),
! ..., xi^(r), through an optional design path array - and,
! independently, occupy the time derivative seats t^(1), ...,
! t^(r) through an optional time path array. An absent seat, in
! either path, means zero: the legacy affine path is the
! design-path-absent specialization of the same law, not a second
! mechanism living beside it.
!
! the total degree-s derivative of R(U(q(eps)), xi(eps), t(eps)) = 0
! splits into the unknown's transport term R_u[q_u^(s)] and
! everything else; everything else is B^(s), and it is assembled
! by gti_chain_rule_assemblies - never re-derived here. Each motif
! row m carries the component path
!
!      U_m^(k) = sum_i w_m(i) q_i^(k),        k = 1, ..., s,
!
! its samples' derivatives read from the caller's vertex array -
! an empty seat meaning zero - under the suppression rule that
! defines B^(s):
!
!      the unknown's q^(s) is zero while B^(s) is assembled,
!      the unknown's q^(k), k < s, is the derivative just solved,
!      history q^(k), k <= s, is read from the vertex array,
!      design/time seats are caller-provided or absent, and an
!      absent seat contributes nothing.
!
! The traversal reads the graph and writes only the caller's
! derivative array: no q value and no solution flag is touched,
! no reverse pass, no functional, no controller, and no scheduler
! exists here - relations run in stored order alone.
!
! The driver carries nothing: no graph, no forms, no solver
! state, no design, no derivatives, no map.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_degree_r_directional_drivers

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form, &
       & GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_ARG_TIME
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif, &
       & gti_time_local_residual_evaluator, gti_evaluation_point
  use gti_time_local_newton_drivers, only : gti_time_local_newton_driver
  use gti_time_graphs      , only : gti_time_graph
  use gti_chain_rule_assemblies, only : gti_chain_rule_assembler, &
       & gti_chain_channel, gti_chain_input

  implicit none

  private
  public :: gti_time_degree_r_directional_options
  public :: gti_time_degree_r_relation_result
  public :: gti_time_degree_r_directional_result
  public :: gti_time_degree_r_directional_driver

  !===================================================================!
  ! The two knobs of the degree-r traversal: how far the degree
  ! loop turns, and where a dense pivot stops being trusted.
  !===================================================================!

  type :: gti_time_degree_r_directional_options

     integer  :: max_degree = 1
     real(dp) :: singular_tolerance = 1.0e-14_dp

   contains

     procedure :: validate => options_validate

  end type gti_time_degree_r_directional_options

  !===================================================================!
  ! What one relation's r solves report: per degree, the
  ! right-hand side -B^(s), the derivative q_u^(s), and how
  ! exactly each linear system was satisfied.
  !===================================================================!

  type :: gti_time_degree_r_relation_result

     integer :: relation_index = 0
     integer :: unknown_vertex = 0
     integer :: max_degree = 0

     real(dp), allocatable :: linear_residual_norm(:)

     type(gti_value_buffer), allocatable :: rhs(:)
     type(gti_value_buffer), allocatable :: derivative(:)

  end type gti_time_degree_r_relation_result

  !===================================================================!
  ! What the whole traversal reports: completion, every step's
  ! account, and the per-vertex derivative array indexed
  ! (degree, vertex) - degree 1 stores no seat above itself.
  !===================================================================!

  type :: gti_time_degree_r_directional_result

     logical :: completed = .false.
     integer :: completed_relations = 0
     integer :: max_degree = 0

     type(gti_value_buffer), allocatable :: vertex_derivative(:,:)

     type(gti_time_degree_r_relation_result), allocatable :: step(:)

  end type gti_time_degree_r_directional_result

  !===================================================================!
  ! The stateless verb-pair. The types keep their public singular
  ! names; Fortran denies a type its host module's name, so the
  ! module speaks in the plural.
  !===================================================================!

  type :: gti_time_degree_r_directional_driver

   contains

     procedure :: solve_relation
     procedure :: solve_all

  end type gti_time_degree_r_directional_driver

contains

  pure subroutine options_validate(this)

    class(gti_time_degree_r_directional_options), intent(in) :: this

    if (this % max_degree < 1) then
       error stop 'gti_time_degree_r_directional_driver: degree is supported'
    end if

    if (this % singular_tolerance <= 0.0_dp) then
       error stop 'gti_time_degree_r_directional_driver: singular tolerance is positive'
    end if

  end subroutine options_validate

  !===================================================================!
  ! One relation's degrees 1 through r: the same J_u, built once
  ! before the degree loop, eliminated against each chain-rule
  ! right-hand side in turn.
  !===================================================================!

  subroutine solve_relation(this, residual_form, graph, relation_index, &
       & design, design_direction, vertex_derivative, options, result, &
       & design_path, time_path)

    class(gti_time_degree_r_directional_driver)  , intent(in)    :: this
    class(gti_differentiable_form)               , intent(in)    :: residual_form
    type(gti_time_graph)                         , intent(in)    :: graph
    integer                                      , intent(in)    :: relation_index
    type(gti_design_bundle)                      , intent(in)    :: design
    type(gti_value_buffer)                       , intent(in)    :: design_direction
    type(gti_value_buffer)                       , intent(inout) :: vertex_derivative(:,:)
    type(gti_time_degree_r_directional_options)  , intent(in)    :: options
    type(gti_time_degree_r_relation_result)      , intent(inout) :: result
    type(gti_value_buffer)                       , intent(in), optional :: design_path(:)
    type(gti_value_buffer)                       , intent(in), optional :: time_path(:)

    type(gti_time_local_residual_evaluator) :: point_builder
    type(gti_time_sample), allocatable      :: samples(:)
    type(gti_evaluation_point)              :: point
    type(gti_value_buffer)                  :: q_star, b_buffer

    real(dp), allocatable :: q_values(:)
    real(dp), allocatable :: local_path(:,:,:)
    real(dp), allocatable :: jacobian(:,:)
    real(dp), allocatable :: b_values(:), rhs(:), solution(:)
    integer , allocatable :: signature(:)
    integer :: n, ncomp, arity, unknown, unknown_vertex
    integer :: i, k, s, vertex_index

    call options % validate()
    call graph % validate()

    if (relation_index < 1 .or. relation_index > graph % num_relations()) then
       error stop 'gti_time_degree_r_directional_driver: relation index is in range'
    end if

    if (size(vertex_derivative, 1) /= options % max_degree) then
       error stop 'gti_time_degree_r_directional_driver: derivative array degree matches options'
    end if

    if (size(vertex_derivative, 2) /= graph % num_vertices()) then
       error stop 'gti_time_degree_r_directional_driver: derivative array vertex count matches graph'
    end if

    !----------------------------------------------------------------!
    ! The legacy affine path is the design-path-absent case: it
    ! alone requires design_direction to carry values. A caller
    ! that supplies design_path owns the design seats directly, and
    ! may pass an empty legacy buffer.
    !----------------------------------------------------------------!

    if (.not. present(design_path)) then
       k = 0
       if (allocated(design_direction % rvals)) k = size(design_direction % rvals)
       if (k == 0) then
          error stop 'gti_time_degree_r_directional_driver: design direction has values'
       end if
    end if

    arity          = graph % relation(relation_index) % arity()
    unknown        = graph % relation(relation_index) % unknown_sample
    unknown_vertex = graph % relation(relation_index) % unknown_vertex()

    if (.not. graph % vertex(unknown_vertex) % has_solution) then
       error stop 'gti_time_degree_r_directional_driver: unknown vertex has solution'
    end if

    do i = 1, arity
       if (i == unknown) cycle
       vertex_index = graph % relation(relation_index) % sample_vertex(i)
       if (.not. graph % vertex(vertex_index) % has_solution) then
          error stop 'gti_time_degree_r_directional_driver: history vertex has solution'
       end if
    end do

    !----------------------------------------------------------------!
    ! Materialize the solved local state; the size law dies here,
    ! against the form's declared output, before any Jacobian
    ! column is requested.
    !----------------------------------------------------------------!

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the materialization replaces it wholesale
    allocate(samples(0))

    call graph % build_samples(relation_index, samples)

    call extract_unknown_q(samples(unknown), q_values, ncomp)
    n = size(q_values)
    call q_star % set_real(q_values, ncomp=ncomp)

    signature = residual_form % output_signature()
    if (signature(1) * signature(2) /= n) then
       error stop 'gti_time_degree_r_directional_driver: degree residual size matches unknown size'
    end if

    call point_builder % build_point(graph % relation(relation_index) % motif, &
         & samples, design, graph % relation(relation_index) % evaluation_time, &
         & point)

    !----------------------------------------------------------------!
    ! The samples' derivative seats, degrees 1 through r - empty
    ! seats are zero derivatives.
    !----------------------------------------------------------------!

    allocate(local_path(n, arity, options % max_degree))
    do k = 1, options % max_degree
       do i = 1, arity
          vertex_index = graph % relation(relation_index) % sample_vertex(i)
          call get_derivative_or_zero(vertex_derivative(k, vertex_index), n, &
               & local_path(:, i, k))
       end do
    end do

    !----------------------------------------------------------------!
    ! The dense J_u, built ONCE from Newton's exact action, before
    ! the degree loop - and eliminated r times inside it. Only the
    ! right-hand side changes with the degree.
    !----------------------------------------------------------------!

    call build_unknown_jacobian(residual_form, &
         & graph % relation(relation_index) % motif, samples, unknown, &
         & q_star, design, graph % relation(relation_index) % evaluation_time, &
         & n, ncomp, jacobian)

    call initialize_step_result(result, relation_index, unknown_vertex, &
         & options % max_degree)

    !----------------------------------------------------------------!
    ! The degree loop: B^(s) with the unknown's top seat
    ! suppressed, one elimination, and the solved q^(s) lands both
    ! on the vertex seat and back in the local path for degree
    ! s + 1 to include.
    !----------------------------------------------------------------!

    do s = 1, options % max_degree

       local_path(:, unknown, s) = 0.0_dp

       call assemble_degree_rhs(residual_form, point, &
            & graph % relation(relation_index) % motif, arity, n, ncomp, &
            & local_path, design_direction, s, b_buffer, &
            & design_path, time_path)

       call b_buffer % get_real(b_values)
       if (size(b_values) /= n) then
          error stop 'gti_time_degree_r_directional_driver: degree residual size matches unknown size'
       end if

       rhs = -b_values
       call dense_solve(jacobian, rhs, options % singular_tolerance, solution)

       result % linear_residual_norm(s) = norm2(matmul(jacobian, solution) + b_values)
       call result % rhs(s) % set_real(rhs)
       call result % derivative(s) % set_real(solution, ncomp=ncomp)

       call set_vertex_derivative(vertex_derivative(s, unknown_vertex), &
            & solution, ncomp)

       local_path(:, unknown, s) = solution

    end do

  end subroutine solve_relation

  !===================================================================!
  ! The whole degree-r traversal: relations in stored forward
  ! order, each unknown's derivatives landing in the result's
  ! vertex array for later relations to read as history.
  !===================================================================!

  subroutine solve_all(this, residual_form, graph, design, design_direction, &
       & options, result, design_path, time_path)

    class(gti_time_degree_r_directional_driver) , intent(in)    :: this
    class(gti_differentiable_form)              , intent(in)    :: residual_form
    type(gti_time_graph)                        , intent(in)    :: graph
    type(gti_design_bundle)                     , intent(in)    :: design
    type(gti_value_buffer)                      , intent(in)    :: design_direction
    type(gti_time_degree_r_directional_options) , intent(in)    :: options
    type(gti_time_degree_r_directional_result)  , intent(inout) :: result
    type(gti_value_buffer)                      , intent(in), optional :: design_path(:)
    type(gti_value_buffer)                      , intent(in), optional :: time_path(:)

    integer :: r

    call options % validate()
    call graph % validate()

    if (graph % num_relations() < 1) then
       error stop 'gti_time_degree_r_directional_driver: at least one relation is required'
    end if

    result % completed = .false.
    result % completed_relations = 0
    result % max_degree = 0

    if (allocated(result % vertex_derivative)) deallocate(result % vertex_derivative)
    if (allocated(result % step))              deallocate(result % step)

    allocate(result % vertex_derivative(options % max_degree, &
         & graph % num_vertices()))
    allocate(result % step(graph % num_relations()))

    do r = 1, graph % num_relations()
       call this % solve_relation(residual_form, graph, r, design, &
            & design_direction, result % vertex_derivative, &
            & options, result % step(r), design_path, time_path)
       result % completed_relations = result % completed_relations + 1
    end do

    result % completed  = .true.
    result % max_degree = options % max_degree

  end subroutine solve_all

  !===================================================================!
  ! One degree's right-hand side, assembled by the chain-rule
  ! seat: each motif row m becomes one state channel carrying the
  ! component path derivatives
  !
  !      U_m^(k) = sum_i w_m(i) q_i^(k),        k = 1, ..., s,
  !
  ! read from the local path tensor - the caller has already
  ! applied the suppression rule there. A design channel and a
  ! time channel are added AFTER the state channels, each only
  ! when it carries at least one occupied seat up to this degree:
  ! with no design_path, the legacy affine channel occupies only
  ! x^(1) = eta; with a design_path, occupied design_path(k)
  ! entries occupy x^(k); with a time_path, occupied time_path(k)
  ! entries occupy t^(k). An absent seat, in every channel, is
  ! zero. This is the only place composition terms are assembled -
  ! gti_chain_rule_assemblies alone generates every pattern.
  !===================================================================!

  subroutine assemble_degree_rhs(residual_form, point, motif, arity, &
       & n, ncomp, local_path, design_direction, degree, b_buffer, &
       & design_path, time_path)

    class(gti_differentiable_form), intent(in)    :: residual_form
    type(gti_evaluation_point)    , intent(in)    :: point
    type(gti_time_motif)          , intent(in)    :: motif
    integer                       , intent(in)    :: arity, n, ncomp
    real(dp)                      , intent(in)    :: local_path(:,:,:)
    type(gti_value_buffer)        , intent(in)    :: design_direction
    integer                       , intent(in)    :: degree
    type(gti_value_buffer)        , intent(inout) :: b_buffer
    type(gti_value_buffer)        , intent(in), optional :: design_path(:)
    type(gti_value_buffer)        , intent(in), optional :: time_path(:)

    type(gti_chain_rule_assembler) :: assembler
    type(gti_chain_input)          :: input

    real(dp), allocatable :: component_path(:)
    integer :: m, i, k, nrules, nchannels, design_slot, time_slot
    logical :: add_design, add_time

    nrules = motif % size()

    !----------------------------------------------------------------!
    ! Which optional channels this degree earns: never an empty
    ! channel, only one that carries an occupied seat.
    !----------------------------------------------------------------!

    if (present(design_path)) then
       add_design = path_has_occupied(design_path, degree)
    else
       add_design = .true.
    end if

    add_time = .false.
    if (present(time_path)) add_time = path_has_occupied(time_path, degree)

    nchannels   = nrules
    design_slot = 0
    time_slot   = 0
    if (add_design) then
       nchannels   = nchannels + 1
       design_slot = nchannels
    end if
    if (add_time) then
       nchannels = nchannels + 1
       time_slot = nchannels
    end if

    allocate(input % channel(nchannels))
    allocate(component_path(n))

    do m = 1, nrules

       allocate(input % channel(m) % derivative(degree))

       do k = 1, degree

          component_path = 0.0_dp
          do i = 1, arity
             component_path = component_path + &
                  & motif % rule(m) % weight(i) * local_path(:, i, k)
          end do

          input % channel(m) % derivative(k) % argument_kind   = GTI_ARG_STATE
          input % channel(m) % derivative(k) % state_component = &
               & motif % rule(m) % state_component
          call input % channel(m) % derivative(k) % values % &
               & set_real(component_path, ncomp=ncomp)

       end do

    end do

    !----------------------------------------------------------------!
    ! The design channel: the legacy affine seat x^(1) = eta when
    ! no design_path was supplied, or the caller-supplied design
    ! seats otherwise.
    !----------------------------------------------------------------!

    if (add_design) then
       if (present(design_path)) then
          call fill_path_channel(input % channel(design_slot), &
               & GTI_ARG_DESIGN, design_path, degree)
       else
          call fill_legacy_design_channel(input % channel(design_slot), &
               & design_direction, degree)
       end if
    end if

    !----------------------------------------------------------------!
    ! The time channel: caller-supplied time seats only - there is
    ! no legacy time path to default to.
    !----------------------------------------------------------------!

    if (add_time) then
       call fill_path_channel(input % channel(time_slot), &
            & GTI_ARG_TIME, time_path, degree)
    end if

    call assembler % assemble(residual_form, point, degree, input, b_buffer)

  end subroutine assemble_degree_rhs

  !===================================================================!
  ! Does a path array carry at least one occupied seat at or below
  ! the requested degree? A seat beyond the array's size is absent,
  ! and an occupied entry is one holding real values.
  !===================================================================!

  pure function path_has_occupied(path, degree) result(has)

    type(gti_value_buffer), intent(in) :: path(:)
    integer                , intent(in) :: degree
    logical :: has

    integer :: k

    has = .false.
    do k = 1, min(degree, size(path))
       if (allocated(path(k) % rvals)) then
          if (size(path(k) % rvals) > 0) then
             has = .true.
             return
          end if
       end if
    end do

  end function path_has_occupied

  !===================================================================!
  ! Fill one channel from a caller-supplied path: seat k takes
  ! path(k) where it is occupied; a seat beyond the path's size, or
  ! an empty entry within it, stays absent - the channel's default,
  ! read as zero.
  !===================================================================!

  pure subroutine fill_path_channel(channel, argument_kind, path, degree)

    type(gti_chain_channel), intent(inout) :: channel
    integer                 , intent(in)    :: argument_kind
    type(gti_value_buffer)  , intent(in)    :: path(:)
    integer                 , intent(in)    :: degree

    real(dp), allocatable :: values(:)
    integer :: k

    allocate(channel % derivative(degree))

    do k = 1, min(degree, size(path))
       call path(k) % get_real(values)
       if (size(values) == 0) cycle
       channel % derivative(k) % argument_kind = argument_kind
       call channel % derivative(k) % values % set_real(values)
    end do

  end subroutine fill_path_channel

  !===================================================================!
  ! Fill the legacy affine design channel: x^(1) = eta, every
  ! higher seat absent - exactly the path xi + eps eta. Called only
  ! when the caller supplied no design_path.
  !===================================================================!

  pure subroutine fill_legacy_design_channel(channel, design_direction, degree)

    type(gti_chain_channel), intent(inout) :: channel
    type(gti_value_buffer) , intent(in)    :: design_direction
    integer                 , intent(in)    :: degree

    real(dp), allocatable :: eta_values(:)

    allocate(channel % derivative(degree))

    call design_direction % get_real(eta_values)
    channel % derivative(1) % argument_kind = GTI_ARG_DESIGN
    call channel % derivative(1) % values % set_real(eta_values)

  end subroutine fill_legacy_design_channel

  !===================================================================!
  ! The dense J_u from Newton's exact action, one basis direction
  ! per column. Built once per relation; no degree ever rebuilds
  ! it.
  !===================================================================!

  subroutine build_unknown_jacobian(residual_form, motif, samples, unknown, &
       & q_star, design, evaluation_time, n, ncomp, jacobian)

    class(gti_differentiable_form), intent(in)  :: residual_form
    type(gti_time_motif)          , intent(in)  :: motif
    type(gti_time_sample)         , intent(in)  :: samples(:)
    integer                       , intent(in)  :: unknown
    type(gti_value_buffer)        , intent(in)  :: q_star
    type(gti_design_bundle)       , intent(in)  :: design
    real(dp)                      , intent(in)  :: evaluation_time
    integer                       , intent(in)  :: n, ncomp
    real(dp), allocatable         , intent(out) :: jacobian(:,:)

    type(gti_time_local_newton_driver) :: newton
    type(gti_value_buffer)             :: dq_basis, column

    real(dp), allocatable :: basis(:), column_values(:)
    integer :: j

    allocate(jacobian(n, n), basis(n))

    do j = 1, n
       basis    = 0.0_dp
       basis(j) = 1.0_dp
       call dq_basis % set_real(basis, ncomp=ncomp)
       call newton % jacobian_action(residual_form, motif, samples, unknown, &
            & q_star, dq_basis, design, evaluation_time, column)
       call column % get_real(column_values)
       jacobian(:, j) = column_values
    end do

  end subroutine build_unknown_jacobian

  !===================================================================!
  ! The unknown sample's solved q, as values: a q that holds no
  ! real values cannot seed a derivative traversal.
  !===================================================================!

  subroutine extract_unknown_q(sample, q_values, ncomp)

    type(gti_time_sample), intent(in)  :: sample
    real(dp), allocatable , intent(out) :: q_values(:)
    integer               , intent(out) :: ncomp

    call sample % state % component(1 + GTI_STATE_Q) % value % &
         & get_real_vector(q_values)

    if (size(q_values) == 0) then
       error stop 'gti_time_degree_r_directional_driver: unknown q has values'
    end if

    ncomp = sample % state % component(1 + GTI_STATE_Q) % value % &
         & num_components()

  end subroutine extract_unknown_q

  !===================================================================!
  ! One vertex's derivative seat, read as values: an empty seat is
  ! a zero derivative; an occupied seat must wear the unknown's
  ! size.
  !===================================================================!

  subroutine get_derivative_or_zero(buffer, n, values)

    type(gti_value_buffer), intent(in)  :: buffer
    integer               , intent(in)  :: n
    real(dp)              , intent(out) :: values(:)

    real(dp), allocatable :: stored(:)

    call buffer % get_real(stored)

    if (size(stored) == 0) then
       values = 0.0_dp
       return
    end if

    if (size(stored) /= n) then
       error stop 'gti_time_degree_r_directional_driver: derivative shape matches vertex q'
    end if

    values = stored

  end subroutine get_derivative_or_zero

  !===================================================================!
  ! Land one solved derivative on its vertex seat.
  !===================================================================!

  subroutine set_vertex_derivative(buffer, values, ncomp)

    type(gti_value_buffer), intent(inout) :: buffer
    real(dp)              , intent(in)    :: values(:)
    integer               , intent(in)    :: ncomp

    call buffer % set_real(values, ncomp=ncomp)

  end subroutine set_vertex_derivative

  !===================================================================!
  ! One relation result, opened for r degrees: identity first,
  ! then the per-degree seats, sized once.
  !===================================================================!

  subroutine initialize_step_result(result, relation_index, unknown_vertex, &
       & max_degree)

    type(gti_time_degree_r_relation_result), intent(inout) :: result
    integer                                , intent(in)    :: relation_index
    integer                                , intent(in)    :: unknown_vertex
    integer                                , intent(in)    :: max_degree

    result % relation_index = relation_index
    result % unknown_vertex = unknown_vertex
    result % max_degree     = max_degree

    if (allocated(result % linear_residual_norm)) then
       deallocate(result % linear_residual_norm)
    end if
    if (allocated(result % rhs))        deallocate(result % rhs)
    if (allocated(result % derivative)) deallocate(result % derivative)

    allocate(result % linear_residual_norm(max_degree))
    allocate(result % rhs(max_degree))
    allocate(result % derivative(max_degree))

    result % linear_residual_norm = huge(1.0_dp)

  end subroutine initialize_step_result

  !===================================================================!
  ! Solve A x = b by Gaussian elimination with partial pivoting,
  ! on private copies. A pivot below the singular tolerance is
  ! refused. A deliberate local duplicate of the driver helpers -
  ! serving this driver alone, not an abstraction.
  !===================================================================!

  pure subroutine dense_solve(a, b, singular_tolerance, x)

    real(dp)             , intent(in)  :: a(:,:), b(:)
    real(dp)             , intent(in)  :: singular_tolerance
    real(dp), allocatable, intent(out) :: x(:)

    real(dp), allocatable :: m(:,:), r(:), row(:)
    real(dp) :: swap_value, factor
    integer :: n, k, p, i

    n = size(b)
    m = a
    r = b
    allocate(row(n))

    do k = 1, n

       p = k - 1 + maxloc(abs(m(k:n, k)), dim=1)

       if (abs(m(p, k)) <= singular_tolerance) then
          error stop 'gti_time_degree_r_directional_driver: dense Jacobian pivot is nonsingular'
       end if

       if (p /= k) then
          row        = m(k, :)
          m(k, :)    = m(p, :)
          m(p, :)    = row
          swap_value = r(k)
          r(k)       = r(p)
          r(p)       = swap_value
       end if

       do i = k + 1, n
          factor    = m(i, k) / m(k, k)
          m(i, k:n) = m(i, k:n) - factor * m(k, k:n)
          r(i)      = r(i) - factor * r(k)
       end do

    end do

    allocate(x(n))
    do k = n, 1, -1
       x(k) = (r(k) - dot_product(m(k, k+1:n), x(k+1:n))) / m(k, k)
    end do

  end subroutine dense_solve

end module gti_time_degree_r_directional_drivers
