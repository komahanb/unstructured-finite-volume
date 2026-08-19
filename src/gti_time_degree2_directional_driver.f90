!=====================================================================!
! GTI TIME DEGREE-2 DIRECTIONAL TRAVERSAL
!
! The first r >= 2 traversal over the time graph: given a solved
! G_time and one design direction eta, walk relations in stored
! forward order and compute the first and second directional
! derivatives of every unknown vertex,
!
!      J_u q_u^(1) = -B1        J_u q_u^(2) = -B2,
!
! with the SAME J_u Newton used - built once per relation from
! the local Jacobian action and eliminated twice. Only the
! chain-rule right-hand side changes with the degree; that is the
! architecture's core invariant, and this seat is where it earns
! its keep.
!
! B1 and B2 are assembled by the higher-order chain-rule seat,
! gti_chain_rule_assemblies, never re-derived here: each motif
! row m carries the component path
!
!      U_m^(s) = sum_i w_m(i) q_i^(s),
!
! its samples' derivatives read from the caller's vertex arrays -
! empty meaning zero - and the suppression rules select what the
! unknown contributes:
!
!      B1:  channel x^(1) with the unknown's q^(1) suppressed,
!           design x^(1) = eta
!      B2:  channel x^(1) FULL - the just-solved q_u^(1) included,
!           channel x^(2) with the unknown's q^(2) suppressed,
!           design x^(1) = eta, design x^(2) absent (zero)
!
! so degree-1 assembly is B1 and degree-2 assembly - transport of
! x^(2) plus curvature over ordered pairs of x^(1) - is B2,
! exactly. No derivative is finite-differenced, and no
! combinatorics is duplicated: what the assembler owns stays
! owned.
!
! The traversal reads the graph and writes only the caller's
! derivative arrays: no q value and no solution flag is touched,
! no reverse pass, no functional, no controller exists here.
!
! The driver carries nothing: no graph, no forms, no solver
! state, no design, no derivatives, no map.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_degree2_directional_drivers

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_direction_bundle, GTI_ARG_STATE, GTI_ARG_DESIGN
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_time_local_schemes, only : gti_time_sample, gti_time_motif, &
       & gti_time_local_residual_evaluator, gti_evaluation_point
  use gti_time_local_newton_drivers, only : gti_time_local_newton_driver
  use gti_time_graphs      , only : gti_time_graph
  use gti_chain_rule_assemblies, only : gti_chain_rule_assembler, &
       & gti_chain_channel, gti_chain_input

  implicit none

  private
  public :: gti_time_degree2_directional_options
  public :: gti_time_degree2_relation_result
  public :: gti_time_degree2_directional_result
  public :: gti_time_degree2_directional_driver

  !===================================================================!
  ! The one knob of the degree-2 traversal.
  !===================================================================!

  type :: gti_time_degree2_directional_options

     real(dp) :: singular_tolerance = 1.0e-14_dp

   contains

     procedure :: validate => options_validate

  end type gti_time_degree2_directional_options

  !===================================================================!
  ! What one relation's two solves report: the right-hand sides,
  ! the derivatives, and how exactly each linear system was
  ! satisfied.
  !===================================================================!

  type :: gti_time_degree2_relation_result

     integer  :: relation_index = 0
     integer  :: unknown_vertex = 0
     real(dp) :: first_linear_residual_norm  = huge(1.0_dp)
     real(dp) :: second_linear_residual_norm = huge(1.0_dp)
     type(gti_value_buffer) :: first_rhs
     type(gti_value_buffer) :: second_rhs
     type(gti_value_buffer) :: q_first
     type(gti_value_buffer) :: q_second

  end type gti_time_degree2_relation_result

  !===================================================================!
  ! What the whole traversal reports: completion, every step's
  ! account, and the per-vertex derivative arrays.
  !===================================================================!

  type :: gti_time_degree2_directional_result

     logical :: completed = .false.
     integer :: completed_relations = 0
     type(gti_value_buffer), allocatable :: vertex_first(:)
     type(gti_value_buffer), allocatable :: vertex_second(:)
     type(gti_time_degree2_relation_result), allocatable :: step(:)

  end type gti_time_degree2_directional_result

  !===================================================================!
  ! The stateless verb-pair. The types keep their public singular
  ! names; Fortran denies a type its host module's name, so the
  ! module speaks in the plural.
  !===================================================================!

  type :: gti_time_degree2_directional_driver

   contains

     procedure :: solve_relation
     procedure :: solve_all

  end type gti_time_degree2_directional_driver

contains

  pure subroutine options_validate(this)

    class(gti_time_degree2_directional_options), intent(in) :: this

    if (this % singular_tolerance <= 0.0_dp) then
       error stop 'gti_time_degree2_directional_driver: singular tolerance is positive'
    end if

  end subroutine options_validate

  !===================================================================!
  ! One relation's degree-1 and degree-2 solves: the same J_u,
  ! eliminated twice, against the chain-rule right-hand sides.
  !===================================================================!

  subroutine solve_relation(this, residual_form, graph, relation_index, &
       & design, design_direction, vertex_first, vertex_second, options, result)

    class(gti_time_degree2_directional_driver), intent(in)    :: this
    class(gti_differentiable_form)            , intent(in)    :: residual_form
    type(gti_time_graph)                      , intent(in)    :: graph
    integer                                   , intent(in)    :: relation_index
    type(gti_design_bundle)                   , intent(in)    :: design
    type(gti_value_buffer)                    , intent(in)    :: design_direction
    type(gti_value_buffer)                    , intent(inout) :: vertex_first(:)
    type(gti_value_buffer)                    , intent(inout) :: vertex_second(:)
    type(gti_time_degree2_directional_options), intent(in)    :: options
    type(gti_time_degree2_relation_result)    , intent(inout) :: result

    type(gti_time_local_newton_driver)      :: newton
    type(gti_time_local_residual_evaluator) :: point_builder
    type(gti_time_sample), allocatable      :: samples(:)
    type(gti_evaluation_point)              :: point
    type(gti_value_buffer)                  :: q_star, dq_basis, column, b_buffer

    real(dp), allocatable :: q_values(:), eta_values(:)
    real(dp), allocatable :: first_local(:,:), second_local(:,:)
    real(dp), allocatable :: jacobian(:,:), column_values(:), basis(:)
    real(dp), allocatable :: b_values(:), rhs(:), solution(:)
    integer :: n, ncomp, arity, unknown, unknown_vertex
    integer :: i, j, k, vertex_index

    call options % validate()
    call graph % validate()

    if (relation_index < 1 .or. relation_index > graph % num_relations()) then
       error stop 'gti_time_degree2_directional_driver: relation index is in range'
    end if

    if (size(vertex_first) /= graph % num_vertices() .or. &
         & size(vertex_second) /= graph % num_vertices()) then
       error stop 'gti_time_degree2_directional_driver: derivative array size matches graph vertices'
    end if

    k = 0
    if (allocated(design_direction % rvals)) k = size(design_direction % rvals)
    if (k == 0) then
       error stop 'gti_time_degree2_directional_driver: design direction has values'
    end if
    call design_direction % get_real(eta_values)

    arity          = graph % relation(relation_index) % arity()
    unknown        = graph % relation(relation_index) % unknown_sample
    unknown_vertex = graph % relation(relation_index) % unknown_vertex()

    if (.not. graph % vertex(unknown_vertex) % has_solution) then
       error stop 'gti_time_degree2_directional_driver: unknown vertex has solution'
    end if

    do i = 1, arity
       if (i == unknown) cycle
       vertex_index = graph % relation(relation_index) % sample_vertex(i)
       if (.not. graph % vertex(vertex_index) % has_solution) then
          error stop 'gti_time_degree2_directional_driver: history vertex has solution'
       end if
    end do

    !----------------------------------------------------------------!
    ! Materialize the solved local state and the samples'
    ! derivative seats - empty seats are zero derivatives.
    !----------------------------------------------------------------!

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the materialization replaces it wholesale
    allocate(samples(0))

    call graph % build_samples(relation_index, samples)

    call samples(unknown) % state % component(1 + GTI_STATE_Q) % value % &
         & get_real_vector(q_values)
    n = size(q_values)
    if (n == 0) then
       error stop 'gti_time_degree2_directional_driver: unknown q has values'
    end if
    ncomp = samples(unknown) % state % component(1 + GTI_STATE_Q) % value % &
         & num_components()
    call q_star % set_real(q_values, ncomp=ncomp)

    call point_builder % build_point(graph % relation(relation_index) % motif, &
         & samples, design, graph % relation(relation_index) % evaluation_time, &
         & point)

    allocate(first_local(n, arity), second_local(n, arity))
    do i = 1, arity
       vertex_index = graph % relation(relation_index) % sample_vertex(i)
       call get_seed_or_zero(vertex_first(vertex_index), n, first_local(:, i))
       call get_seed_or_zero(vertex_second(vertex_index), n, second_local(:, i))
    end do

    !----------------------------------------------------------------!
    ! B1: the total first derivative with the unknown's q^(1)
    ! suppressed - assembled by the chain-rule seat at degree 1.
    !----------------------------------------------------------------!

    first_local(:, unknown) = 0.0_dp

    call assemble_degree_rhs(residual_form, point, &
         & graph % relation(relation_index) % motif, unknown, arity, n, ncomp, &
         & first_local, second_local, eta_values, 1, b_buffer)

    call b_buffer % get_real(b_values)
    if (size(b_values) /= n) then
       error stop 'gti_time_degree2_directional_driver: degree residual size matches unknown size'
    end if

    !----------------------------------------------------------------!
    ! The dense J_u, built ONCE from Newton's exact action - and
    ! eliminated twice below. Only the right-hand side changes
    ! with the degree.
    !----------------------------------------------------------------!

    allocate(jacobian(n, n), basis(n))

    do j = 1, n
       basis    = 0.0_dp
       basis(j) = 1.0_dp
       call dq_basis % set_real(basis, ncomp=ncomp)
       call newton % jacobian_action(residual_form, &
            & graph % relation(relation_index) % motif, samples, unknown, &
            & q_star, dq_basis, design, &
            & graph % relation(relation_index) % evaluation_time, column)
       call column % get_real(column_values)
       jacobian(:, j) = column_values
    end do

    rhs = -b_values
    call dense_solve(jacobian, rhs, options % singular_tolerance, solution)

    result % relation_index = relation_index
    result % unknown_vertex = unknown_vertex
    result % first_linear_residual_norm = norm2(matmul(jacobian, solution) + b_values)
    call result % first_rhs % set_real(rhs)
    call result % q_first % set_real(solution, ncomp=ncomp)

    call set_vertex_derivative(vertex_first(unknown_vertex), solution, ncomp)

    !----------------------------------------------------------------!
    ! B2: the total second derivative with the unknown's q^(2)
    ! suppressed and the just-solved q^(1) INCLUDED - assembled by
    ! the chain-rule seat at degree 2.
    !----------------------------------------------------------------!

    first_local(:, unknown)  = solution
    second_local(:, unknown) = 0.0_dp

    call assemble_degree_rhs(residual_form, point, &
         & graph % relation(relation_index) % motif, unknown, arity, n, ncomp, &
         & first_local, second_local, eta_values, 2, b_buffer)

    call b_buffer % get_real(b_values)
    if (size(b_values) /= n) then
       error stop 'gti_time_degree2_directional_driver: degree residual size matches unknown size'
    end if

    rhs = -b_values
    call dense_solve(jacobian, rhs, options % singular_tolerance, solution)

    result % second_linear_residual_norm = norm2(matmul(jacobian, solution) + b_values)
    call result % second_rhs % set_real(rhs)
    call result % q_second % set_real(solution, ncomp=ncomp)

    call set_vertex_derivative(vertex_second(unknown_vertex), solution, ncomp)

  end subroutine solve_relation

  !===================================================================!
  ! The whole degree-2 traversal: relations in stored forward
  ! order, each unknown's first and second derivatives landing in
  ! the result's vertex arrays for later relations to read as
  ! history.
  !===================================================================!

  subroutine solve_all(this, residual_form, graph, design, design_direction, &
       & options, result)

    class(gti_time_degree2_directional_driver), intent(in)    :: this
    class(gti_differentiable_form)            , intent(in)    :: residual_form
    type(gti_time_graph)                      , intent(in)    :: graph
    type(gti_design_bundle)                   , intent(in)    :: design
    type(gti_value_buffer)                    , intent(in)    :: design_direction
    type(gti_time_degree2_directional_options), intent(in)    :: options
    type(gti_time_degree2_directional_result) , intent(inout) :: result

    integer :: r

    call options % validate()
    call graph % validate()

    if (graph % num_relations() < 1) then
       error stop 'gti_time_degree2_directional_driver: at least one relation is required'
    end if

    result % completed = .false.
    result % completed_relations = 0

    if (allocated(result % vertex_first))  deallocate(result % vertex_first)
    if (allocated(result % vertex_second)) deallocate(result % vertex_second)
    if (allocated(result % step))          deallocate(result % step)

    allocate(result % vertex_first(graph % num_vertices()))
    allocate(result % vertex_second(graph % num_vertices()))
    allocate(result % step(graph % num_relations()))

    do r = 1, graph % num_relations()
       call this % solve_relation(residual_form, graph, r, design, &
            & design_direction, result % vertex_first, result % vertex_second, &
            & options, result % step(r))
       result % completed_relations = result % completed_relations + 1
    end do

    result % completed = .true.

  end subroutine solve_all

  !===================================================================!
  ! One degree's right-hand side, assembled by the chain-rule
  ! seat: each motif row m becomes one state channel carrying the
  ! component path derivatives
  !
  !      U_m^(s) = sum_i w_m(i) q_i^(s),
  !
  ! read from the local derivative matrices - the caller has
  ! already applied the suppression rules there - and the design
  ! channel carries eta as x^(1) with x^(2) absent. Degree 1
  ! yields B1; degree 2 yields B2.
  !===================================================================!

  subroutine assemble_degree_rhs(residual_form, point, motif, unknown, arity, &
       & n, ncomp, first_local, second_local, eta_values, degree, b_buffer)

    class(gti_differentiable_form), intent(in)    :: residual_form
    type(gti_evaluation_point)    , intent(in)    :: point
    type(gti_time_motif)          , intent(in)    :: motif
    integer                       , intent(in)    :: unknown, arity, n, ncomp
    real(dp)                      , intent(in)    :: first_local(:,:)
    real(dp)                      , intent(in)    :: second_local(:,:)
    real(dp)                      , intent(in)    :: eta_values(:)
    integer                       , intent(in)    :: degree
    type(gti_value_buffer)        , intent(inout) :: b_buffer

    type(gti_chain_rule_assembler) :: assembler
    type(gti_chain_input)          :: input

    real(dp), allocatable :: component_path(:)
    integer :: m, i, nrules, seats

    ! the unknown index participates through the caller-prepared
    ! local matrices; it is not consulted again here
    associate(unread => unknown)
    end associate

    nrules = motif % size()
    seats  = degree

    allocate(input % channel(nrules + 1))
    allocate(component_path(n))

    do m = 1, nrules

       allocate(input % channel(m) % derivative(seats))

       component_path = 0.0_dp
       do i = 1, arity
          component_path = component_path + &
               & motif % rule(m) % weight(i) * first_local(:, i)
       end do
       input % channel(m) % derivative(1) % argument_kind   = GTI_ARG_STATE
       input % channel(m) % derivative(1) % state_component = &
            & motif % rule(m) % state_component
       call input % channel(m) % derivative(1) % values % &
            & set_real(component_path, ncomp=ncomp)

       if (degree >= 2) then
          component_path = 0.0_dp
          do i = 1, arity
             component_path = component_path + &
                  & motif % rule(m) % weight(i) * second_local(:, i)
          end do
          input % channel(m) % derivative(2) % argument_kind   = GTI_ARG_STATE
          input % channel(m) % derivative(2) % state_component = &
               & motif % rule(m) % state_component
          call input % channel(m) % derivative(2) % values % &
               & set_real(component_path, ncomp=ncomp)
       end if

    end do

    !----------------------------------------------------------------!
    ! The design channel: eta as x^(1); x^(2) stays an empty seat,
    ! and absence is zero.
    !----------------------------------------------------------------!

    allocate(input % channel(nrules + 1) % derivative(seats))
    input % channel(nrules + 1) % derivative(1) % argument_kind = GTI_ARG_DESIGN
    call input % channel(nrules + 1) % derivative(1) % values % &
         & set_real(eta_values)

    call assembler % assemble(residual_form, point, degree, input, b_buffer)

  end subroutine assemble_degree_rhs

  !===================================================================!
  ! One vertex's derivative seat, read as values: an empty seat is
  ! a zero derivative; an occupied seat must wear the unknown's
  ! size.
  !===================================================================!

  subroutine get_seed_or_zero(buffer, n, values)

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
       error stop 'gti_time_degree2_directional_driver: derivative shape matches vertex q'
    end if

    values = stored

  end subroutine get_seed_or_zero

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
          error stop 'gti_time_degree2_directional_driver: dense Jacobian pivot is nonsingular'
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

end module gti_time_degree2_directional_drivers
