!=====================================================================!
! GTI TIME FORWARD DRIVER
!
! The first executor over the explicit time graph: a relation-
! order forward march,
!
!      for r = 1, ..., num_relations:
!         samples = build_samples(r)
!         q0      = the unknown vertex's current q
!         q*      = local Newton solve on the relation's motif
!         write_unknown_q(r, q*)
!
! Each relation is one local implicit problem R_r(q_unknown) = 0
! over fixed history; the graph couples relations only through
! shared vertices - relation r writes vertex v, a later relation
! reads v as history - so the march is nothing but local solve,
! write shared vertex, and the next local solve sees the written
! vertex. Stored relation order IS the execution order: no
! topological scheduler, no adaptivity, no traversal abstraction
! exists here, and neither tangent nor adjoint is ever named.
!
! Before any relation is solved, its history must be solved: every
! non-unknown vertex the relation references carries has_solution.
! The unknown vertex's own q is the Newton initial guess, taken as
! a flat real vector - values required, one component per entry.
!
! Non-convergence is an answer, not an error: the failed step is
! reported, the march stops there, and the non-converged q is
! NEVER written back - the graph keeps the initial guess and the
! vertex stays unsolved.
!
! The driver carries nothing: no form, no graph, no solver state,
! no design, no schedule, no map.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_time_forward_drivers

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_design_bundles   , only : gti_design_bundle
  use gti_form_interface   , only : gti_differentiable_form
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_time_local_schemes, only : gti_time_sample
  use gti_time_local_newton_drivers, only : gti_time_local_newton_driver, &
       & gti_time_local_newton_options, gti_time_local_newton_result
  use gti_time_graphs      , only : gti_time_graph

  implicit none

  private
  public :: gti_time_forward_driver
  public :: gti_time_forward_options
  public :: gti_time_forward_step_result
  public :: gti_time_forward_result

  !===================================================================!
  ! The knobs of one forward march: exactly the local Newton
  ! knobs, and nothing of its own yet.
  !===================================================================!

  type :: gti_time_forward_options

     type(gti_time_local_newton_options) :: newton

   contains

     procedure :: validate => options_validate

  end type gti_time_forward_options

  !===================================================================!
  ! What one relation's solve reports: which relation, which
  ! global vertex it filled, and the local Newton account.
  !===================================================================!

  type :: gti_time_forward_step_result

     integer  :: relation_index = 0
     integer  :: unknown_vertex = 0
     logical  :: converged = .false.
     integer  :: iterations = 0
     real(dp) :: residual_norm = huge(1.0_dp)
     type(gti_value_buffer) :: q
     type(gti_value_buffer) :: residual

  end type gti_time_forward_step_result

  !===================================================================!
  ! What the whole march reports: how far it went, where it
  ! stopped if it stopped, and every step's account.
  !===================================================================!

  type :: gti_time_forward_result

     logical :: converged = .false.
     integer :: completed_relations = 0
     integer :: failed_relation = 0
     type(gti_time_forward_step_result), allocatable :: step(:)

  end type gti_time_forward_result

  !===================================================================!
  ! The stateless verb-pair. The types keep their public singular
  ! names; Fortran denies a type its host module's name, so the
  ! module speaks in the plural.
  !===================================================================!

  type :: gti_time_forward_driver

   contains

     procedure :: solve_relation
     procedure :: solve_all

  end type gti_time_forward_driver

contains

  !===================================================================!
  ! The march's knobs are lawful when Newton's are.
  !===================================================================!

  pure subroutine options_validate(this)

    class(gti_time_forward_options), intent(in) :: this

    call this % newton % validate()

  end subroutine options_validate

  !===================================================================!
  ! Solve exactly one relation: admit the graph, prove the history
  ! solved, materialize the local samples, take the unknown's
  ! current q as the initial guess, run the local Newton solve,
  ! and - only on convergence - write q* back to the shared
  ! vertex. A failed solve reports and writes nothing.
  !===================================================================!

  subroutine solve_relation(this, residual_form, graph, relation_index, &
       & design, options, result)

    class(gti_time_forward_driver)    , intent(in)    :: this
    class(gti_differentiable_form)    , intent(in)    :: residual_form
    type(gti_time_graph)              , intent(inout) :: graph
    integer                           , intent(in)    :: relation_index
    type(gti_design_bundle)           , intent(in)    :: design
    type(gti_time_forward_options)    , intent(in)    :: options
    type(gti_time_forward_step_result), intent(inout) :: result

    type(gti_time_local_newton_driver) :: newton
    type(gti_time_local_newton_result) :: newton_result
    type(gti_time_sample), allocatable :: samples(:)
    type(gti_value_buffer)             :: q_initial

    real(dp), allocatable :: q_values(:)
    integer :: k, unknown, vertex_index

    call options % validate()

    if (relation_index < 1 .or. relation_index > graph % num_relations()) then
       error stop 'gti_time_forward_driver: relation index is in range'
    end if

    call graph % validate()

    !----------------------------------------------------------------!
    ! The history law: every non-unknown vertex this relation
    ! reads must already carry a solution.
    !----------------------------------------------------------------!

    unknown = graph % relation(relation_index) % unknown_sample

    do k = 1, graph % relation(relation_index) % arity()
       if (k == unknown) cycle
       vertex_index = graph % relation(relation_index) % sample_vertex(k)
       if (.not. graph % vertex(vertex_index) % has_solution) then
          error stop 'gti_time_forward_driver: history vertex has solution'
       end if
    end do

    !----------------------------------------------------------------!
    ! Materialize, and take the unknown's current q as the Newton
    ! initial guess: values required, one component per entry.
    !----------------------------------------------------------------!

    ! a defined descriptor quiets a gfortran maybe-uninitialized
    ! false positive; the materialization replaces it wholesale
    allocate(samples(0))

    call graph % build_samples(relation_index, samples)

    call samples(unknown) % state % component(1 + GTI_STATE_Q) % value % &
         & get_real_vector(q_values)
    if (size(q_values) == 0) then
       error stop 'gti_time_forward_driver: unknown initial q has values'
    end if
    if (samples(unknown) % state % component(1 + GTI_STATE_Q) % value % &
         & num_components() /= 1) then
       error stop 'gti_time_forward_driver: unknown initial q is a vector'
    end if

    call q_initial % set_real(q_values)

    call newton % solve(residual_form, graph % relation(relation_index) % motif, &
         & samples, unknown, q_initial, design, &
         & graph % relation(relation_index) % evaluation_time, &
         & options % newton, newton_result)

    result % relation_index = relation_index
    result % unknown_vertex = graph % relation(relation_index) % unknown_vertex()
    result % converged      = newton_result % converged
    result % iterations     = newton_result % iterations
    result % residual_norm  = newton_result % residual_norm
    result % q              = newton_result % q
    result % residual       = newton_result % residual

    if (newton_result % converged) then
       call graph % write_unknown_q(relation_index, newton_result % q)
    end if

  end subroutine solve_relation

  !===================================================================!
  ! March every relation in stored order. The first failed step
  ! ends the march with its account on record; the graph keeps
  ! everything the converged steps wrote and nothing the failed
  ! step computed.
  !===================================================================!

  subroutine solve_all(this, residual_form, graph, design, options, result)

    class(gti_time_forward_driver), intent(in)    :: this
    class(gti_differentiable_form), intent(in)    :: residual_form
    type(gti_time_graph)          , intent(inout) :: graph
    type(gti_design_bundle)       , intent(in)    :: design
    type(gti_time_forward_options), intent(in)    :: options
    type(gti_time_forward_result) , intent(inout) :: result

    integer :: r

    call options % validate()
    call graph % validate()

    if (graph % num_relations() < 1) then
       error stop 'gti_time_forward_driver: at least one relation is required'
    end if

    result % converged = .false.
    result % completed_relations = 0
    result % failed_relation = 0

    if (allocated(result % step)) deallocate(result % step)
    allocate(result % step(graph % num_relations()))

    do r = 1, graph % num_relations()

       call this % solve_relation(residual_form, graph, r, design, options, &
            & result % step(r))

       if (.not. result % step(r) % converged) then
          result % converged = .false.
          result % failed_relation = r
          result % completed_relations = r - 1
          return
       end if

    end do

    result % converged = .true.
    result % completed_relations = graph % num_relations()
    result % failed_relation = 0

  end subroutine solve_all

end module gti_time_forward_drivers
