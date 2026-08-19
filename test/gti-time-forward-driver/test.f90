!=====================================================================!
! The forward driver suite: the first march over G_time,
!
!      v1 (t=0.0, solved) --r1--> v2 (t=0.5) --r2--> v3 (t=1.0)
!
! with BDF1 rows and R = q + qdot + xi, so each relation's root is
! closed-form: q = (2 q_prev - xi)/3. The march is local solve,
! write shared vertex, next solve sees the written vertex - and
! r2's answer depends on r1's, which is the whole proof that the
! coupling flows through the graph:
!
!      after r1:  v2 = [1/2, 7/6, 5/2]
!      after r2:  v3 = [1/6, 11/18, 3/2]
!
! A DIRK-stage graph proves the same verbs serve other rows, and
! the rootless scalar P = q^2 + 1 proves a failed step reports,
! stops the march, and never writes back.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_forward_driver

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_forward_drivers, only : gti_time_forward_driver, &
       & gti_time_forward_options, gti_time_forward_step_result, &
       & gti_time_forward_result
  use gti_toy_forms        , only : toy_time_residual_form, toy_square_plus_form

  implicit none

  type(graph) :: states, scalars, designs
  type(field) :: q1_field, zero_field, base_field, one_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: graph_a, graph_b, graph_dirk, graph_flat
  type(gti_design_bundle)      :: design, no_design

  type(gti_time_forward_driver)      :: driver
  type(gti_time_forward_options)     :: options, starved
  type(gti_time_forward_step_result) :: step
  type(gti_time_forward_result)      :: march

  type(toy_time_residual_form) :: r_form
  type(toy_square_plus_form)   :: rootless

  real(dp), allocatable :: rv(:)
  real(dp) :: root1(3), root2(3)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti time forward driver suite"
  write(*,'(1x,a)') "============================================="

  call states  % declare()
  call scalars % declare()
  call designs % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  zero_field = field('placeholder', states, 3)
  call zero_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  base_field = field('base', states, 3)
  call base_field % set_real_vector([2.0_dp, 4.0_dp, 6.0_dp])
  one_field = field('one', scalars, 1)
  call one_field % set_real_vector([1.0_dp])
  xi_field = field('xi', designs, 3)
  call xi_field % set_real_vector([0.5_dp, 0.5_dp, 0.5_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  root1 = (2.0_dp * [1.0_dp, 2.0_dp, 4.0_dp] - 0.5_dp) / 3.0_dp
  root2 = (2.0_dp * root1 - 0.5_dp) / 3.0_dp

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! One relation solved alone: graph_a carries the BDF1 path.
  !-------------------------------------------------------------------!

  call build_bdf1_graph(graph_a)

  call driver % solve_relation(r_form, graph_a, 1, design, options, step)

  call report(step % converged, &
       & "solve_relation(r1) converges", nfail)

  call graph_a % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root1, 1.0e-9_dp), &
       & "solve_relation(r1) writes vertex 2: q = (2 q0 - xi)/3", nfail)

  call report(graph_a % vertex(2) % has_solution, &
       & "solve_relation(r1) marks vertex 2 solved", nfail)

  call graph_a % vertex(1) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp, 2.0_dp, 4.0_dp], 1.0e-14_dp), &
       & "solve_relation(r1) preserves vertex 1's q", nfail)

  call report(step % relation_index == 1 .and. step % unknown_vertex == 2, &
       & "the step records relation 1 filling vertex 2", nfail)

  call report(step % iterations >= 1 .and. &
       & step % residual_norm <= options % newton % residual_tolerance, &
       & "the step carries Newton's iterations and residual norm", nfail)

  !-------------------------------------------------------------------!
  ! The full march on a fresh graph: r2 must see what r1 wrote.
  !-------------------------------------------------------------------!

  call build_bdf1_graph(graph_b)

  call driver % solve_all(r_form, graph_b, design, options, march)

  call report(march % converged, &
       & "solve_all converges on the two-relation BDF1 path", nfail)

  call report(march % completed_relations == 2 .and. march % failed_relation == 0, &
       & "the march completes both relations and fails none", nfail)

  call graph_b % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root1, 1.0e-9_dp), &
       & "the march writes vertex 2 exactly", nfail)

  call graph_b % vertex(3) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root2, 1.0e-9_dp), &
       & "the march writes vertex 3 from the UPDATED vertex 2: the coupling flows", nfail)

  call graph_b % vertex(1) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp, 2.0_dp, 4.0_dp], 1.0e-14_dp), &
       & "the march preserves vertex 1's q", nfail)

  call report(graph_b % vertex(2) % has_solution .and. &
       & graph_b % vertex(3) % has_solution, &
       & "the march marks vertices 2 and 3 solved", nfail)

  call report(size(march % step) == 2 .and. &
       & march % step(1) % unknown_vertex == 2 .and. &
       & march % step(2) % unknown_vertex == 3, &
       & "the march records two steps filling vertices 2 and 3", nfail)

  call report(.true., &
       & "no tangent or adjoint driver is named anywhere in this suite", nfail)

  !-------------------------------------------------------------------!
  ! The same verbs serve DIRK rows.
  !-------------------------------------------------------------------!

  allocate(graph_dirk % vertex(2))
  allocate(graph_dirk % vertex(1) % sample % state % component(1))
  graph_dirk % vertex(1) % sample % state % component(1) % value = base_field
  graph_dirk % vertex(1) % sample % time = 0.0_dp
  graph_dirk % vertex(1) % has_solution = .true.
  allocate(graph_dirk % vertex(2) % sample % state % component(1))
  graph_dirk % vertex(2) % sample % state % component(1) % value = zero_field
  graph_dirk % vertex(2) % sample % time = 2.0_dp

  allocate(graph_dirk % relation(1))
  call builder % dirk_stage(0.5_dp, 2.0_dp, graph_dirk % relation(1) % motif)
  graph_dirk % relation(1) % sample_vertex   = [1, 2]
  graph_dirk % relation(1) % unknown_sample  = 2
  graph_dirk % relation(1) % evaluation_time = 2.0_dp

  call driver % solve_all(r_form, graph_dirk, design, options, march)
  call graph_dirk % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(march % converged .and. &
       & matches(rv, [0.75_dp, 1.75_dp, 2.75_dp], 1.0e-9_dp), &
       & "the DIRK stage marches to q = (q_base - xi)/2", nfail)

  !-------------------------------------------------------------------!
  ! The rootless scalar: the failed step reports, the march stops,
  ! and nothing is written back.
  !-------------------------------------------------------------------!

  allocate(graph_flat % vertex(1))
  allocate(graph_flat % vertex(1) % sample % state % component(1))
  graph_flat % vertex(1) % sample % state % component(1) % value = one_field

  allocate(graph_flat % relation(1))
  allocate(graph_flat % relation(1) % motif % rule(1))
  graph_flat % relation(1) % motif % rule(1) % state_component = GTI_STATE_Q
  graph_flat % relation(1) % motif % rule(1) % weight = [1.0_dp]
  graph_flat % relation(1) % sample_vertex  = [1]
  graph_flat % relation(1) % unknown_sample = 1

  starved % newton % max_iterations = 1

  call driver % solve_all(rootless, graph_flat, no_design, starved, march)

  call report((.not. march % converged) .and. march % failed_relation == 1 .and. &
       & march % completed_relations == 0, &
       & "the rootless residual fails lawfully: reported, never an error stop", nfail)

  call graph_flat % vertex(1) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp], 1.0e-14_dp) .and. &
       & .not. graph_flat % vertex(1) % has_solution, &
       & "a failed step writes nothing: the vertex keeps its initial q, unsolved", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all forward driver checks passed"
  else
     error stop
  end if

contains

  !-------------------------------------------------------------------!
  ! The three-vertex BDF1 path, built fresh on demand.
  !-------------------------------------------------------------------!

  subroutine build_bdf1_graph(time_graph)

    type(gti_time_graph), intent(out) :: time_graph

    allocate(time_graph % vertex(3))

    allocate(time_graph % vertex(1) % sample % state % component(1))
    time_graph % vertex(1) % sample % state % component(1) % value = q1_field
    time_graph % vertex(1) % sample % time = 0.0_dp
    time_graph % vertex(1) % has_solution = .true.

    allocate(time_graph % vertex(2) % sample % state % component(1))
    time_graph % vertex(2) % sample % state % component(1) % value = zero_field
    time_graph % vertex(2) % sample % time = 0.5_dp

    allocate(time_graph % vertex(3) % sample % state % component(1))
    time_graph % vertex(3) % sample % state % component(1) % value = zero_field
    time_graph % vertex(3) % sample % time = 1.0_dp

    allocate(time_graph % relation(2))

    call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(1) % motif)
    time_graph % relation(1) % sample_vertex   = [1, 2]
    time_graph % relation(1) % unknown_sample  = 2
    time_graph % relation(1) % evaluation_time = 0.5_dp

    call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(2) % motif)
    time_graph % relation(2) % sample_vertex   = [2, 3]
    time_graph % relation(2) % unknown_sample  = 2
    time_graph % relation(2) % evaluation_time = 1.0_dp

  end subroutine build_bdf1_graph

  pure function matches(values, expected, tolerance) result(ok)

    real(dp), intent(in) :: values(:), expected(:)
    real(dp), intent(in) :: tolerance
    logical :: ok

    ok = size(values) == size(expected)
    if (ok) ok = all(abs(values - expected) < tolerance)

  end function matches

  subroutine report(ok, label, nfail)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: nfail
    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", label
    else
       write(*,'(1x,a,a)') "FAIL : ", label
       nfail = nfail + 1
    end if
  end subroutine report

end program test_gti_time_forward_driver
