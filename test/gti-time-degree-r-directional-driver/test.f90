!=====================================================================!
! The degree-r directional traversal suite: over the solved
! nonlinear BDF1 chain
!
!      R = (q_u - q_h) + q_u^2 - xi,   q = 1 everywhere, xi = 1,
!
! ONE driver walks the graph at r = 1, 2, 3, 4 - and then at
! r = 8, past every retired table - computing every directional
! derivative along eta = 1 from the same J_u = 3 per relation,
! built once and eliminated r times:
!
!      r1:  q2^(1) =  1/3          (rhs      1)
!           q2^(2) = -2/27         (rhs  -2/9)
!           q2^(3) =  4/81         (rhs   4/27)
!           q2^(4) = -40/729       (rhs -40/243)
!      r2:  q3^(1) =  4/9          (rhs      4/3)
!           q3^(2) = -38/243       (rhs    -38/81)
!           q3^(3) =  340/2187     (rhs    340/729)
!           q3^(4) = -14848/59049  (rhs -14848/19683)
!
! Though the residual's own calculus stops at the quadratic, the
! third and fourth derivatives are alive: the chain-rule
! right-hand side transports products of lower-degree curvature,
! assembled by gti_chain_rule_assemblies at every degree. Degree
! 3 and degree 4 are verification cases of the SAME driver, not
! drivers of their own; at r = 2 the general loop must reproduce
! the dedicated degree-2 seat exactly. The graph's primal content
! survives every traversal untouched, and every number arrives
! through partial actions alone.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_degree_r_directional_driver

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_degree_r_directional_drivers, only : &
       & gti_time_degree_r_directional_driver, &
       & gti_time_degree_r_directional_options, &
       & gti_time_degree_r_directional_result, &
       & gti_time_degree_r_relation_result
  use gti_time_degree2_directional_drivers, only : &
       & gti_time_degree2_directional_driver, &
       & gti_time_degree2_directional_options, &
       & gti_time_degree2_directional_result
  use gti_toy_forms        , only : toy_qdot_square_form, &
       & toy_qdot_square_time_form

  implicit none

  type(graph) :: states, designs
  type(field) :: one_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph
  type(gti_design_bundle)      :: design

  type(gti_time_degree_r_directional_driver)  :: driver
  type(gti_time_degree_r_directional_options) :: options
  type(gti_time_degree_r_directional_result)  :: result

  type(gti_time_degree2_directional_driver)  :: old_driver
  type(gti_time_degree2_directional_options) :: old_options
  type(gti_time_degree2_directional_result)  :: old_result

  type(toy_qdot_square_form) :: r_form
  type(gti_value_buffer)     :: eta

  real(dp), allocatable :: rv(:), ov(:)
  logical :: norms_ok, agrees
  integer :: v, s8, nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti degree-r directional traversal suite"
  write(*,'(1x,a)') "============================================="

  call states  % declare()
  call designs % declare()

  one_field = field('one', states, 1)
  call one_field % set_real_vector([1.0_dp])
  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([1.0_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  call eta % set_real([1.0_dp])

  !-------------------------------------------------------------------!
  ! The solved chain: q = 1 at every vertex IS the root of
  ! (q_u - q_h) + q_u^2 - xi at xi = 1.
  !-------------------------------------------------------------------!

  allocate(time_graph % vertex(3))
  do v = 1, 3
     allocate(time_graph % vertex(v) % sample % state % component(1))
     time_graph % vertex(v) % sample % state % component(1) % value = one_field
     time_graph % vertex(v) % sample % time = real(v - 1, dp)
     time_graph % vertex(v) % has_solution = .true.
  end do

  allocate(time_graph % relation(2))
  call builder % bdf_uniform(1, 1.0_dp, time_graph % relation(1) % motif)
  time_graph % relation(1) % sample_vertex   = [1, 2]
  time_graph % relation(1) % unknown_sample  = 2
  time_graph % relation(1) % evaluation_time = 1.0_dp
  call builder % bdf_uniform(1, 1.0_dp, time_graph % relation(2) % motif)
  time_graph % relation(2) % sample_vertex   = [2, 3]
  time_graph % relation(2) % unknown_sample  = 2
  time_graph % relation(2) % evaluation_time = 2.0_dp

  norms_ok = .true.

  !-------------------------------------------------------------------!
  ! r = 1: the traversal is the tangent march.
  !-------------------------------------------------------------------!

  options % max_degree = 1
  call driver % solve_all(r_form, time_graph, design, eta, options, result)
  call fold_norms(result, norms_ok)

  call report(result % completed, &
       & "r=1: solve_all completes", nfail)

  call report(result % max_degree == 1, &
       & "r=1: the result records max_degree = 1", nfail)

  call report(result % completed_relations == 2, &
       & "r=1: both relations completed, vertices 2 then 3", nfail)

  call result % vertex_derivative(1, 2) % get_real(rv)
  call report(matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "r=1: q2^(1) = 1/3", nfail)

  call result % vertex_derivative(1, 3) % get_real(rv)
  call report(matches(rv, [4.0_dp / 9.0_dp], 1.0e-12_dp), &
       & "r=1: q3^(1) = 4/9, v2's derivative fed r2 as history", nfail)

  call report(size(result % vertex_derivative, 1) == 1 .and. &
       & size(result % step(1) % rhs) == 1 .and. &
       & size(result % step(1) % derivative) == 1 .and. &
       & size(result % step(1) % linear_residual_norm) == 1, &
       & "r=1: no degree-2/3/4 storage exists in the result", nfail)

  !-------------------------------------------------------------------!
  ! r = 2: the general loop must BE the dedicated degree-2 seat.
  !-------------------------------------------------------------------!

  options % max_degree = 2
  call driver % solve_all(r_form, time_graph, design, eta, options, result)
  call fold_norms(result, norms_ok)

  call result % vertex_derivative(1, 2) % get_real(rv)
  call report(matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "r=2: q2^(1) = 1/3", nfail)

  call result % vertex_derivative(2, 2) % get_real(rv)
  call report(matches(rv, [-2.0_dp / 27.0_dp], 1.0e-12_dp), &
       & "r=2: q2^(2) = -2/27", nfail)

  call result % vertex_derivative(1, 3) % get_real(rv)
  call report(matches(rv, [4.0_dp / 9.0_dp], 1.0e-12_dp), &
       & "r=2: q3^(1) = 4/9", nfail)

  call result % vertex_derivative(2, 3) % get_real(rv)
  call report(matches(rv, [-38.0_dp / 243.0_dp], 1.0e-12_dp), &
       & "r=2: q3^(2) = -38/243", nfail)

  call old_driver % solve_all(r_form, time_graph, design, eta, &
       & old_options, old_result)

  agrees = .true.
  call result % vertex_derivative(1, 2) % get_real(rv)
  call old_result % vertex_first(2) % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % vertex_derivative(2, 2) % get_real(rv)
  call old_result % vertex_second(2) % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % vertex_derivative(1, 3) % get_real(rv)
  call old_result % vertex_first(3) % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % vertex_derivative(2, 3) % get_real(rv)
  call old_result % vertex_second(3) % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % step(1) % rhs(1) % get_real(rv)
  call old_result % step(1) % first_rhs % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % step(1) % rhs(2) % get_real(rv)
  call old_result % step(1) % second_rhs % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % step(2) % rhs(1) % get_real(rv)
  call old_result % step(2) % first_rhs % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)
  call result % step(2) % rhs(2) % get_real(rv)
  call old_result % step(2) % second_rhs % get_real(ov)
  agrees = agrees .and. matches(rv, ov, 1.0e-14_dp)

  call report(agrees, &
       & "r=2: the general loop matches the dedicated degree-2 driver", nfail)

  !-------------------------------------------------------------------!
  ! r = 3: alive though the residual's third partials are zero.
  !-------------------------------------------------------------------!

  options % max_degree = 3
  call driver % solve_all(r_form, time_graph, design, eta, options, result)
  call fold_norms(result, norms_ok)

  call result % vertex_derivative(3, 2) % get_real(rv)
  call report(matches(rv, [4.0_dp / 81.0_dp], 1.0e-12_dp), &
       & "r=3: q2^(3) = 4/81, lower-degree curvature transported", nfail)

  call result % vertex_derivative(3, 3) % get_real(rv)
  call report(matches(rv, [340.0_dp / 2187.0_dp], 1.0e-12_dp), &
       & "r=3: q3^(3) = 340/2187", nfail)

  call result % step(1) % rhs(3) % get_real(rv)
  call report(matches(rv, [4.0_dp / 27.0_dp], 1.0e-12_dp), &
       & "r=3: r1's degree-3 right-hand side is -B3 = [4/27]", nfail)

  call result % step(2) % rhs(3) % get_real(rv)
  call report(matches(rv, [340.0_dp / 729.0_dp], 1.0e-12_dp), &
       & "r=3: r2's degree-3 right-hand side is -B3 = [340/729]", nfail)

  !-------------------------------------------------------------------!
  ! r = 4: the top of the assembler's vocabulary.
  !-------------------------------------------------------------------!

  options % max_degree = 4
  call driver % solve_all(r_form, time_graph, design, eta, options, result)
  call fold_norms(result, norms_ok)

  call result % vertex_derivative(4, 2) % get_real(rv)
  call report(matches(rv, [-40.0_dp / 729.0_dp], 1.0e-12_dp), &
       & "r=4: q2^(4) = -40/729", nfail)

  call result % vertex_derivative(4, 3) % get_real(rv)
  call report(matches(rv, [-14848.0_dp / 59049.0_dp], 1.0e-12_dp), &
       & "r=4: q3^(4) = -14848/59049", nfail)

  call result % step(1) % rhs(4) % get_real(rv)
  call report(matches(rv, [-40.0_dp / 243.0_dp], 1.0e-12_dp), &
       & "r=4: r1's degree-4 right-hand side is -B4 = [-40/243]", nfail)

  call result % step(2) % rhs(4) % get_real(rv)
  call report(matches(rv, [-14848.0_dp / 19683.0_dp], 1.0e-12_dp), &
       & "r=4: r2's degree-4 right-hand side is -B4 = [-14848/19683]", nfail)

  !-------------------------------------------------------------------!
  ! r = 8: four turns past every retired table. The residual's own
  ! calculus still stops at the quadratic; every value below is
  ! transported curvature, and the same J_u = 3 is eliminated
  ! eight times per relation. Exact values by rational Taylor
  ! composition of q_u(eps) = (-1 + sqrt(9 + 4 eps))/2 and its
  ! feed-forward into relation 2.
  !-------------------------------------------------------------------!

  options % max_degree = 8
  call driver % solve_all(r_form, time_graph, design, eta, options, result)
  call fold_norms(result, norms_ok)

  call report(result % completed .and. result % max_degree == 8 .and. &
       & result % completed_relations == 2, &
       & "r=8: solve_all completes and records max_degree = 8", nfail)

  call report(size(result % vertex_derivative, 1) == 8 .and. &
       & size(result % vertex_derivative, 2) == 3 .and. &
       & size(result % step(1) % rhs) == 8 .and. &
       & size(result % step(1) % derivative) == 8 .and. &
       & size(result % step(1) % linear_residual_norm) == 8, &
       & "r=8: every result array carries eight degree seats", nfail)

  call result % vertex_derivative(5, 2) % get_real(rv)
  call report(matches(rv, [560.0_dp / 6561.0_dp], 1.0e-10_dp), &
       & "r=8: q2^(5) = 560/6561", nfail)

  call result % vertex_derivative(6, 2) % get_real(rv)
  call report(matches(rv, [-1120.0_dp / 6561.0_dp], 1.0e-10_dp), &
       & "r=8: q2^(6) = -1120/6561", nfail)

  call result % vertex_derivative(7, 2) % get_real(rv)
  call report(matches(rv, [24640.0_dp / 59049.0_dp], 1.0e-10_dp), &
       & "r=8: q2^(7) = 24640/59049", nfail)

  call result % vertex_derivative(8, 2) % get_real(rv)
  call report(matches(rv, [-640640.0_dp / 531441.0_dp], 1.0e-10_dp), &
       & "r=8: q2^(8) = -640640/531441", nfail)

  call result % vertex_derivative(5, 3) % get_real(rv)
  call report(matches(rv, [897680.0_dp / 1594323.0_dp], 1.0e-10_dp), &
       & "r=8: q3^(5) = 897680/1594323", nfail)

  call result % vertex_derivative(6, 3) % get_real(rv)
  call report(matches(rv, [-95200.0_dp / 59049.0_dp], 1.0e-10_dp), &
       & "r=8: q3^(6) = -95200/59049", nfail)

  call result % vertex_derivative(7, 3) % get_real(rv)
  call report(matches(rv, [726772480.0_dp / 129140163.0_dp], 1.0e-10_dp), &
       & "r=8: q3^(7) = 726772480/129140163", nfail)

  call result % vertex_derivative(8, 3) % get_real(rv)
  call report(matches(rv, [-80862642560.0_dp / 3486784401.0_dp], 1.0e-10_dp), &
       & "r=8: q3^(8) = -80862642560/3486784401", nfail)

  agrees = .true.
  do v = 1, 2
     do s8 = 1, 8
        call result % step(v) % rhs(s8) % get_real(rv)
        call result % step(v) % derivative(s8) % get_real(ov)
        agrees = agrees .and. matches(rv, 3.0_dp * ov, 1.0e-10_dp)
     end do
  end do
  call report(agrees, &
       & "r=8: every stored rhs equals J_u q^(s) = 3 q^(s), both relations", nfail)

  !-------------------------------------------------------------------!
  ! Design-path law: the hard-coded xi^(1) = eta, xi^(s>=2) absent
  ! occupancy is retired. A caller-supplied design_path can occupy
  ! xi^(2), which the legacy path could never reach.
  !
  !      R = (q_u - q_h) + q_u^2 - xi,   q_h = 1, xi = 1, q_u = 1
  !
  ! One relation. design_path: xi^(1) = 1, xi^(2) = 6, xi^(s>=3)
  ! absent - degree 1 reproduces the legacy answer since xi^(1)
  ! matches eta; degree 2 reaches a seat the legacy path never
  ! occupies.
  !
  !      degree 1: 3 q1 - xi1 = 0             => q1 = 1/3
  !      degree 2: 3 q2 + 2 q1^2 - xi2 = 0     => q2 = 52/27
  !-------------------------------------------------------------------!

  block
    type(gti_time_graph) :: path_graph
    type(gti_time_degree_r_directional_options) :: path_options
    type(gti_time_degree_r_relation_result)     :: path_step
    type(gti_value_buffer) :: vd(2, 2)
    type(gti_value_buffer) :: design_path(2), empty_design_path(1)
    type(gti_value_buffer) :: no_eta
    integer :: pv

    allocate(path_graph % vertex(2))
    do pv = 1, 2
       allocate(path_graph % vertex(pv) % sample % state % component(1))
       path_graph % vertex(pv) % sample % state % component(1) % value = one_field
       path_graph % vertex(pv) % has_solution = .true.
    end do
    allocate(path_graph % relation(1))
    call builder % bdf_uniform(1, 1.0_dp, path_graph % relation(1) % motif)
    path_graph % relation(1) % sample_vertex   = [1, 2]
    path_graph % relation(1) % unknown_sample  = 2
    path_graph % relation(1) % evaluation_time = 1.0_dp

    call design_path(1) % set_real([1.0_dp])
    call design_path(2) % set_real([6.0_dp])

    path_options % max_degree = 2
    call driver % solve_relation(r_form, path_graph, 1, design, eta, vd, &
         & path_options, path_step, design_path=design_path)

    call path_step % derivative(1) % get_real(rv)
    call report(matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp), &
         & "design-path: q1 = 1/3 with xi1 = 1", nfail)

    call path_step % derivative(2) % get_real(rv)
    call report(matches(rv, [52.0_dp / 27.0_dp], 1.0e-12_dp), &
         & "design-path: q2 = 52/27 with xi2 = 6, no longer hard-coded absent", nfail)

    !----------------------------------------------------------------!
    ! The pass-through law: solve_all forwards design_path to
    ! solve_relation untouched - proven by reaching the same seat
    ! through the whole-graph verb, and by the vertex_derivative
    ! write-back the direct solve_relation checks above do not
    ! exercise.
    !----------------------------------------------------------------!

    block
      type(gti_time_degree_r_directional_result) :: path_result
      call driver % solve_all(r_form, path_graph, design, eta, &
           & path_options, path_result, design_path=design_path)
      call path_result % vertex_derivative(1, 2) % get_real(rv)
      call report(matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp), &
           & "solve_all forwards design_path: q1 = 1/3", nfail)
      call path_result % vertex_derivative(2, 2) % get_real(rv)
      call report(matches(rv, [52.0_dp / 27.0_dp], 1.0e-12_dp), &
           & "solve_all forwards design_path: q2 = 52/27, through the whole-graph verb", nfail)
    end block

    !----------------------------------------------------------------!
    ! With design_path absent, the legacy mechanism is untouched:
    ! q2 = -2/27 still holds for the same relation, same eta.
    !----------------------------------------------------------------!

    call driver % solve_relation(r_form, path_graph, 1, design, eta, vd, &
         & path_options, path_step)

    call path_step % derivative(1) % get_real(rv)
    call report(matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp), &
         & "design-path absent: legacy q1 = 1/3 unaffected", nfail)

    call path_step % derivative(2) % get_real(rv)
    call report(matches(rv, [-2.0_dp / 27.0_dp], 1.0e-12_dp), &
         & "design-path absent: legacy q2 = -2/27 recovered exactly", nfail)

    !----------------------------------------------------------------!
    ! Edge law: design_path shorter than max_degree - only the
    ! provided seat is occupied; the missing degree-2 seat is
    ! absent, exactly as the legacy path leaves it.
    !----------------------------------------------------------------!

    call driver % solve_relation(r_form, path_graph, 1, design, eta, vd, &
         & path_options, path_step, design_path=design_path(1:1))

    call path_step % derivative(2) % get_real(rv)
    call report(matches(rv, [-2.0_dp / 27.0_dp], 1.0e-12_dp), &
         & "design-path shorter than max_degree: missing xi2 is absent", nfail)

    !----------------------------------------------------------------!
    ! Edge law: design_path present but empty - no design channel
    ! is added, and the design contribution is exactly zero.
    !----------------------------------------------------------------!

    path_options % max_degree = 1
    call driver % solve_relation(r_form, path_graph, 1, design, eta, vd(1:1, :), &
         & path_options, path_step, design_path=empty_design_path)

    call path_step % derivative(1) % get_real(rv)
    call report(matches(rv, [0.0_dp], 1.0e-12_dp), &
         & "design-path present but empty: no design channel, q1 = 0", nfail)

    !----------------------------------------------------------------!
    ! Edge law: design_path present allows an empty legacy
    ! design_direction - the caller need not carry both.
    !----------------------------------------------------------------!

    call no_eta % clear()
    path_options % max_degree = 2
    call driver % solve_relation(r_form, path_graph, 1, design, no_eta, vd, &
         & path_options, path_step, design_path=design_path)

    call path_step % derivative(1) % get_real(rv)
    agrees = matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp)
    call path_step % derivative(2) % get_real(rv)
    agrees = agrees .and. matches(rv, [52.0_dp / 27.0_dp], 1.0e-12_dp)
    call report(agrees, &
         & "design-path present: an empty legacy design_direction is lawful", nfail)

  end block

  !-------------------------------------------------------------------!
  ! Time-path law: a time derivative channel now reaches the chain
  ! assembler, alongside a design channel, without conflating the
  ! two.
  !
  !      R(q_u, q_h, xi, t) = (q_u - q_h) + q_u^2 - xi - t
  !      q_h = 1, xi = 1, t = 0, q_u = 1
  !
  ! One relation. design_path: xi^(1) = 1, xi^(s>=2) absent.
  ! time_path: t^(1) = 2, t^(2) = 12, t^(s>=3) absent.
  !
  !      degree 1: 3 q1 - xi1 - t1 = 0        => q1 = (1 + 2) / 3 = 1
  !      degree 2: 3 q2 + 2 q1^2 - t2 = 0      => q2 = (12 - 2) / 3 = 10/3
  !-------------------------------------------------------------------!

  block
    type(gti_time_graph)            :: t_graph
    type(toy_qdot_square_time_form) :: t_form
    type(gti_time_degree_r_directional_options) :: t_options
    type(gti_time_degree_r_relation_result)     :: t_step
    type(gti_value_buffer) :: vd(2, 2)
    type(gti_value_buffer) :: design_path(1), time_path(2), empty_time_path(1)
    integer :: pv

    allocate(t_graph % vertex(2))
    do pv = 1, 2
       allocate(t_graph % vertex(pv) % sample % state % component(1))
       t_graph % vertex(pv) % sample % state % component(1) % value = one_field
       t_graph % vertex(pv) % has_solution = .true.
    end do
    allocate(t_graph % relation(1))
    call builder % bdf_uniform(1, 1.0_dp, t_graph % relation(1) % motif)
    t_graph % relation(1) % sample_vertex   = [1, 2]
    t_graph % relation(1) % unknown_sample  = 2
    t_graph % relation(1) % evaluation_time = 0.0_dp

    call design_path(1) % set_real([1.0_dp])
    call time_path(1) % set_real([2.0_dp])
    call time_path(2) % set_real([12.0_dp])

    !----------------------------------------------------------------!
    ! Without time_path: the time-aware form recovers exactly the
    ! design-only answer - the time channel is simply never added.
    !----------------------------------------------------------------!

    t_options % max_degree = 1
    call driver % solve_relation(t_form, t_graph, 1, design, eta, vd(1:1, :), &
         & t_options, t_step, design_path=design_path)

    call t_step % derivative(1) % get_real(rv)
    call report(matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp), &
         & "time-path absent: the design-only answer q1 = 1/3 is recovered", nfail)

    !----------------------------------------------------------------!
    ! With time_path: t1 = 2 and t2 = 12 both reach the assembler.
    !----------------------------------------------------------------!

    t_options % max_degree = 2
    call driver % solve_relation(t_form, t_graph, 1, design, eta, vd, &
         & t_options, t_step, design_path=design_path, time_path=time_path)

    call t_step % derivative(1) % get_real(rv)
    call report(matches(rv, [1.0_dp], 1.0e-12_dp), &
         & "time-path: q1 = 1 with xi1 = 1, t1 = 2", nfail)

    call t_step % derivative(2) % get_real(rv)
    call report(matches(rv, [10.0_dp / 3.0_dp], 1.0e-12_dp), &
         & "time-path: q2 = 10/3 with t2 = 12, xi2 absent", nfail)

    !----------------------------------------------------------------!
    ! The pass-through law: solve_all forwards both design_path and
    ! time_path to solve_relation untouched.
    !----------------------------------------------------------------!

    block
      type(gti_time_degree_r_directional_result) :: t_result
      call driver % solve_all(t_form, t_graph, design, eta, &
           & t_options, t_result, design_path=design_path, time_path=time_path)
      call t_result % vertex_derivative(1, 2) % get_real(rv)
      call report(matches(rv, [1.0_dp], 1.0e-12_dp), &
           & "solve_all forwards design_path and time_path: q1 = 1", nfail)
      call t_result % vertex_derivative(2, 2) % get_real(rv)
      call report(matches(rv, [10.0_dp / 3.0_dp], 1.0e-12_dp), &
           & "solve_all forwards design_path and time_path: q2 = 10/3", nfail)
    end block

    !----------------------------------------------------------------!
    ! Edge law: time_path shorter than max_degree - only t1 is
    ! occupied; t2 is absent, so degree 2 sees no design and no
    ! time seat at all: 3 q2 + 2 q1^2 = 0, q1 = 1 => q2 = -2/3.
    !----------------------------------------------------------------!

    call driver % solve_relation(t_form, t_graph, 1, design, eta, vd, &
         & t_options, t_step, design_path=design_path, &
         & time_path=time_path(1:1))

    call t_step % derivative(1) % get_real(rv)
    call report(matches(rv, [1.0_dp], 1.0e-12_dp), &
         & "time-path shorter than max_degree: q1 = 1 from the occupied t1", nfail)

    call t_step % derivative(2) % get_real(rv)
    call report(matches(rv, [-2.0_dp / 3.0_dp], 1.0e-12_dp), &
         & "time-path shorter than max_degree: missing t2 is absent", nfail)

    !----------------------------------------------------------------!
    ! Edge law: time_path present but empty - no time channel is
    ! added, and the time contribution is exactly zero.
    !----------------------------------------------------------------!

    t_options % max_degree = 1
    call driver % solve_relation(t_form, t_graph, 1, design, eta, vd(1:1, :), &
         & t_options, t_step, design_path=design_path, &
         & time_path=empty_time_path)

    call t_step % derivative(1) % get_real(rv)
    call report(matches(rv, [1.0_dp / 3.0_dp], 1.0e-12_dp), &
         & "time-path present but empty: no time channel, design-only answer holds", nfail)

  end block

  !-------------------------------------------------------------------!
  ! Shared laws: exact eliminations, and a graph that never felt
  ! the traversals.
  !-------------------------------------------------------------------!

  call report(norms_ok, &
       & "every linear residual norm at every degree is below 1e-12", nfail)

  do v = 1, 3
     call time_graph % vertex(v) % sample % state % component(1 + GTI_STATE_Q) % &
          & value % get_real_vector(rv)
     if (.not. matches(rv, [1.0_dp], 1.0e-14_dp)) exit
  end do
  call report(v == 4, &
       & "the graph's q values survive all five traversals untouched", nfail)

  call report(time_graph % vertex(1) % has_solution .and. &
       & time_graph % vertex(2) % has_solution .and. &
       & time_graph % vertex(3) % has_solution, &
       & "the graph's solution flags survive all five traversals untouched", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all degree-r traversal checks passed"
  else
     error stop
  end if

contains

  pure function matches(values, expected, tolerance) result(ok)

    real(dp), intent(in) :: values(:), expected(:)
    real(dp), intent(in) :: tolerance
    logical :: ok

    ok = size(values) == size(expected)
    if (ok) ok = all(abs(values - expected) < tolerance)

  end function matches

  subroutine fold_norms(result, norms_ok)
    type(gti_time_degree_r_directional_result), intent(in)    :: result
    logical                                    , intent(inout) :: norms_ok
    integer :: r, s
    do r = 1, size(result % step)
       do s = 1, size(result % step(r) % linear_residual_norm)
          norms_ok = norms_ok .and. &
               & result % step(r) % linear_residual_norm(s) <= 1.0e-12_dp
       end do
    end do
  end subroutine fold_norms

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

end program test_gti_time_degree_r_directional_driver
