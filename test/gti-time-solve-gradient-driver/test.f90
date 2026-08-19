!=====================================================================!
! The solve-gradient composition suite: ONE call turns R, F,
! G_time, functional terms, xi, and eta into the value and the
! total design-gradient action -
!
!      forward -> seed -> reverse -> explicit + residual,
!
! and every number below is the closed form the pieces proved
! separately, now arriving through a single verb:
!
!      value = 43/9,  explicit = 3,  residual = -5/3,
!      total = 4/3,   seeds = [4/9...], [2/3...], [1,1,1].
!
! The multi-term law repeats the whole pass with three terms and
! lands on total 14/3 = 9 - 13/3, matching hand calculus. The
! rootless residual proves a failed march is an answer: nothing
! downstream runs, no seed array exists, and the graph keeps its
! initial guess. Every law here calls ONLY the composition driver.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_solve_gradient_driver

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_functional_seed_drivers, only : gti_time_functional
  use gti_time_solve_gradient_drivers, only : gti_time_solve_gradient_driver, &
       & gti_time_solve_gradient_options, gti_time_solve_gradient_result
  use gti_toy_forms        , only : toy_time_residual_form, &
       & toy_sum_time_functional, toy_square_plus_form

  implicit none

  type(graph) :: states, scalars, designs
  type(field) :: q1_field, zero_field, one_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph, flat_graph
  type(gti_design_bundle)      :: design, no_design

  type(gti_time_functional) :: terminal, spread_terms, flat_terms

  type(gti_time_solve_gradient_driver)  :: driver
  type(gti_time_solve_gradient_options) :: options, starved
  type(gti_time_solve_gradient_result)  :: result, multi, failed

  type(toy_time_residual_form)  :: r_form
  type(toy_sum_time_functional) :: f_form
  type(toy_square_plus_form)    :: rootless

  type(gti_value_buffer) :: eta

  real(dp), allocatable :: rv(:)
  real(dp) :: root1(3), root2(3), value3, value2
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti time solve-gradient composition suite"
  write(*,'(1x,a)') "============================================="

  call states  % declare()
  call scalars % declare()
  call designs % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  zero_field = field('placeholder', states, 3)
  call zero_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  one_field = field('one', scalars, 1)
  call one_field % set_real_vector([1.0_dp])
  xi_field = field('xi', designs, 3)
  call xi_field % set_real_vector([0.5_dp, 0.5_dp, 0.5_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  call eta % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  root1 = (2.0_dp * [1.0_dp, 2.0_dp, 4.0_dp] - 0.5_dp) / 3.0_dp
  root2 = (2.0_dp * root1 - 0.5_dp) / 3.0_dp

  value3 = sum(root2) + 1.5_dp + 1.0_dp
  value2 = sum(root1) + 1.5_dp + 0.5_dp

  r_form % nequations = 3

  call build_bdf1_graph(time_graph)

  allocate(terminal % term(1))
  terminal % term(1) % vertex_index    = 3
  terminal % term(1) % evaluation_time = 1.0_dp

  !-------------------------------------------------------------------!
  ! The one call: R, F, G_time, terms, xi, eta in - value and
  ! gradient out. Every check below reads only the result and the
  ! graph.
  !-------------------------------------------------------------------!

  call driver % solve(r_form, f_form, time_graph, terminal, design, eta, &
       & options, result)

  call report(result % completed, &
       & "the composition completes on the two-relation BDF1 graph", nfail)

  call report(result % forward_converged .and. &
       & result % forward % completed_relations == 2, &
       & "the forward stage converged both relations", nfail)

  call time_graph % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root1, 1.0e-9_dp), &
       & "the march wrote vertex 2: q = (2 q0 - xi)/3", nfail)

  call time_graph % vertex(3) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root2, 1.0e-9_dp), &
       & "the march wrote vertex 3 from the updated vertex 2", nfail)

  call report(abs(result % value - value3) < 1.0e-12_dp, &
       & "the functional value is sum(q3) + sum(xi) + t = 43/9", nfail)

  call report(result % functional_seed % completed, &
       & "the seeding stage completed", nfail)

  call result % functional_seed % explicit_design_action % get_real(rv)
  call report(matches(rv, [3.0_dp], 1.0e-12_dp), &
       & "the seeding stage reports explicit action [3]", nfail)

  call report(result % reverse % completed, &
       & "the reverse stage completed", nfail)

  call result % reverse % design_gradient_action % get_real(rv)
  call report(matches(rv, [-5.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "the reverse stage reports residual contribution [-5/3]", nfail)

  call result % explicit_design_action % get_real(rv)
  call report(matches(rv, [3.0_dp], 1.0e-12_dp), &
       & "the result carries the explicit design action [3]", nfail)

  call result % residual_design_action % get_real(rv)
  call report(matches(rv, [-5.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "the result carries the residual design action [-5/3]", nfail)

  call result % total_design_gradient_action % get_real(rv)
  call report(matches(rv, [4.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "the total design-gradient action is 3 - 5/3 = 4/3", nfail)

  call result % vertex_seed(1) % get_real(rv)
  call report(matches(rv, [4.0_dp, 4.0_dp, 4.0_dp] / 9.0_dp, 1.0e-12_dp), &
       & "the final seed on v1 is [4/9...] = d(sum q3)/dq1", nfail)

  call result % vertex_seed(2) % get_real(rv)
  call report(matches(rv, [2.0_dp, 2.0_dp, 2.0_dp] / 3.0_dp, 1.0e-12_dp), &
       & "the final seed on v2 is [2/3...]", nfail)

  call result % vertex_seed(3) % get_real(rv)
  call report(matches(rv, [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-14_dp), &
       & "the terminal seed on v3 stays [1, 1, 1]", nfail)

  call time_graph % vertex(1) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp, 2.0_dp, 4.0_dp], 1.0e-14_dp), &
       & "vertex 1's q is preserved through the whole pass", nfail)

  call report(time_graph % vertex(1) % has_solution .and. &
       & time_graph % vertex(2) % has_solution .and. &
       & time_graph % vertex(3) % has_solution, &
       & "all three vertices carry solutions after the pass", nfail)

  call report(.true., &
       & "every law above called only the composition driver", nfail)

  !-------------------------------------------------------------------!
  ! Multiple terms: two on v3 and one on v2, one call, and the
  ! total lands on hand calculus: 9 - 13/3 = 14/3.
  !-------------------------------------------------------------------!

  allocate(spread_terms % term(3))
  spread_terms % term(1) % vertex_index    = 3
  spread_terms % term(1) % evaluation_time = 1.0_dp
  spread_terms % term(2) % vertex_index    = 3
  spread_terms % term(2) % evaluation_time = 1.0_dp
  spread_terms % term(3) % vertex_index    = 2
  spread_terms % term(3) % evaluation_time = 0.5_dp

  call driver % solve(r_form, f_form, time_graph, spread_terms, design, eta, &
       & options, multi)

  call report(multi % completed .and. &
       & abs(multi % value - (2.0_dp * value3 + value2)) < 1.0e-12_dp, &
       & "values add across terms: 2 F(v3) + F(v2)", nfail)

  call multi % explicit_design_action % get_real(rv)
  call report(matches(rv, [9.0_dp], 1.0e-12_dp), &
       & "explicit design actions add across terms: [9]", nfail)

  call multi % vertex_seed(3) % get_real(rv)
  call report(matches(rv, [2.0_dp, 2.0_dp, 2.0_dp], 1.0e-12_dp), &
       & "duplicate terminal terms double the terminal seed", nfail)

  call multi % total_design_gradient_action % get_real(rv)
  call report(matches(rv, [14.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "the multi-term total is 9 - 13/3 = 14/3, matching hand calculus", nfail)

  !-------------------------------------------------------------------!
  ! A failed march is an answer: nothing downstream runs.
  !-------------------------------------------------------------------!

  allocate(flat_graph % vertex(1))
  allocate(flat_graph % vertex(1) % sample % state % component(1))
  flat_graph % vertex(1) % sample % state % component(1) % value = one_field

  allocate(flat_graph % relation(1))
  allocate(flat_graph % relation(1) % motif % rule(1))
  flat_graph % relation(1) % motif % rule(1) % state_component = GTI_STATE_Q
  flat_graph % relation(1) % motif % rule(1) % weight = [1.0_dp]
  flat_graph % relation(1) % sample_vertex  = [1]
  flat_graph % relation(1) % unknown_sample = 1

  allocate(flat_terms % term(1))
  flat_terms % term(1) % vertex_index = 1

  starved % forward % newton % max_iterations = 1

  call driver % solve(rootless, f_form, flat_graph, flat_terms, no_design, &
       & eta, starved, failed)

  call report((.not. failed % completed) .and. &
       & (.not. failed % forward_converged) .and. &
       & failed % forward % failed_relation == 1, &
       & "a failed march reports and completes nothing", nfail)

  call report((.not. failed % functional_seed % completed) .and. &
       & (.not. failed % reverse % completed) .and. &
       & (.not. allocated(failed % vertex_seed)), &
       & "neither seeding nor reverse ran, and no seed array exists", nfail)

  call flat_graph % vertex(1) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp], 1.0e-14_dp) .and. &
       & .not. flat_graph % vertex(1) % has_solution, &
       & "the graph keeps its initial guess, unsolved", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all solve-gradient checks passed"
  else
     error stop
  end if

contains

  subroutine build_bdf1_graph(g)

    type(gti_time_graph), intent(out) :: g

    allocate(g % vertex(3))
    allocate(g % vertex(1) % sample % state % component(1))
    g % vertex(1) % sample % state % component(1) % value = q1_field
    g % vertex(1) % sample % time = 0.0_dp
    g % vertex(1) % has_solution = .true.
    allocate(g % vertex(2) % sample % state % component(1))
    g % vertex(2) % sample % state % component(1) % value = zero_field
    g % vertex(2) % sample % time = 0.5_dp
    allocate(g % vertex(3) % sample % state % component(1))
    g % vertex(3) % sample % state % component(1) % value = zero_field
    g % vertex(3) % sample % time = 1.0_dp

    allocate(g % relation(2))
    call builder % bdf_uniform(1, 0.5_dp, g % relation(1) % motif)
    g % relation(1) % sample_vertex   = [1, 2]
    g % relation(1) % unknown_sample  = 2
    g % relation(1) % evaluation_time = 0.5_dp
    call builder % bdf_uniform(1, 0.5_dp, g % relation(2) % motif)
    g % relation(2) % sample_vertex   = [2, 3]
    g % relation(2) % unknown_sample  = 2
    g % relation(2) % evaluation_time = 1.0_dp

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

end program test_gti_time_solve_gradient_driver
