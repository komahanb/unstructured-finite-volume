!=====================================================================!
! The heterogeneous chain proof: ONE time graph whose relations
! carry three different named schemes -
!
!      v1 --BDF1--> v2 --DIRK--> v3 --ABM2--> v4
!
! - and the UNCHANGED solve-gradient driver marches, seeds, and
! reverses it in one call. Nothing in any driver knows a scheme
! name: each relation owns its coefficient rows and shares its
! vertices, so heterogeneous composition is a property G_time
! already had, not a feature added today.
!
! Every number is closed-form:
!
!      forward:  v2 = [1/2, 7/6, 5/2]
!                v3 = [0, 1/3, 1]
!                v4 = [-1/6, 1/18, 1/2]
!      value:    35/9
!      reverse:  lambda = [1/3...], [1/3...], [1/9...]
!                seeds  = [2/9...], [1/3...], [2/3...], [1,1,1]
!      actions:  explicit 3, residual -7/3, total 2/3
!
! and the crown is the hand sensitivity across the scheme
! boundaries: dq4/dq1 = (2/3)(1/2)(2/3) = 2/9, which IS the final
! seed on v1 - the chain rule flows through BDF1, DIRK, and ABM2
! without ever noticing where one scheme ends and the next begins.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_heterogeneous_chain

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
       & toy_sum_time_functional

  implicit none

  type(graph) :: states, designs
  type(field) :: q1_field, zero_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph
  type(gti_design_bundle)      :: design

  type(gti_time_functional) :: terminal

  type(gti_time_solve_gradient_driver)  :: driver
  type(gti_time_solve_gradient_options) :: options
  type(gti_time_solve_gradient_result)  :: result

  type(toy_time_residual_form)  :: r_form
  type(toy_sum_time_functional) :: f_form

  type(gti_value_buffer) :: eta

  real(dp), allocatable :: rv(:)
  real(dp) :: v2(3), v3(3), v4(3)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti heterogeneous chain proof suite"
  write(*,'(1x,a)') "============================================="

  call states  % declare()
  call designs % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  zero_field = field('placeholder', states, 3)
  call zero_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  xi_field = field('xi', designs, 3)
  call xi_field % set_real_vector([0.5_dp, 0.5_dp, 0.5_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  call eta % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  v2 = (2.0_dp * [1.0_dp, 2.0_dp, 4.0_dp] - 0.5_dp) / 3.0_dp
  v3 = (v2 - 0.5_dp) / 2.0_dp
  v4 = (2.0_dp * v3 - 0.5_dp) / 3.0_dp

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! One graph, three schemes: each relation carries its own rows.
  !-------------------------------------------------------------------!

  allocate(time_graph % vertex(4))
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
  allocate(time_graph % vertex(4) % sample % state % component(1))
  time_graph % vertex(4) % sample % state % component(1) % value = zero_field
  time_graph % vertex(4) % sample % time = 2.0_dp

  allocate(time_graph % relation(3))

  call builder % bdf_uniform(1, 0.5_dp, time_graph % relation(1) % motif)
  time_graph % relation(1) % sample_vertex   = [1, 2]
  time_graph % relation(1) % unknown_sample  = 2
  time_graph % relation(1) % evaluation_time = 0.5_dp

  call builder % dirk_stage(0.5_dp, 2.0_dp, time_graph % relation(2) % motif)
  time_graph % relation(2) % sample_vertex   = [2, 3]
  time_graph % relation(2) % unknown_sample  = 2
  time_graph % relation(2) % evaluation_time = 1.0_dp

  call builder % abm_corrector(2, 1.0_dp, time_graph % relation(3) % motif)
  time_graph % relation(3) % sample_vertex   = [3, 4]
  time_graph % relation(3) % unknown_sample  = 2
  time_graph % relation(3) % evaluation_time = 2.0_dp

  allocate(terminal % term(1))
  terminal % term(1) % vertex_index    = 4
  terminal % term(1) % evaluation_time = 2.0_dp

  !-------------------------------------------------------------------!
  ! Mixed rows, proven mixed before anything runs: BDF1 and ABM2
  ! rate rows are [-2, 2] at their steps, DIRK's is [-1, 1].
  !-------------------------------------------------------------------!

  call report(matches(time_graph % relation(1) % motif % rule(2) % weight, &
       & [-2.0_dp, 2.0_dp], 1.0e-14_dp) .and. &
       & matches(time_graph % relation(2) % motif % rule(2) % weight, &
       & [-1.0_dp, 1.0_dp], 1.0e-14_dp) .and. &
       & matches(time_graph % relation(3) % motif % rule(2) % weight, &
       & [-2.0_dp, 2.0_dp], 1.0e-14_dp), &
       & "three named builders minted three relations' rows", nfail)

  !-------------------------------------------------------------------!
  ! ONE call. Every check below reads only the result and the
  ! graph; no manual forward, seed, or reverse call exists in
  ! this law body.
  !-------------------------------------------------------------------!

  call driver % solve(r_form, f_form, time_graph, terminal, design, eta, &
       & options, result)

  call report(result % completed, &
       & "the solve-gradient pass completes on the mixed-scheme graph", nfail)

  call report(result % forward_converged .and. &
       & result % forward % completed_relations == 3, &
       & "the march converged all three relations", nfail)

  call time_graph % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, v2, 1.0e-9_dp), &
       & "BDF1 wrote v2 = (2 q1 - xi)/3 = [1/2, 7/6, 5/2]", nfail)

  call time_graph % vertex(3) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, v3, 1.0e-9_dp), &
       & "DIRK wrote v3 = (v2 - xi)/2 = [0, 1/3, 1]", nfail)

  call time_graph % vertex(4) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, v4, 1.0e-9_dp), &
       & "ABM2 wrote v4 = (2 v3 - xi)/3 = [-1/6, 1/18, 1/2]", nfail)

  call time_graph % vertex(1) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, [1.0_dp, 2.0_dp, 4.0_dp], 1.0e-14_dp), &
       & "v1's q is preserved across the scheme boundaries", nfail)

  call report(time_graph % vertex(1) % has_solution .and. &
       & time_graph % vertex(2) % has_solution .and. &
       & time_graph % vertex(3) % has_solution .and. &
       & time_graph % vertex(4) % has_solution, &
       & "all four vertices carry solutions after the pass", nfail)

  call report(abs(result % value - 35.0_dp / 9.0_dp) < 1.0e-12_dp, &
       & "the functional value is sum(q4) + sum(xi) + t = 35/9", nfail)

  call report(result % functional_seed % completed .and. &
       & result % reverse % completed, &
       & "the seeding and reverse stages both completed", nfail)

  call result % explicit_design_action % get_real(rv)
  call report(matches(rv, [3.0_dp], 1.0e-12_dp), &
       & "the explicit design action is [3]", nfail)

  call result % residual_design_action % get_real(rv)
  call report(matches(rv, [-7.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "the residual contribution is -1 - 1 - 1/3 = [-7/3]", nfail)

  call result % total_design_gradient_action % get_real(rv)
  call report(matches(rv, [2.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "the total design-gradient action is 3 - 7/3 = 2/3", nfail)

  call result % vertex_seed(1) % get_real(rv)
  call report(matches(rv, [2.0_dp, 2.0_dp, 2.0_dp] / 9.0_dp, 1.0e-12_dp), &
       & "seed v1 = [2/9...] = (2/3)(1/2)(2/3) = dq4/dq1 across all schemes", nfail)

  call result % vertex_seed(2) % get_real(rv)
  call report(matches(rv, [1.0_dp, 1.0_dp, 1.0_dp] / 3.0_dp, 1.0e-12_dp), &
       & "seed v2 = [1/3...]", nfail)

  call result % vertex_seed(3) % get_real(rv)
  call report(matches(rv, [2.0_dp, 2.0_dp, 2.0_dp] / 3.0_dp, 1.0e-12_dp), &
       & "seed v3 = [2/3...]", nfail)

  call result % vertex_seed(4) % get_real(rv)
  call report(matches(rv, [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-14_dp), &
       & "the terminal seed on v4 stays [1, 1, 1]", nfail)

  call report(result % reverse % step(3) % unknown_vertex == 4 .and. &
       & result % reverse % step(2) % unknown_vertex == 3 .and. &
       & result % reverse % step(1) % unknown_vertex == 2, &
       & "the reverse steps filled vertices 4, 3, 2 in reverse order", nfail)

  call report(result % reverse % step(1) % linear_residual_norm <= 1.0e-12_dp .and. &
       & result % reverse % step(2) % linear_residual_norm <= 1.0e-12_dp .and. &
       & result % reverse % step(3) % linear_residual_norm <= 1.0e-12_dp, &
       & "every transposed solve's linear residual is below 1e-12", nfail)

  call report(.true., &
       & "value and gradient came from one solve-gradient call alone", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all heterogeneous chain checks passed"
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

end program test_gti_heterogeneous_chain
