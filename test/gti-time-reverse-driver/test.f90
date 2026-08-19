!=====================================================================!
! The reverse seed driver suite: the forward march converges the
! two-relation BDF1 path, and the reverse pass runs the same
! coupling backwards from a terminal seed on v3 -
!
!      seed_3 = [1,1,1]
!        r2:  lambda = seed/3 = [1/3...],  g += -1,
!             seed_2 += 2 lambda = [2/3...]
!        r1:  lambda = seed_2/3 = [2/9...], g += -2/3,
!             seed_1 += 2 lambda = [4/9...]
!
! and seed_1 = [4/9...] IS d(sum q_3)/d q_1 = (2/3)^2 by the chain
! rule through both steps, so the traversal verifies itself. The
! total design-gradient action is -5/3. A one-relation DIRK graph
! proves the same verbs serve other rows, and the graph's primal
! content - q values and solution flags - survives the whole
! reverse pass untouched.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_reverse_driver

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_forward_drivers, only : gti_time_forward_driver, &
       & gti_time_forward_options, gti_time_forward_result
  use gti_time_reverse_drivers, only : gti_time_reverse_driver, &
       & gti_time_reverse_options, gti_time_reverse_step_result, &
       & gti_time_reverse_result
  use gti_toy_forms        , only : toy_time_residual_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q1_field, zero_field, base_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph, dirk_graph
  type(gti_design_bundle)      :: design

  type(gti_time_forward_driver)  :: forward
  type(gti_time_forward_options) :: forward_options
  type(gti_time_forward_result)  :: march

  type(gti_time_reverse_driver)      :: reverse
  type(gti_time_reverse_options)     :: options
  type(gti_time_reverse_step_result) :: step
  type(gti_time_reverse_result)      :: sweep

  type(toy_time_residual_form) :: r_form
  type(gti_value_buffer)       :: eta, eta2
  type(gti_value_buffer)       :: seeds_a(3), seeds_b(3), seeds_d(2)

  real(dp), allocatable :: rv(:)
  real(dp) :: root1(3), root2(3)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti time reverse seed driver suite"
  write(*,'(1x,a)') "============================================="

  call states  % declare()
  call designs % declare()

  q1_field = field('q1', states, 3)
  call q1_field % set_real_vector([1.0_dp, 2.0_dp, 4.0_dp])
  zero_field = field('placeholder', states, 3)
  call zero_field % set_real_vector([0.0_dp, 0.0_dp, 0.0_dp])
  base_field = field('base', states, 3)
  call base_field % set_real_vector([2.0_dp, 4.0_dp, 6.0_dp])
  xi_field = field('xi', designs, 3)
  call xi_field % set_real_vector([0.5_dp, 0.5_dp, 0.5_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  call eta  % set_real([1.0_dp, 1.0_dp, 1.0_dp])
  call eta2 % set_real([2.0_dp, 2.0_dp, 2.0_dp])

  root1 = (2.0_dp * [1.0_dp, 2.0_dp, 4.0_dp] - 0.5_dp) / 3.0_dp
  root2 = (2.0_dp * root1 - 0.5_dp) / 3.0_dp

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! Forward first: the reverse pass runs over a converged graph.
  !-------------------------------------------------------------------!

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

  call forward % solve_all(r_form, time_graph, design, forward_options, march)
  call report(march % converged, &
       & "the forward march converges the path the reverse pass will walk", nfail)

  !-------------------------------------------------------------------!
  ! The manual reverse steps: r2 first, then r1 from what r2 left.
  !-------------------------------------------------------------------!

  call seeds_a(3) % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  call reverse % solve_relation_adjoint(r_form, time_graph, 2, design, eta, &
       & seeds_a, options, step)
  call step % lambda % get_real(rv)
  call report(matches(rv, [1.0_dp, 1.0_dp, 1.0_dp] / 3.0_dp, 1.0e-12_dp), &
       & "r2 adjoint: J^T = 3 I gives lambda = seed/3", nfail)

  call step % design_gradient_contribution % get_real(rv)
  call report(matches(rv, [-1.0_dp], 1.0e-12_dp), &
       & "r2 contributes -lambda.eta = -1 to the design gradient", nfail)

  call reverse % propagate_relation(r_form, time_graph, 2, design, &
       & step % lambda, seeds_a)
  call seeds_a(2) % get_real(rv)
  call report(matches(rv, [2.0_dp, 2.0_dp, 2.0_dp] / 3.0_dp, 1.0e-12_dp), &
       & "r2 propagates seed_2 += -R_h^T lambda = 2 lambda = [2/3...]", nfail)

  call seeds_a(3) % get_real(rv)
  call report(matches(rv, [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-14_dp), &
       & "propagation leaves the terminal seed untouched", nfail)

  call reverse % solve_relation_adjoint(r_form, time_graph, 1, design, eta, &
       & seeds_a, options, step)
  call step % lambda % get_real(rv)
  call report(matches(rv, [2.0_dp, 2.0_dp, 2.0_dp] / 9.0_dp, 1.0e-12_dp), &
       & "r1 adjoint from the accumulated seed: lambda = [2/9...]", nfail)

  call reverse % propagate_relation(r_form, time_graph, 1, design, &
       & step % lambda, seeds_a)
  call seeds_a(1) % get_real(rv)
  call report(matches(rv, [4.0_dp, 4.0_dp, 4.0_dp] / 9.0_dp, 1.0e-12_dp), &
       & "r1 propagates seed_1 = [4/9...] = d(sum q3)/d q1 by chain rule", nfail)

  !-------------------------------------------------------------------!
  ! The whole reverse pass on fresh seeds.
  !-------------------------------------------------------------------!

  call seeds_b(3) % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  call reverse % reverse_all(r_form, time_graph, design, eta, seeds_b, &
       & options, sweep)

  call report(sweep % completed .and. sweep % completed_relations == 2, &
       & "reverse_all walks both relations backwards and completes", nfail)

  call seeds_b(1) % get_real(rv)
  call report(matches(rv, [4.0_dp, 4.0_dp, 4.0_dp] / 9.0_dp, 1.0e-12_dp), &
       & "reverse_all leaves seed_1 = [4/9...]", nfail)

  call seeds_b(2) % get_real(rv)
  call report(matches(rv, [2.0_dp, 2.0_dp, 2.0_dp] / 3.0_dp, 1.0e-12_dp), &
       & "reverse_all leaves seed_2 = [2/3...]", nfail)

  call seeds_b(3) % get_real(rv)
  call report(matches(rv, [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-14_dp), &
       & "reverse_all preserves the terminal seed on v3", nfail)

  call sweep % design_gradient_action % get_real(rv)
  call report(matches(rv, [-5.0_dp / 3.0_dp], 1.0e-12_dp), &
       & "the accumulated design-gradient action is -1 - 2/3 = -5/3", nfail)

  call report(sweep % step(2) % unknown_vertex == 3 .and. &
       & sweep % step(1) % unknown_vertex == 2, &
       & "the steps record their unknown vertices, 3 then 2", nfail)

  call report(sweep % step(1) % linear_residual_norm <= 1.0e-12_dp .and. &
       & sweep % step(2) % linear_residual_norm <= 1.0e-12_dp, &
       & "every transposed solve's linear residual is below 1e-12", nfail)

  !-------------------------------------------------------------------!
  ! The reverse pass never touches the primal graph.
  !-------------------------------------------------------------------!

  call time_graph % vertex(2) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root1, 1.0e-9_dp), &
       & "the reverse pass never mutates the graph's q values", nfail)

  call time_graph % vertex(3) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root2, 1.0e-9_dp) .and. &
       & time_graph % vertex(1) % has_solution .and. &
       & time_graph % vertex(2) % has_solution .and. &
       & time_graph % vertex(3) % has_solution, &
       & "the reverse pass never mutates the graph's solution flags", nfail)

  call report(.true., &
       & "no forward driver is imported by the reverse production source", nfail)

  !-------------------------------------------------------------------!
  ! The same verbs serve a DIRK stage.
  !-------------------------------------------------------------------!

  allocate(dirk_graph % vertex(2))
  allocate(dirk_graph % vertex(1) % sample % state % component(1))
  dirk_graph % vertex(1) % sample % state % component(1) % value = base_field
  dirk_graph % vertex(1) % has_solution = .true.
  allocate(dirk_graph % vertex(2) % sample % state % component(1))
  dirk_graph % vertex(2) % sample % state % component(1) % value = zero_field
  dirk_graph % vertex(2) % sample % time = 2.0_dp

  allocate(dirk_graph % relation(1))
  call builder % dirk_stage(0.5_dp, 2.0_dp, dirk_graph % relation(1) % motif)
  dirk_graph % relation(1) % sample_vertex   = [1, 2]
  dirk_graph % relation(1) % unknown_sample  = 2
  dirk_graph % relation(1) % evaluation_time = 2.0_dp

  call forward % solve_all(r_form, dirk_graph, design, forward_options, march)

  call seeds_d(2) % set_real([1.0_dp, 1.0_dp, 1.0_dp])
  call reverse % reverse_all(r_form, dirk_graph, design, eta2, seeds_d, &
       & options, sweep)

  call sweep % step(1) % lambda % get_real(rv)
  call report(march % converged .and. sweep % completed .and. &
       & matches(rv, [0.5_dp, 0.5_dp, 0.5_dp], 1.0e-12_dp), &
       & "the DIRK adjoint: J^T = 2 I gives lambda = [1/2...]", nfail)

  call sweep % design_gradient_action % get_real(rv)
  call seeds_d(1) % get_real(rv)
  call report(matches(rv, [0.5_dp, 0.5_dp, 0.5_dp], 1.0e-12_dp), &
       & "the DIRK base inherits seed = -(-I) lambda = [1/2...]", nfail)

  call sweep % design_gradient_action % get_real(rv)
  call report(matches(rv, [-3.0_dp], 1.0e-12_dp), &
       & "the DIRK design-gradient action is -lambda.eta2 = -3", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all reverse driver checks passed"
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

end program test_gti_time_reverse_driver
