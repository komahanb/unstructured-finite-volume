!=====================================================================!
! The functional seed driver suite: scalar terms on time vertices
! become the three pieces the reverse pass needs -
!
!      F_time = sum_k F_k(q_{v_k}, xi, t_k)
!        -> scalar value
!        -> vertex q-seeds
!        -> explicit F_xi[eta]
!
! over the solved two-relation BDF1 path with the clocked sum
! functional F = sum q + sum xi + t. The terminal term at v3
! gives value 41/18 + 3/2 + 1 = 43/9, seed [1,1,1], and explicit
! design action 3.
!
! (The task sheet's expected value reads 34/9, which omits the +t
! term its own formula and time-law include; the true terminal
! value at t = 1 is 43/9, used here.)
!
! The crown is the composition this seat was built for: seed_all
! then the EXISTING reverse driver, and the total design action
!
!      3 + (-5/3) = 4/3
!
! equals d(sum q3 + sum xi)/dxi[eta] by hand chain rule - the full
! gradient of a time functional over the marched graph, assembled
! from parts that never saw each other.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_functional_seed_driver

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_design_bundles   , only : gti_design_bundle
  use gti_time_local_schemes, only : gti_evaluation_point
  use gti_time_motif_builders, only : gti_time_motif_builder
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_forward_drivers, only : gti_time_forward_driver, &
       & gti_time_forward_options, gti_time_forward_result
  use gti_time_reverse_drivers, only : gti_time_reverse_driver, &
       & gti_time_reverse_options, gti_time_reverse_result
  use gti_time_functional_seed_drivers, only : gti_time_functional_term, &
       & gti_time_functional, gti_time_functional_seed_result, &
       & gti_time_functional_seed_driver
  use gti_toy_forms        , only : toy_time_residual_form, &
       & toy_sum_time_functional

  implicit none

  type(graph) :: states, designs
  type(field) :: q1_field, zero_field, xi_field

  type(gti_time_motif_builder) :: builder
  type(gti_time_graph)         :: time_graph
  type(gti_design_bundle)      :: design

  type(gti_time_forward_driver)  :: forward
  type(gti_time_forward_options) :: forward_options
  type(gti_time_forward_result)  :: march

  type(gti_time_reverse_driver)  :: reverse
  type(gti_time_reverse_options) :: reverse_options
  type(gti_time_reverse_result)  :: sweep

  type(gti_time_functional)             :: terminal, spread_terms
  type(gti_time_functional_seed_driver) :: seeder
  type(gti_time_functional_seed_result) :: seeding

  type(toy_time_residual_form)   :: r_form
  type(toy_sum_time_functional)  :: f_form

  type(gti_evaluation_point) :: point
  type(gti_value_buffer)     :: eta, out
  type(gti_value_buffer)     :: seeds_a(3), seeds_b(3), seeds_c(3)

  real(dp), allocatable :: rv(:), total_contrib(:)
  real(dp) :: root1(3), root2(3), value3, value2
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti time functional seed driver suite"
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

  root1 = (2.0_dp * [1.0_dp, 2.0_dp, 4.0_dp] - 0.5_dp) / 3.0_dp
  root2 = (2.0_dp * root1 - 0.5_dp) / 3.0_dp

  value3 = sum(root2) + 1.5_dp + 1.0_dp
  value2 = sum(root1) + 1.5_dp + 0.5_dp

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! Forward first: seeds are built on a solved trajectory.
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
       & "the forward march solves the trajectory the terms will read", nfail)

  !-------------------------------------------------------------------!
  ! The terminal functional: one term on v3 at t = 1.
  !-------------------------------------------------------------------!

  allocate(terminal % term(1))
  terminal % term(1) % vertex_index    = 3
  terminal % term(1) % evaluation_time = 1.0_dp

  call report(terminal % num_terms() == 1, &
       & "the terminal functional counts one term", nfail)

  call terminal % validate(time_graph)
  call report(.true., &
       & "validate accepts the terminal functional on the solved graph", nfail)

  call seeder % build_point(time_graph, terminal % term(1), design, point)
  call point % state % component(1 + GTI_STATE_Q) % value % get_real_vector(rv)
  call report(abs(point % time - 1.0_dp) < 1.0e-14_dp .and. &
       & matches(rv, root2, 1.0e-9_dp), &
       & "build_point carries the vertex state and the term time", nfail)

  call seeder % value(f_form, time_graph, terminal % term(1), design, out)
  call out % get_real(rv)
  call report(matches(rv, [value3], 1.0e-12_dp), &
       & "the terminal value is sum(q3) + sum(xi) + t = 43/9", nfail)

  call seeder % q_gradient(f_form, time_graph, terminal % term(1), design, out)
  call out % get_real(rv)
  call report(matches(rv, [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-12_dp), &
       & "the q-gradient of the sum functional is [1, 1, 1]", nfail)

  call seeder % design_action(f_form, time_graph, terminal % term(1), design, &
       & eta, out)
  call out % get_real(rv)
  call report(matches(rv, [3.0_dp], 1.0e-12_dp), &
       & "the explicit design action is F_xi[eta] = 3", nfail)

  !-------------------------------------------------------------------!
  ! seed_all: value gathered, seed landed, explicit part summed.
  !-------------------------------------------------------------------!

  call seeder % seed_all(f_form, time_graph, terminal, design, eta, &
       & seeds_a, seeding)

  call report(seeding % completed, &
       & "seed_all completes", nfail)

  call report(abs(seeding % value - value3) < 1.0e-12_dp, &
       & "seed_all gathers the scalar value 43/9", nfail)

  call seeds_a(3) % get_real(rv)
  call report(matches(rv, [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-12_dp), &
       & "seed_all lands the terminal seed [1, 1, 1] on v3", nfail)

  call seeds_a(1) % get_real(rv)
  call report(size(rv) == 0, &
       & "seed_all leaves the unseeded v1 empty", nfail)

  call seeds_a(2) % get_real(rv)
  call report(size(rv) == 0, &
       & "seed_all leaves the unseeded v2 empty", nfail)

  call seeding % explicit_design_action % get_real(rv)
  call report(matches(rv, [3.0_dp], 1.0e-12_dp), &
       & "seed_all sums the explicit design action to [3]", nfail)

  call time_graph % vertex(3) % sample % state % component(1 + GTI_STATE_Q) % &
       & value % get_real_vector(rv)
  call report(matches(rv, root2, 1.0e-9_dp), &
       & "seeding never mutates the graph's q values", nfail)

  call report(time_graph % vertex(1) % has_solution .and. &
       & time_graph % vertex(2) % has_solution .and. &
       & time_graph % vertex(3) % has_solution, &
       & "seeding never mutates the graph's solution flags", nfail)

  !-------------------------------------------------------------------!
  ! Accumulation: two terms on v3 and one on v2 - duplicates add,
  ! vertices seed independently, values and design actions sum.
  !-------------------------------------------------------------------!

  allocate(spread_terms % term(3))
  spread_terms % term(1) % vertex_index    = 3
  spread_terms % term(1) % evaluation_time = 1.0_dp
  spread_terms % term(2) % vertex_index    = 3
  spread_terms % term(2) % evaluation_time = 1.0_dp
  spread_terms % term(3) % vertex_index    = 2
  spread_terms % term(3) % evaluation_time = 0.5_dp

  call seeder % seed_all(f_form, time_graph, spread_terms, design, eta, &
       & seeds_b, seeding)

  call seeds_b(3) % get_real(rv)
  call report(matches(rv, [2.0_dp, 2.0_dp, 2.0_dp], 1.0e-12_dp), &
       & "duplicate terms on one vertex add: v3 seed is [2, 2, 2]", nfail)

  call seeds_b(2) % get_real(rv)
  call report(matches(rv, [1.0_dp, 1.0_dp, 1.0_dp], 1.0e-12_dp), &
       & "a term on v2 seeds v2 independently", nfail)

  call seeding % explicit_design_action % get_real(rv)
  call report(matches(rv, [9.0_dp], 1.0e-12_dp), &
       & "explicit design actions add across terms: [9]", nfail)

  call report(abs(seeding % value - (2.0_dp * value3 + value2)) < 1.0e-12_dp, &
       & "values add across terms: 2 F(v3) + F(v2)", nfail)

  !-------------------------------------------------------------------!
  ! The composition this seat exists for: seed, then reverse, and
  ! the total design action equals hand calculus.
  !-------------------------------------------------------------------!

  call seeder % seed_all(f_form, time_graph, terminal, design, eta, &
       & seeds_c, seeding)

  call reverse % reverse_all(r_form, time_graph, design, eta, seeds_c, &
       & reverse_options, sweep)

  call seeding % explicit_design_action % get_real(rv)
  call sweep % design_gradient_action % get_real(total_contrib)
  call report(sweep % completed .and. &
       & abs(rv(1) + total_contrib(1) - 4.0_dp / 3.0_dp) < 1.0e-12_dp, &
       & "seed_all + reverse_all: total dF/dxi[eta] = 3 - 5/3 = 4/3", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all functional seed checks passed"
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

end program test_gti_time_functional_seed_driver
