!=====================================================================!
! The adaptive controller suite: the policy layer's first
! inhabitant, proven over the scalar nonlinear chain
!
!      R = qdot + q^2 - xi,   xi = 1,   q0 = 2,
!
! whose BDF1 step has the closed root
!
!      q_u(q_h, h) = ( -1 + sqrt(1 + 4 h (q_h + h)) ) / (2 h).
!
! Every decision below is the policy's and every execution is the
! mechanisms': a generous tolerance accepts the first try; a tuned
! tolerance rejects h and accepts h/2, the rejected step leaving
! no trace; a rootless residual spends the whole budget and fails
! LAWFULLY with the graph exactly as it was; the halving policy
! clamps its preferred order to the history and promotes to BDF2 -
! the variable-step rows pricing genuinely nonuniform spacing -
! the moment two vertices exist; and three advances in a row
! march the chain, proving adaptivity is composition. Not one
! line of any traversal or growth driver changed to make this
! seat exist.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_time_adaptive_controller

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_design_bundles   , only : gti_design_bundle
  use gti_state_bundles    , only : GTI_STATE_Q
  use gti_time_graphs      , only : gti_time_graph
  use gti_time_forward_drivers, only : gti_time_forward_options
  use gti_time_adaptive_controllers, only : &
       & gti_time_adaptive_controller, gti_time_adaptive_controller_options, &
       & gti_time_adaptive_advance_result, gti_time_halving_step_policy
  use gti_toy_forms        , only : toy_qdot_square_form, toy_square_plus_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q0_field, one_field, xi_field

  type(gti_time_graph)    :: time_graph
  type(gti_design_bundle) :: design

  type(gti_time_adaptive_controller)         :: controller
  type(gti_time_adaptive_controller_options) :: options
  type(gti_time_forward_options)             :: forward_options
  type(gti_time_halving_step_policy)         :: policy
  type(gti_time_adaptive_advance_result)     :: result

  type(toy_qdot_square_form) :: r_form
  type(toy_square_plus_form) :: rootless

  real(dp), allocatable :: rv(:)
  real(dp) :: expected1, expected2, expected3
  integer  :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti adaptive controller suite"
  write(*,'(1x,a)') "============================================="

  call states  % declare()
  call designs % declare()

  q0_field = field('q0', states, 1)
  call q0_field % set_real_vector([2.0_dp])
  one_field = field('one', states, 1)
  call one_field % set_real_vector([1.0_dp])
  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([1.0_dp])

  allocate(design % component(1))
  design % component(1) % value = xi_field

  !-------------------------------------------------------------------!
  ! A generous tolerance: the first proposal is accepted, the
  ! graph grows by one vertex, and both the root and the measured
  ! evidence match their closed forms.
  !-------------------------------------------------------------------!

  call seed_base_graph(time_graph, q0_field, 0.0_dp)

  policy % initial_step    = 0.5_dp
  policy % error_tolerance = 1.0e9_dp
  policy % preferred_order = 1

  call controller % advance(r_form, time_graph, design, policy, &
       & forward_options, options, result)

  expected1 = bdf1_root(2.0_dp, 0.5_dp)

  call report(result % accepted .and. result % attempts == 1 .and. &
       & result % appended_vertex == 2, &
       & "a generous tolerance accepts the first proposal", nfail)

  call read_vertex_q(time_graph, 2, rv)
  call report(time_graph % num_vertices() == 2 .and. &
       & matches(rv, [expected1], 1.0e-9_dp) .and. &
       & matches([time_graph % vertex(2) % sample % time], [0.5_dp], 1.0e-14_dp), &
       & "the accepted step stands: t = 1/2, q = the closed BDF1 root", nfail)

  call report(result % attempt(1) % kept .and. &
       & result % attempt(1) % evidence % converged .and. &
       & result % attempt(1) % evidence % has_estimate .and. &
       & matches([result % attempt(1) % evidence % error_estimate], &
       &         [abs(expected1 - 2.0_dp)], 1.0e-9_dp), &
       & "the evidence is the distance from the constant predictor", nfail)

  !-------------------------------------------------------------------!
  ! A tuned tolerance: h = 1/2 is rejected, h = 1/4 is accepted,
  ! and the rejected step leaves no trace - only the accepted
  ! vertex exists, at the accepted time.
  !-------------------------------------------------------------------!

  call seed_base_graph(time_graph, q0_field, 0.0_dp)

  policy % initial_step    = 0.5_dp
  policy % error_tolerance = 0.45_dp

  call controller % advance(r_form, time_graph, design, policy, &
       & forward_options, options, result)

  expected2 = bdf1_root(2.0_dp, 0.25_dp)

  call report(result % accepted .and. result % attempts == 2 .and. &
       & .not. result % attempt(1) % kept .and. result % attempt(2) % kept, &
       & "the tuned tolerance rejects h and accepts h/2", nfail)

  call read_vertex_q(time_graph, 2, rv)
  call report(time_graph % num_vertices() == 2 .and. &
       & matches(rv, [expected2], 1.0e-9_dp) .and. &
       & matches([time_graph % vertex(2) % sample % time], [0.25_dp], 1.0e-14_dp), &
       & "the rejected step left no trace: one new vertex, at t = 1/4", nfail)

  call report(matches([result % attempt(1) % evidence % error_estimate], &
       &              [abs(bdf1_root(2.0_dp, 0.5_dp) - 2.0_dp)], 1.0e-9_dp) .and. &
       & matches([result % attempt(1) % evidence % step], [0.5_dp], 1.0e-14_dp) .and. &
       & matches([result % attempt(2) % evidence % step], [0.25_dp], 1.0e-14_dp), &
       & "both attempts' evidence is on the record, exact", nfail)

  !-------------------------------------------------------------------!
  ! Order promotion: on one history vertex the policy's preferred
  ! order 2 is clamped to 1; the moment the graph carries two, the
  ! next advance builds the BDF2 candidate.
  !-------------------------------------------------------------------!

  call seed_base_graph(time_graph, q0_field, 0.0_dp)

  policy % initial_step    = 0.5_dp
  policy % error_tolerance = 1.0e9_dp
  policy % preferred_order = 2

  call controller % advance(r_form, time_graph, design, policy, &
       & forward_options, options, result)
  call report(result % accepted .and. &
       & result % attempt(1) % evidence % order == 1, &
       & "one history vertex clamps the preferred order to BDF1", nfail)

  call controller % advance(r_form, time_graph, design, policy, &
       & forward_options, options, result)

  expected3 = bdf2_uniform_root(2.0_dp, expected1)

  call read_vertex_q(time_graph, 3, rv)
  call report(result % accepted .and. &
       & result % attempt(1) % evidence % order == 2 .and. &
       & matches(rv, [expected3], 1.0e-9_dp), &
       & "two history vertices promote to BDF2, and its root is exact", nfail)

  !-------------------------------------------------------------------!
  ! Nonuniform order 2 at equilibrium: history spaced by h = 1,
  ! candidate step h = 2/5 - the variable-step rows price the
  ! spacing, the equilibrium root is exactly one, and the linear
  ! predictor makes the estimate exactly zero, accepted under the
  ! tightest of tolerances.
  !-------------------------------------------------------------------!

  call seed_two_vertex_graph(time_graph, one_field)

  policy % initial_step    = 0.4_dp
  policy % error_tolerance = 1.0e-12_dp
  policy % preferred_order = 2

  call controller % advance(r_form, time_graph, design, policy, &
       & forward_options, options, result)

  call read_vertex_q(time_graph, 3, rv)
  call report(result % accepted .and. &
       & result % attempt(1) % evidence % order == 2 .and. &
       & matches(rv, [1.0_dp], 1.0e-12_dp) .and. &
       & matches([time_graph % vertex(3) % sample % time], [1.4_dp], 1.0e-14_dp) .and. &
       & result % attempt(1) % evidence % error_estimate <= 1.0e-12_dp, &
       & "nonuniform BDF2 holds the equilibrium under a tight tolerance", nfail)

  !-------------------------------------------------------------------!
  ! A rootless residual: every attempt fails to converge, the
  ! budget is spent, and failure is a lawful answer - the graph
  ! stands exactly as it was, and nothing error-stopped.
  !-------------------------------------------------------------------!

  call seed_base_graph(time_graph, q0_field, 0.0_dp)

  policy % initial_step    = 0.5_dp
  policy % error_tolerance = 1.0e9_dp
  policy % preferred_order = 1
  options % max_attempts   = 3

  call controller % advance(rootless, time_graph, design, policy, &
       & forward_options, options, result)

  call read_vertex_q(time_graph, 1, rv)
  call report(.not. result % accepted .and. result % attempts == 3 .and. &
       & time_graph % num_vertices() == 1 .and. &
       & matches(rv, [2.0_dp], 1.0e-14_dp), &
       & "a spent budget is a lawful answer: the graph stands untouched", nfail)

  call report(.not. result % attempt(1) % evidence % converged .and. &
       & .not. result % attempt(3) % evidence % converged .and. &
       & .not. result % attempt(1) % evidence % has_estimate .and. &
       & .not. result % attempt(2) % kept, &
       & "every failed attempt is on the record: unconverged, unmeasured", nfail)

  options % max_attempts = 8

  !-------------------------------------------------------------------!
  ! Marching is composition: three advances in a row walk the
  ! BDF1 chain through its closed roots.
  !-------------------------------------------------------------------!

  call seed_base_graph(time_graph, q0_field, 0.0_dp)

  policy % initial_step    = 0.5_dp
  policy % error_tolerance = 1.0e9_dp
  policy % preferred_order = 1

  call controller % advance(r_form, time_graph, design, policy, &
       & forward_options, options, result)
  call controller % advance(r_form, time_graph, design, policy, &
       & forward_options, options, result)
  call controller % advance(r_form, time_graph, design, policy, &
       & forward_options, options, result)

  expected1 = bdf1_root(2.0_dp, 0.5_dp)
  expected2 = bdf1_root(expected1, 0.5_dp)
  expected3 = bdf1_root(expected2, 0.5_dp)

  call read_vertex_q(time_graph, 4, rv)
  call report(time_graph % num_vertices() == 4 .and. &
       & time_graph % num_relations() == 3 .and. &
       & matches(rv, [expected3], 1.0e-9_dp) .and. &
       & matches([time_graph % vertex(4) % sample % time], [1.5_dp], 1.0e-14_dp), &
       & "three advances march the chain through its closed roots", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all adaptive controller checks passed"
  else
     error stop
  end if

contains

  subroutine seed_base_graph(g, q_field, t0)

    type(gti_time_graph), intent(out) :: g
    type(field)         , intent(in)  :: q_field
    real(dp)            , intent(in)  :: t0

    allocate(g % vertex(1))
    allocate(g % vertex(1) % sample % state % component(1))
    g % vertex(1) % sample % state % component(1) % value = q_field
    g % vertex(1) % sample % time = t0
    g % vertex(1) % has_solution = .true.

  end subroutine seed_base_graph

  subroutine seed_two_vertex_graph(g, q_field)

    type(gti_time_graph), intent(out) :: g

    type(field), intent(in) :: q_field
    integer :: v

    allocate(g % vertex(2))
    do v = 1, 2
       allocate(g % vertex(v) % sample % state % component(1))
       g % vertex(v) % sample % state % component(1) % value = q_field
       g % vertex(v) % sample % time = real(v - 1, dp)
       g % vertex(v) % has_solution = .true.
    end do

  end subroutine seed_two_vertex_graph

  subroutine read_vertex_q(g, vertex_index, values)

    type(gti_time_graph) , intent(in)  :: g
    integer              , intent(in)  :: vertex_index
    real(dp), allocatable, intent(out) :: values(:)

    call g % vertex(vertex_index) % sample % state % &
         & component(1 + GTI_STATE_Q) % value % get_real_vector(values)

  end subroutine read_vertex_q

  !-------------------------------------------------------------------!
  ! The closed roots: BDF1 from history q_h over h, and uniform
  ! BDF2 from [q0, q1] over h = 1/2 - both at xi = 1, both by the
  ! quadratic formula, independent of every seat under test.
  !-------------------------------------------------------------------!

  pure function bdf1_root(q_h, h) result(root)

    real(dp), intent(in) :: q_h, h
    real(dp) :: root

    root = (-1.0_dp + sqrt(1.0_dp + 4.0_dp * h * (q_h + h))) / (2.0_dp * h)

  end function bdf1_root

  pure function bdf2_uniform_root(q0, q1) result(root)

    real(dp), intent(in) :: q0, q1
    real(dp) :: root

    ! h = 1/2: qdot = q0 - 4 q1 + 3 q_u, so
    ! q_u^2 + 3 q_u - (1 + 4 q1 - q0) = 0
    root = (-3.0_dp + sqrt(9.0_dp + 4.0_dp * (1.0_dp + 4.0_dp * q1 - q0))) / 2.0_dp

  end function bdf2_uniform_root

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

end program test_gti_time_adaptive_controller
