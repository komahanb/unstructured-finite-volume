!=====================================================================!
! The GTI form suite (phase 1): one residual form and one
! functional form answer the same abstract contract -
!
!      name, input_signature, output_signature, max_degree,
!      value, partial_action
!
! - and their partial actions to order 2 reproduce the exact
! calculus of the toy pair
!
!      R_i(q, xi) = q_i^2 + xi
!      F(q, xi)   = 1/2 sum q_i^2 + xi
!
! with the bundles reporting shape and occupancy separately, and
! the value buffer obeying the one shape law of the library.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_form_interface

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_state_bundles    , only : gti_state_bundle
  use gti_design_bundles   , only : gti_design_bundle
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, &
       & GTI_ARG_STATE, GTI_ARG_DESIGN, &
       & GTI_STATE_Q, GTI_STATE_QDOT, GTI_STATE_QDDOT
  use gti_toy_forms        , only : toy_residual_form, toy_functional_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q_field, qdot_field, xi_field

  type(gti_evaluation_point) :: point
  type(gti_state_bundle)     :: loose_state
  type(gti_design_bundle)    :: loose_design, seated_design

  type(toy_residual_form)   :: r_form
  type(toy_functional_form) :: f_form

  type(gti_value_buffer)     :: buffer, out, out0
  type(gti_partial_request)  :: request
  type(gti_direction_bundle) :: dir_q_v1, dir_q_v2, dir_qdot, dir_xi
  type(gti_direction_bundle) :: no_directions(0)

  real(dp), allocatable :: rv(:), rv0(:)
  integer , allocatable :: sig(:)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti differentiable form suite (phase 1)"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! The value buffer: one shape law, absence answered by zero
  ! length.
  !-------------------------------------------------------------------!

  call buffer % get_real(rv)
  call report(size(rv) == 0 .and. buffer % nentries == 0, &
       & "a fresh buffer answers with nothing, zero length", nfail)

  call buffer % set_real([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], ncomp=2)
  call buffer % get_real(rv)
  call report(buffer % nentries == 3 .and. buffer % ncomp == 2 .and. &
       & matches(rv, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]), &
       & "set_real replaces values and shape together", nfail)

  call buffer % clear()
  call buffer % get_real(rv)
  call report(size(rv) == 0 .and. buffer % nentries == 0 .and. buffer % ncomp == 1, &
       & "clear forgets the values and empties the shape", nfail)

  !-------------------------------------------------------------------!
  ! The evaluation point: q = [1,2,3] and qdot = [4,5,6] on the
  ! state set, xi = [0.5] on the design set, and one empty seat
  ! for qddot.
  !-------------------------------------------------------------------!

  call states  % declare()
  call designs % declare()

  q_field = field('q', states, 3)
  call q_field % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp])

  qdot_field = field('qdot', states, 3)
  call qdot_field % set_real_vector([4.0_dp, 5.0_dp, 6.0_dp])

  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([0.5_dp])

  allocate(point % state % component(3))
  point % state % component(1) % value = q_field
  point % state % component(2) % value = qdot_field

  allocate(point % design % component(1))
  point % design % component(1) % value = xi_field

  !-------------------------------------------------------------------!
  ! The state bundle: shape and occupancy are separate answers.
  !-------------------------------------------------------------------!

  call report(loose_state % differential_degree() == -1 .and. &
       & .not. loose_state % has_component(GTI_STATE_Q), &
       & "an unallocated state bundle has degree -1 and no q", nfail)

  call report(point % state % differential_degree() == 2, &
       & "three seats make degree 2, occupied or not", nfail)

  call report(point % state % has_component(GTI_STATE_Q) .and. &
       &      point % state % has_component(GTI_STATE_QDOT) .and. &
       & .not. point % state % has_component(GTI_STATE_QDDOT), &
       & "q and qdot sit, qddot's seat stays empty", nfail)

  !-------------------------------------------------------------------!
  ! The design bundle: a seat is not an entry.
  !-------------------------------------------------------------------!

  call report(loose_design % size() == 0 .and. .not. loose_design % has_entries(), &
       & "an unallocated design bundle has size 0 and no entries", nfail)

  allocate(seated_design % component(1))
  call report(seated_design % size() == 1 .and. .not. seated_design % has_entries(), &
       & "a seat is not an entry: size 1, no entries", nfail)

  call report(point % design % size() == 1 .and. point % design % has_entries(), &
       & "one design field seated: size 1, entries present", nfail)

  !-------------------------------------------------------------------!
  ! One abstract contract serves both forms: the same class dummy
  ! answers name, signatures, and degree for R and for F.
  !-------------------------------------------------------------------!

  r_form % nequations = 3

  call check_shared_contract(r_form, 'toy residual',   [3, 1], nfail)
  call check_shared_contract(f_form, 'toy functional', [1, 1], nfail)

  !-------------------------------------------------------------------!
  ! Values: residual-sized for R, scalar for F, each filling the
  ! shape its output_signature declares.
  !-------------------------------------------------------------------!

  call r_form % value(point, out)
  call out % get_real(rv)
  call report(matches(rv, [1.5_dp, 4.5_dp, 9.5_dp]), &
       & "R value: R_i = q_i^2 + xi, residual-sized", nfail)

  sig = r_form % output_signature()
  call report(out % nentries == sig(1) .and. out % ncomp == sig(2), &
       & "R fills the shape its output_signature declares", nfail)

  call f_form % value(point, out)
  call out % get_real(rv)
  call report(matches(rv, [7.5_dp]), &
       & "F value: F = 1/2 sum q^2 + xi, scalar", nfail)

  sig = f_form % output_signature()
  call report(out % nentries == sig(1) .and. out % ncomp == sig(2), &
       & "F fills the shape its output_signature declares", nfail)

  !-------------------------------------------------------------------!
  ! Order 0: the empty contraction is the value itself.
  !-------------------------------------------------------------------!

  request = gti_partial_request()

  call r_form % value(point, out0)
  call out0 % get_real(rv0)
  call r_form % partial_action(point, request, no_directions, out)
  call out % get_real(rv)
  call report(matches(rv, rv0), &
       & "order 0 partial action is the value itself (R)", nfail)

  call f_form % value(point, out0)
  call out0 % get_real(rv0)
  call f_form % partial_action(point, request, no_directions, out)
  call out % get_real(rv)
  call report(matches(rv, rv0), &
       & "order 0 partial action is the value itself (F)", nfail)

  !-------------------------------------------------------------------!
  ! Order 1: directional derivatives along v = [1,0,2] in q and
  ! w = [0.25] in xi. A direction names the argument it perturbs.
  !-------------------------------------------------------------------!

  call dir_q_v1 % values % set_real([1.0_dp, 0.0_dp, 2.0_dp])

  dir_xi % argument_kind = GTI_ARG_DESIGN
  call dir_xi % values % set_real([0.25_dp])

  dir_qdot % state_component = GTI_STATE_QDOT
  call dir_qdot % values % set_real([1.0_dp, 1.0_dp, 1.0_dp])

  request = gti_partial_request(order=1, argument_kind=[GTI_ARG_STATE], &
       & state_component=[GTI_STATE_Q])
  call r_form % partial_action(point, request, [dir_q_v1], out)
  call out % get_real(rv)
  call report(matches(rv, [2.0_dp, 0.0_dp, 12.0_dp]), &
       & "order 1 state action on R: D R [v] = 2 q v", nfail)

  request = gti_partial_request(order=1, argument_kind=[GTI_ARG_DESIGN])
  call r_form % partial_action(point, request, [dir_xi], out)
  call out % get_real(rv)
  call report(matches(rv, [0.25_dp, 0.25_dp, 0.25_dp]), &
       & "order 1 design action on R: D R [w] = w on every equation", nfail)

  request = gti_partial_request(order=1, argument_kind=[GTI_ARG_STATE], &
       & state_component=[GTI_STATE_Q])
  call f_form % partial_action(point, request, [dir_q_v1], out)
  call out % get_real(rv)
  call report(matches(rv, [7.0_dp]), &
       & "order 1 state action on F: D F [v] = sum q v", nfail)

  request = gti_partial_request(order=1, argument_kind=[GTI_ARG_DESIGN])
  call f_form % partial_action(point, request, [dir_xi], out)
  call out % get_real(rv)
  call report(matches(rv, [0.25_dp]), &
       & "order 1 design action on F: D F [w] = w", nfail)

  request = gti_partial_request(order=1, argument_kind=[GTI_ARG_STATE], &
       & state_component=[GTI_STATE_QDOT])
  call f_form % partial_action(point, request, [dir_qdot], out)
  call out % get_real(rv)
  call report(matches(rv, [0.0_dp]), &
       & "a qdot direction contracts F to zero: F never reads qdot", nfail)

  !-------------------------------------------------------------------!
  ! Order 2: the toy curvature lives in (q, q) alone; every mixed
  ! slot contracts to zero.
  !-------------------------------------------------------------------!

  call dir_q_v2 % values % set_real([3.0_dp, 1.0_dp, 1.0_dp])

  request = gti_partial_request(order=2, &
       & argument_kind=[GTI_ARG_STATE, GTI_ARG_STATE], &
       & state_component=[GTI_STATE_Q, GTI_STATE_Q])
  call r_form % partial_action(point, request, [dir_q_v1, dir_q_v2], out)
  call out % get_real(rv)
  call report(matches(rv, [6.0_dp, 0.0_dp, 4.0_dp]), &
       & "order 2 action on R: D2 R [v1,v2] = 2 v1 v2", nfail)

  request = gti_partial_request(order=2, &
       & argument_kind=[GTI_ARG_STATE, GTI_ARG_DESIGN], &
       & state_component=[GTI_STATE_Q, GTI_STATE_Q])
  call r_form % partial_action(point, request, [dir_q_v1, dir_xi], out)
  call out % get_real(rv)
  call report(matches(rv, [0.0_dp, 0.0_dp, 0.0_dp]), &
       & "order 2 mixed state-design action on R is zero", nfail)

  request = gti_partial_request(order=2, &
       & argument_kind=[GTI_ARG_STATE, GTI_ARG_STATE], &
       & state_component=[GTI_STATE_Q, GTI_STATE_Q])
  call f_form % partial_action(point, request, [dir_q_v1, dir_q_v2], out)
  call out % get_real(rv)
  call report(matches(rv, [5.0_dp]), &
       & "order 2 action on F: D2 F [v1,v2] = sum v1 v2", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all gti form checks passed"
  else
     error stop
  end if

contains

  !-------------------------------------------------------------------!
  ! The shared-contract probe: one class(gti_differentiable_form)
  ! dummy answers for whichever form arrives - that single dummy
  ! IS the proof that R and F share the abstract interface.
  !-------------------------------------------------------------------!

  subroutine check_shared_contract(phi, expected_name, expected_out, nfail)

    class(gti_differentiable_form), intent(in) :: phi
    character(len=*), intent(in)    :: expected_name
    integer         , intent(in)    :: expected_out(2)
    integer         , intent(inout) :: nfail

    integer, allocatable :: insig(:), outsig(:)

    insig  = phi % input_signature()
    outsig = phi % output_signature()

    call report(phi % name() == expected_name .and. &
         & size(insig) == 2 .and. &
         & all(insig == [GTI_ARG_STATE, GTI_ARG_DESIGN]) .and. &
         & size(outsig) == 2 .and. &
         & all(outsig == expected_out) .and. &
         & phi % max_degree() == 2, &
         & expected_name // " answers the one abstract contract", nfail)

  end subroutine check_shared_contract

  pure function matches(values, expected) result(ok)

    real(dp), intent(in) :: values(:), expected(:)
    logical :: ok

    ok = size(values) == size(expected)
    if (ok) ok = all(abs(values - expected) < 1.0e-14_dp)

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

end program test_gti_form_interface
