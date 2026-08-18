!=====================================================================!
! The GTI evaluator suite: every call to R and F goes through the
! stateless gti_form_evaluator, whose form dummy is
! class(gti_differentiable_form) - the test never touches a form's
! bindings directly, so what passes here passes for any form.
!
! The toy pair under the driver:
!
!      R_i(q, xi) = q_i^2 + xi
!      F(q, xi)   = 1/2 sum q_i^2 + xi
!
! and every answer - value, order 0, 1, 2 - is held to the output
! law [nentries, ncomp] the form declares: residual-sized [3,1]
! for R, scalar [1,1] for F.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_gti_form_evaluator

  use iso_fortran_env      , only : dp => REAL64
  use fractal_graph        , only : graph
  use class_graph_field    , only : field
  use gti_value_buffers    , only : gti_value_buffer
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_interface   , only : gti_partial_request, gti_direction_bundle, &
       & GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_STATE_Q
  use gti_form_evaluators  , only : gti_form_evaluator
  use gti_toy_forms        , only : toy_residual_form, toy_functional_form

  implicit none

  type(graph) :: states, designs
  type(field) :: q_field, xi_field

  type(gti_evaluation_point) :: point
  type(gti_form_evaluator)   :: evaluator

  type(toy_residual_form)   :: r_form
  type(toy_functional_form) :: f_form

  type(gti_value_buffer)     :: out, out0
  type(gti_partial_request)  :: request
  type(gti_direction_bundle) :: dir_q_v1, dir_q_v2
  type(gti_direction_bundle) :: no_directions(0)

  real(dp), allocatable :: rv(:), rv0(:)
  integer :: nfail

  nfail = 0
  write(*,'(1x,a)') "============================================="
  write(*,'(1x,a)') "gti form evaluator suite (driver seat)"
  write(*,'(1x,a)') "============================================="

  !-------------------------------------------------------------------!
  ! The point: q = [1,2,3] on the state set, xi = [0.5] on the
  ! design set.
  !-------------------------------------------------------------------!

  call states  % declare()
  call designs % declare()

  q_field = field('q', states, 3)
  call q_field % set_real_vector([1.0_dp, 2.0_dp, 3.0_dp])

  xi_field = field('xi', designs, 1)
  call xi_field % set_real_vector([0.5_dp])

  allocate(point % state % component(1))
  point % state % component(1) % value = q_field

  allocate(point % design % component(1))
  point % design % component(1) % value = xi_field

  r_form % nequations = 3

  !-------------------------------------------------------------------!
  ! Values through the driver.
  !-------------------------------------------------------------------!

  call evaluator % value(r_form, point, out)
  call out % get_real(rv)
  call report(matches(rv, [1.5_dp, 4.5_dp, 9.5_dp]), &
       & "the evaluator answers R value: R_i = q_i^2 + xi", nfail)

  call report(out % nentries == 3 .and. out % ncomp == 1, &
       & "the evaluator admits the residual shape [n,1]", nfail)

  call evaluator % value(f_form, point, out)
  call out % get_real(rv)
  call report(matches(rv, [7.5_dp]), &
       & "the evaluator answers F value: F = 1/2 sum q^2 + xi", nfail)

  call report(out % nentries == 1 .and. out % ncomp == 1, &
       & "the evaluator admits the functional shape [1,1]", nfail)

  !-------------------------------------------------------------------!
  ! Order 0 through the driver is the value itself.
  !-------------------------------------------------------------------!

  request = gti_partial_request()

  call evaluator % value(r_form, point, out0)
  call out0 % get_real(rv0)
  call evaluator % partial_action(r_form, point, request, no_directions, out)
  call out % get_real(rv)
  call report(matches(rv, rv0), &
       & "order 0 through the evaluator is R's value", nfail)

  call evaluator % value(f_form, point, out0)
  call out0 % get_real(rv0)
  call evaluator % partial_action(f_form, point, request, no_directions, out)
  call out % get_real(rv)
  call report(matches(rv, rv0), &
       & "order 0 through the evaluator is F's value", nfail)

  !-------------------------------------------------------------------!
  ! Order 1 and order 2 in the q direction, through the driver.
  !-------------------------------------------------------------------!

  call dir_q_v1 % values % set_real([1.0_dp, 0.0_dp, 2.0_dp])
  call dir_q_v2 % values % set_real([3.0_dp, 1.0_dp, 1.0_dp])

  request = gti_partial_request(order=1, argument_kind=[GTI_ARG_STATE], &
       & state_component=[GTI_STATE_Q])
  call evaluator % partial_action(r_form, point, request, [dir_q_v1], out)
  call out % get_real(rv)
  call report(matches(rv, [2.0_dp, 0.0_dp, 12.0_dp]), &
       & "order 1 q action on R through the evaluator: 2 q v", nfail)

  call evaluator % partial_action(f_form, point, request, [dir_q_v1], out)
  call out % get_real(rv)
  call report(matches(rv, [7.0_dp]), &
       & "order 1 q action on F through the evaluator: sum q v", nfail)

  request = gti_partial_request(order=2, &
       & argument_kind=[GTI_ARG_STATE, GTI_ARG_STATE], &
       & state_component=[GTI_STATE_Q, GTI_STATE_Q])
  call evaluator % partial_action(r_form, point, request, [dir_q_v1, dir_q_v2], out)
  call out % get_real(rv)
  call report(matches(rv, [6.0_dp, 0.0_dp, 4.0_dp]), &
       & "order 2 q,q action on R through the evaluator: 2 v1 v2", nfail)

  call evaluator % partial_action(f_form, point, request, [dir_q_v1, dir_q_v2], out)
  call out % get_real(rv)
  call report(matches(rv, [5.0_dp]), &
       & "order 2 q,q action on F through the evaluator: sum v1 v2", nfail)

  write(*,'(1x,a)') "============================================="
  if (nfail .eq. 0) then
     write(*,'(1x,a)') "all gti evaluator checks passed"
  else
     error stop
  end if

contains

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

end program test_gti_form_evaluator
