!=====================================================================!
! The refusals that must die, one per invocation:
!
!      rdegree    R is asked for order 3, one past its declared
!                 degree
!      fdegree    F is asked the same
!      dircount   an order-1 request arrives with no direction
!      mismatch   a state request arrives with a design direction
!      badkind    a slot speaks an argument kind outside the four
!                 GTI_ARG_* words, request and direction agreeing
!                 on the nonsense
!      badcomp    a state slot names a state component outside the
!                 three GTI_STATE_* orders, again in agreement
!      bshape     a value vector that does not fill entries times
!                 components
!
! Every case must error stop before any arithmetic; a case that
! returns is a failure of the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_interface   , only : gti_partial_request, gti_direction_bundle, &
       & GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_STATE_Q
  use gti_toy_forms        , only : toy_residual_form, toy_functional_form

  implicit none

  type(toy_residual_form)    :: r_form
  type(toy_functional_form)  :: f_form
  type(gti_evaluation_point) :: point
  type(gti_partial_request)  :: request
  type(gti_direction_bundle) :: three(3), one(1)
  type(gti_direction_bundle) :: none(0)
  type(gti_value_buffer)     :: out
  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

  case ('rdegree')

     request = gti_partial_request(order=3, &
          & argument_kind=[GTI_ARG_STATE, GTI_ARG_STATE, GTI_ARG_STATE], &
          & state_component=[GTI_STATE_Q, GTI_STATE_Q, GTI_STATE_Q])
     call r_form % partial_action(point, request, three, out)

  case ('fdegree')

     request = gti_partial_request(order=3, &
          & argument_kind=[GTI_ARG_STATE, GTI_ARG_STATE, GTI_ARG_STATE], &
          & state_component=[GTI_STATE_Q, GTI_STATE_Q, GTI_STATE_Q])
     call f_form % partial_action(point, request, three, out)

  case ('dircount')

     request = gti_partial_request(order=1, argument_kind=[GTI_ARG_STATE], &
          & state_component=[GTI_STATE_Q])
     call r_form % partial_action(point, request, none, out)

  case ('mismatch')

     request = gti_partial_request(order=1, argument_kind=[GTI_ARG_STATE], &
          & state_component=[GTI_STATE_Q])
     one(1) % argument_kind = GTI_ARG_DESIGN
     call r_form % partial_action(point, request, one, out)

  case ('badkind')

     request = gti_partial_request(order=1, argument_kind=[999])
     one(1) % argument_kind = 999
     call r_form % partial_action(point, request, one, out)

  case ('badcomp')

     request = gti_partial_request(order=1, argument_kind=[GTI_ARG_STATE], &
          & state_component=[999])
     one(1) % argument_kind   = GTI_ARG_STATE
     one(1) % state_component = 999
     call r_form % partial_action(point, request, one, out)

  case ('bshape')

     call out % set_real([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], ncomp=2)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
