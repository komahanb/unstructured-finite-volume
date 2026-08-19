!=====================================================================!
! THE CONTROLLER TOYS: the scalar nonlinear residual
!
!      R = qdot + q^2 - xi
!
! whose BDF1 step from history q_h over h at xi = 1 has the closed
! root
!
!      q_u = ( -1 + sqrt(1 + 4 h (q_h + h)) ) / (2 h),
!
! the deterministic specimen every acceptance test prices exactly -
! plus P = q^2 + 1, strictly positive, so Newton has nothing to
! find and the controller's budget must fail lawfully. All partials
! exact; nothing here estimates or decides.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_toy_forms

  use iso_fortran_env  , only : dp => REAL64
  use gti_value_buffers, only : gti_value_buffer
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_interface, only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, &
       & GTI_ARG_STATE, GTI_ARG_DESIGN, &
       & GTI_STATE_Q, GTI_STATE_QDOT

  implicit none

  private
  public :: toy_qdot_square_form
  public :: toy_square_plus_form

  type, extends(gti_differentiable_form) :: toy_qdot_square_form
   contains
     procedure :: name             => qs_name
     procedure :: input_signature  => qs_input_signature
     procedure :: output_signature => qs_output_signature
     procedure :: max_degree       => qs_max_degree
     procedure :: value            => qs_value
     procedure :: partial_action   => qs_partial_action
  end type toy_qdot_square_form

  type, extends(gti_differentiable_form) :: toy_square_plus_form
   contains
     procedure :: name             => sp_name
     procedure :: input_signature  => sp_input_signature
     procedure :: output_signature => sp_output_signature
     procedure :: max_degree       => sp_max_degree
     procedure :: value            => sp_value
     procedure :: partial_action   => sp_partial_action
  end type toy_square_plus_form

contains

  !===================================================================!
  ! Shared readers and slot predicates.
  !===================================================================!

  subroutine read_scalars(point, q1, qdot1, xi1)
    type(gti_evaluation_point), intent(in)  :: point
    real(dp)                  , intent(out) :: q1, qdot1, xi1
    real(dp), allocatable :: v(:)
    call point % state % component(1 + GTI_STATE_Q) % value % get_real_vector(v)
    q1 = v(1)
    if (point % state % has_component(GTI_STATE_QDOT)) then
       call point % state % component(1 + GTI_STATE_QDOT) % value % get_real_vector(v)
       qdot1 = v(1)
    else
       qdot1 = 0.0_dp
    end if
    if (point % design % has_entries()) then
       call point % design % component(1) % value % get_real_vector(v)
       xi1 = v(1)
    else
       xi1 = 0.0_dp
    end if
  end subroutine read_scalars

  pure function slot_is(request, slot, kind, component) result(does)
    type(gti_partial_request), intent(in) :: request
    integer                  , intent(in) :: slot, kind, component
    logical :: does
    does = .false.
    if (request % argument_kind(slot) /= kind) return
    if (kind == GTI_ARG_STATE) then
       does = request % state_component(slot) == component
    else
       does = .true.
    end if
  end function slot_is

  !===================================================================!
  ! R = qdot + q^2 - xi.
  !===================================================================!

  pure function qs_name(this) result(name)
    class(toy_qdot_square_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy qdot square'
  end function qs_name

  pure function qs_input_signature(this) result(signature)
    class(toy_qdot_square_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE, GTI_ARG_DESIGN]
  end function qs_input_signature

  pure function qs_output_signature(this) result(signature)
    class(toy_qdot_square_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function qs_output_signature

  pure function qs_max_degree(this) result(degree)
    class(toy_qdot_square_form), intent(in) :: this
    integer :: degree
    degree = 2
  end function qs_max_degree

  subroutine qs_value(this, point, output)
    class(toy_qdot_square_form), intent(in)    :: this
    type(gti_evaluation_point) , intent(in)    :: point
    type(gti_value_buffer)     , intent(inout) :: output
    real(dp) :: q1, qdot1, xi1
    call read_scalars(point, q1, qdot1, xi1)
    call output % set_real([qdot1 + q1 * q1 - xi1])
  end subroutine qs_value

  subroutine qs_partial_action(this, point, request, directions, output)
    class(toy_qdot_square_form), intent(in)    :: this
    type(gti_evaluation_point) , intent(in)    :: point
    type(gti_partial_request)  , intent(in)    :: request
    type(gti_direction_bundle) , intent(in)    :: directions(:)
    type(gti_value_buffer)     , intent(inout) :: output
    real(dp), allocatable :: v1(:), v2(:)
    real(dp) :: q1, qdot1, xi1
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       call read_scalars(point, q1, qdot1, xi1)
       if (slot_is(request, 1, GTI_ARG_STATE, GTI_STATE_Q)) then
          call directions(1) % values % get_real(v1)
          call output % set_real([2.0_dp * q1 * v1(1)])
       else if (slot_is(request, 1, GTI_ARG_STATE, GTI_STATE_QDOT)) then
          call directions(1) % values % get_real(v1)
          call output % set_real([v1(1)])
       else if (request % argument_kind(1) == GTI_ARG_DESIGN) then
          call directions(1) % values % get_real(v1)
          call output % set_real([-v1(1)])
       else
          call output % set_real([0.0_dp])
       end if
    case (2)
       if (slot_is(request, 1, GTI_ARG_STATE, GTI_STATE_Q) .and. &
            & slot_is(request, 2, GTI_ARG_STATE, GTI_STATE_Q)) then
          call directions(1) % values % get_real(v1)
          call directions(2) % values % get_real(v2)
          call output % set_real([2.0_dp * v1(1) * v2(1)])
       else
          call output % set_real([0.0_dp])
       end if
    end select
  end subroutine qs_partial_action

  !===================================================================!
  ! P = q^2 + 1: strictly positive, so Newton has nothing to find
  ! and must fail lawfully - the non-convergence specimen.
  !===================================================================!

  pure function sp_name(this) result(name)
    class(toy_square_plus_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy square plus one'
  end function sp_name

  pure function sp_input_signature(this) result(signature)
    class(toy_square_plus_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE]
  end function sp_input_signature

  pure function sp_output_signature(this) result(signature)
    class(toy_square_plus_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function sp_output_signature

  pure function sp_max_degree(this) result(degree)
    class(toy_square_plus_form), intent(in) :: this
    integer :: degree
    degree = 2
  end function sp_max_degree

  subroutine sp_value(this, point, output)
    class(toy_square_plus_form), intent(in)    :: this
    type(gti_evaluation_point) , intent(in)    :: point
    type(gti_value_buffer)     , intent(inout) :: output
    real(dp) :: q1, qdot1, xi1
    call read_scalars(point, q1, qdot1, xi1)
    call output % set_real([q1 * q1 + 1.0_dp])
  end subroutine sp_value

  subroutine sp_partial_action(this, point, request, directions, output)
    class(toy_square_plus_form), intent(in)    :: this
    type(gti_evaluation_point) , intent(in)    :: point
    type(gti_partial_request)  , intent(in)    :: request
    type(gti_direction_bundle) , intent(in)    :: directions(:)
    type(gti_value_buffer)     , intent(inout) :: output
    real(dp), allocatable :: v1(:)
    real(dp) :: q1, qdot1, xi1
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       call read_scalars(point, q1, qdot1, xi1)
       if (slot_is(request, 1, GTI_ARG_STATE, GTI_STATE_Q)) then
          call directions(1) % values % get_real(v1)
          call output % set_real([2.0_dp * q1 * v1(1)])
       else
          call output % set_real([0.0_dp])
       end if
    case (2)
       call output % set_real([0.0_dp])
    end select
  end subroutine sp_partial_action

end module gti_toy_forms
