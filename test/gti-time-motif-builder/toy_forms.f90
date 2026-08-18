!=====================================================================!
! THE TIME TOY: one residual form over the full state tuple,
!
!      R_i(q, qdot, qddot, xi) = q_i + qdot_i + qddot_i + xi
!
! reading each state component when present and as zero when
! absent - q alone is required. The form is linear, so its exact
! partial actions are immediate:
!
!      D R [v_q]      = v          D R [v_qdot]  = v
!      D R [v_qddot]  = v          D R [w_xi]    = w 1
!      D^2 R [., .]   = 0          every pair
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
       & GTI_STATE_Q, GTI_STATE_QDOT, GTI_STATE_QDDOT

  implicit none

  private
  public :: toy_time_residual_form

  !===================================================================!
  ! R : (q, qdot, qddot) x R -> R^n. The one configuration is the
  ! codomain dimension n, visible through output_signature.
  !===================================================================!

  type, extends(gti_differentiable_form) :: toy_time_residual_form

     integer :: nequations = 0

   contains

     procedure :: name             => residual_name
     procedure :: input_signature  => residual_input_signature
     procedure :: output_signature => residual_output_signature
     procedure :: max_degree       => residual_max_degree
     procedure :: value            => residual_value
     procedure :: partial_action   => residual_partial_action

  end type toy_time_residual_form

contains

  !===================================================================!
  ! What the toy reads of a point: q as required, qdot and qddot
  ! as zero when their seats are absent or empty, and the one
  ! design value.
  !===================================================================!

  subroutine read_time_point(point, qv, qdotv, qddotv, xi)

    type(gti_evaluation_point), intent(in)  :: point
    real(dp), allocatable     , intent(out) :: qv(:), qdotv(:), qddotv(:)
    real(dp)                  , intent(out) :: xi

    real(dp), allocatable :: xv(:)

    if (.not. point % state % has_component(GTI_STATE_Q)) then
       error stop 'gti_toy_forms: the time toy needs the state q'
    end if
    if (.not. point % design % has_entries()) then
       error stop 'gti_toy_forms: the time toy needs a design entry'
    end if

    call point % state % component(1 + GTI_STATE_Q) % value % get_real_vector(qv)

    if (point % state % has_component(GTI_STATE_QDOT)) then
       call point % state % component(1 + GTI_STATE_QDOT) % value % get_real_vector(qdotv)
    else
       qdotv = spread(0.0_dp, dim=1, ncopies=size(qv))
    end if

    if (point % state % has_component(GTI_STATE_QDDOT)) then
       call point % state % component(1 + GTI_STATE_QDDOT) % value % get_real_vector(qddotv)
    else
       qddotv = spread(0.0_dp, dim=1, ncopies=size(qv))
    end if

    call point % design % component(1) % value % get_real_vector(xv)
    if (size(xv) < 1) then
       error stop 'gti_toy_forms: the time toy needs one design value'
    end if
    xi = xv(1)

  end subroutine read_time_point

  pure function residual_name(this) result(name)

    class(toy_time_residual_form), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'toy time residual'

  end function residual_name

  pure function residual_input_signature(this) result(signature)

    class(toy_time_residual_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = [GTI_ARG_STATE, GTI_ARG_DESIGN]

  end function residual_input_signature

  pure function residual_output_signature(this) result(signature)

    class(toy_time_residual_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = [this % nequations, 1]

  end function residual_output_signature

  pure function residual_max_degree(this) result(degree)

    class(toy_time_residual_form), intent(in) :: this
    integer :: degree

    degree = 2

  end function residual_max_degree

  subroutine residual_value(this, point, output)

    class(toy_time_residual_form), intent(in)    :: this
    type(gti_evaluation_point)   , intent(in)    :: point
    type(gti_value_buffer)       , intent(inout) :: output

    real(dp), allocatable :: qv(:), qdotv(:), qddotv(:)
    real(dp) :: xi

    call read_time_point(point, qv, qdotv, qddotv, xi)
    if (size(qv) /= this % nequations) then
       error stop 'gti_toy_forms: the state must fill the residual domain exactly'
    end if

    call output % set_real(qv + qdotv + qddotv + xi)

  end subroutine residual_value

  !===================================================================!
  ! The linear form's partial actions: every first partial is the
  ! direction itself (broadcast for xi), every second partial is
  ! zero.
  !===================================================================!

  subroutine residual_partial_action(this, point, request, directions, output)

    class(toy_time_residual_form), intent(in)    :: this
    type(gti_evaluation_point)   , intent(in)    :: point
    type(gti_partial_request)    , intent(in)    :: request
    type(gti_direction_bundle)   , intent(in)    :: directions(:)
    type(gti_value_buffer)       , intent(inout) :: output

    real(dp), allocatable :: qv(:), qdotv(:), qddotv(:), v1(:)
    real(dp) :: xi

    call this % require_supported(request, directions)

    select case (request % order)

    case (0)

       call this % value(point, output)

    case (1)

       call read_time_point(point, qv, qdotv, qddotv, xi)

       if (request % argument_kind(1) == GTI_ARG_STATE) then
          call directions(1) % values % get_real(v1)
          call output % set_real(v1)
       else if (request % argument_kind(1) == GTI_ARG_DESIGN) then
          call directions(1) % values % get_real(v1)
          call output % set_real(spread(v1(1), dim=1, ncopies=size(qv)))
       else
          call output % set_real(spread(0.0_dp, dim=1, ncopies=size(qv)))
       end if

    case (2)

       call read_time_point(point, qv, qdotv, qddotv, xi)

       call output % set_real(spread(0.0_dp, dim=1, ncopies=size(qv)))

    end select

  end subroutine residual_partial_action

end module gti_toy_forms
