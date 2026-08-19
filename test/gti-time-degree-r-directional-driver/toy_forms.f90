!=====================================================================!
! THE DEGREE-R TOYS: the scalar nonlinear residual
!
!      R = qdot + q^2 - xi
!
! whose exact partials are
!
!      D_q R[v] = 2 q v      D_qdot R[v] = v      D_xi R[v] = -v
!      D_qq R[v,w] = 2 v w   every other second partial zero,
!      every partial of order three through eight zero,
!
! lawful to order 8 - and although its own calculus stops at the
! quadratic, q^(3) through q^(8) do not vanish along the path: the
! chain-rule right-hand side carries products of lower-degree
! curvature, which is exactly what the degree loop must transport.
! Plus the design-only scalar (J_u = 0, the singular specimen) and
! the short form (two equations, the size-law specimen). All
! partials exact; no derivative here is approximated.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_toy_forms

  use iso_fortran_env  , only : dp => REAL64
  use gti_value_buffers, only : gti_value_buffer
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_interface, only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, &
       & GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_ARG_TIME, &
       & GTI_STATE_Q, GTI_STATE_QDOT

  implicit none

  private
  public :: toy_qdot_square_form
  public :: toy_qdot_square_time_form
  public :: toy_design_only_scalar_form
  public :: toy_short_form

  type, extends(gti_differentiable_form) :: toy_qdot_square_form
   contains
     procedure :: name             => qs_name
     procedure :: input_signature  => qs_input_signature
     procedure :: output_signature => qs_output_signature
     procedure :: max_degree       => qs_max_degree
     procedure :: value            => qs_value
     procedure :: partial_action   => qs_partial_action
  end type toy_qdot_square_form

  type, extends(gti_differentiable_form) :: toy_qdot_square_time_form
   contains
     procedure :: name             => qst_name
     procedure :: input_signature  => qst_input_signature
     procedure :: output_signature => qst_output_signature
     procedure :: max_degree       => qst_max_degree
     procedure :: value            => qst_value
     procedure :: partial_action   => qst_partial_action
  end type toy_qdot_square_time_form

  type, extends(gti_differentiable_form) :: toy_design_only_scalar_form
   contains
     procedure :: name             => do_name
     procedure :: input_signature  => do_input_signature
     procedure :: output_signature => do_output_signature
     procedure :: max_degree       => do_max_degree
     procedure :: value            => do_value
     procedure :: partial_action   => do_partial_action
  end type toy_design_only_scalar_form

  type, extends(gti_differentiable_form) :: toy_short_form
   contains
     procedure :: name             => short_name
     procedure :: input_signature  => short_input_signature
     procedure :: output_signature => short_output_signature
     procedure :: max_degree       => short_max_degree
     procedure :: value            => short_value
     procedure :: partial_action   => short_partial_action
  end type toy_short_form

contains

  !===================================================================!
  ! Shared readers and slot predicates.
  !===================================================================!

  subroutine read_scalars(point, q1, qdot1, xi1)
    type(gti_evaluation_point), intent(in)  :: point
    real(dp)                  , intent(out) :: q1, qdot1, xi1
    real(dp), allocatable :: v(:)
    if (.not. point % state % has_component(GTI_STATE_Q)) then
       error stop 'gti_toy_forms: the toys need the state q'
    end if
    call point % state % component(1 + GTI_STATE_Q) % value % get_real_vector(v)
    q1 = v(1)
    if (point % state % has_component(GTI_STATE_QDOT)) then
       call point % state % component(1 + GTI_STATE_QDOT) % value % get_real_vector(v)
       qdot1 = v(1)
    else
       qdot1 = 0.0_dp
    end if
    if (.not. point % design % has_entries()) then
       error stop 'gti_toy_forms: the toys need a design entry'
    end if
    call point % design % component(1) % value % get_real_vector(v)
    xi1 = v(1)
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
  ! R = qdot + q^2 - xi, lawful to order 4: the third and fourth
  ! partials exist and are exactly zero.
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
    degree = 8
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
    case default
       ! orders 3 through 8: a quadratic's calculus ends here, exactly
       call output % set_real([0.0_dp])
    end select
  end subroutine qs_partial_action

  !===================================================================!
  ! R = qdot + q^2 - xi - t: the time-path specimen. Identical to
  ! the qdot-square toy in every partial but one - D_t R[v] = -v,
  ! the same shape as D_xi R[v] - so a time channel and a design
  ! channel are exercised by one form without conflating them.
  !===================================================================!

  pure function qst_name(this) result(name)
    class(toy_qdot_square_time_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy qdot square time'
  end function qst_name

  pure function qst_input_signature(this) result(signature)
    class(toy_qdot_square_time_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_ARG_TIME]
  end function qst_input_signature

  pure function qst_output_signature(this) result(signature)
    class(toy_qdot_square_time_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function qst_output_signature

  pure function qst_max_degree(this) result(degree)
    class(toy_qdot_square_time_form), intent(in) :: this
    integer :: degree
    degree = 8
  end function qst_max_degree

  subroutine qst_value(this, point, output)
    class(toy_qdot_square_time_form), intent(in)    :: this
    type(gti_evaluation_point)      , intent(in)    :: point
    type(gti_value_buffer)          , intent(inout) :: output
    real(dp) :: q1, qdot1, xi1
    call read_scalars(point, q1, qdot1, xi1)
    call output % set_real([qdot1 + q1 * q1 - xi1 - point % time])
  end subroutine qst_value

  subroutine qst_partial_action(this, point, request, directions, output)
    class(toy_qdot_square_time_form), intent(in)    :: this
    type(gti_evaluation_point)      , intent(in)    :: point
    type(gti_partial_request)       , intent(in)    :: request
    type(gti_direction_bundle)      , intent(in)    :: directions(:)
    type(gti_value_buffer)          , intent(inout) :: output
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
       else if (request % argument_kind(1) == GTI_ARG_TIME) then
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
    case default
       ! orders 3 through 8: every remaining partial is exactly zero
       call output % set_real([0.0_dp])
    end select
  end subroutine qst_partial_action

  !===================================================================!
  ! R = xi: J_u = 0, the singular specimen.
  !===================================================================!

  pure function do_name(this) result(name)
    class(toy_design_only_scalar_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy design only scalar'
  end function do_name

  pure function do_input_signature(this) result(signature)
    class(toy_design_only_scalar_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_DESIGN]
  end function do_input_signature

  pure function do_output_signature(this) result(signature)
    class(toy_design_only_scalar_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function do_output_signature

  pure function do_max_degree(this) result(degree)
    class(toy_design_only_scalar_form), intent(in) :: this
    integer :: degree
    degree = 4
  end function do_max_degree

  subroutine do_value(this, point, output)
    class(toy_design_only_scalar_form), intent(in)    :: this
    type(gti_evaluation_point)        , intent(in)    :: point
    type(gti_value_buffer)            , intent(inout) :: output
    real(dp), allocatable :: v(:)
    if (.not. point % design % has_entries()) then
       error stop 'gti_toy_forms: the toys need a design entry'
    end if
    call point % design % component(1) % value % get_real_vector(v)
    call output % set_real([v(1)])
  end subroutine do_value

  subroutine do_partial_action(this, point, request, directions, output)
    class(toy_design_only_scalar_form), intent(in)    :: this
    type(gti_evaluation_point)        , intent(in)    :: point
    type(gti_partial_request)         , intent(in)    :: request
    type(gti_direction_bundle)        , intent(in)    :: directions(:)
    type(gti_value_buffer)            , intent(inout) :: output
    real(dp), allocatable :: v1(:)
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       if (request % argument_kind(1) == GTI_ARG_DESIGN) then
          call directions(1) % values % get_real(v1)
          call output % set_real([v1(1)])
       else
          call output % set_real([0.0_dp])
       end if
    case default
       call output % set_real([0.0_dp])
    end select
  end subroutine do_partial_action

  !===================================================================!
  ! The short form: two equations against one unknown.
  !===================================================================!

  pure function short_name(this) result(name)
    class(toy_short_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy short'
  end function short_name

  pure function short_input_signature(this) result(signature)
    class(toy_short_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE]
  end function short_input_signature

  pure function short_output_signature(this) result(signature)
    class(toy_short_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [2, 1]
  end function short_output_signature

  pure function short_max_degree(this) result(degree)
    class(toy_short_form), intent(in) :: this
    integer :: degree
    degree = 4
  end function short_max_degree

  subroutine short_value(this, point, output)
    class(toy_short_form)     , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output
    associate(unread => point)
    end associate
    call output % set_real([1.0_dp, 1.0_dp])
  end subroutine short_value

  subroutine short_partial_action(this, point, request, directions, output)
    class(toy_short_form)     , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output
    call this % require_supported(request, directions)
    if (request % order == 0) then
       call this % value(point, output)
    else
       call output % set_real([1.0_dp, 1.0_dp])
    end if
  end subroutine short_partial_action

end module gti_toy_forms
