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
  public :: toy_square_form
  public :: toy_constant_form
  public :: toy_wide_form
  public :: toy_short_form

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

  !===================================================================!
  ! The nonlinear scalar specimen: R = [ q(1)^2 - 4 ], reading q
  ! alone, with exact partials - the one toy whose Newton needs
  ! more than a single step.
  !===================================================================!

  type, extends(gti_differentiable_form) :: toy_square_form

   contains

     procedure :: name             => square_name
     procedure :: input_signature  => square_input_signature
     procedure :: output_signature => square_output_signature
     procedure :: max_degree       => square_max_degree
     procedure :: value            => square_value
     procedure :: partial_action   => square_partial_action

  end type toy_square_form

  !===================================================================!
  ! Three refusal specimens. The constant form has a flat residual
  ! and hence a singular Jacobian; the wide form answers a lawful
  ! two-component output no Newton driver can accept as a vector;
  ! the short form answers fewer equations than unknowns.
  !===================================================================!

  type, extends(gti_differentiable_form) :: toy_constant_form

   contains

     procedure :: name             => constant_name
     procedure :: input_signature  => constant_input_signature
     procedure :: output_signature => constant_output_signature
     procedure :: max_degree       => constant_max_degree
     procedure :: value            => constant_value
     procedure :: partial_action   => constant_partial_action

  end type toy_constant_form

  type, extends(gti_differentiable_form) :: toy_wide_form

   contains

     procedure :: name             => wide_name
     procedure :: input_signature  => wide_input_signature
     procedure :: output_signature => wide_output_signature
     procedure :: max_degree       => wide_max_degree
     procedure :: value            => wide_value
     procedure :: partial_action   => wide_partial_action

  end type toy_wide_form

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

  !===================================================================!
  ! The square form: R = [ q(1)^2 - 4 ].
  !===================================================================!

  pure function square_name(this) result(name)
    class(toy_square_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy square'
  end function square_name

  pure function square_input_signature(this) result(signature)
    class(toy_square_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE]
  end function square_input_signature

  pure function square_output_signature(this) result(signature)
    class(toy_square_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function square_output_signature

  pure function square_max_degree(this) result(degree)
    class(toy_square_form), intent(in) :: this
    integer :: degree
    degree = 2
  end function square_max_degree

  subroutine square_read_q(point, q1)
    type(gti_evaluation_point), intent(in)  :: point
    real(dp)                  , intent(out) :: q1
    real(dp), allocatable :: qv(:)
    if (.not. point % state % has_component(GTI_STATE_Q)) then
       error stop 'gti_toy_forms: the square form needs the state q'
    end if
    call point % state % component(1 + GTI_STATE_Q) % value % get_real_vector(qv)
    q1 = qv(1)
  end subroutine square_read_q

  subroutine square_value(this, point, output)
    class(toy_square_form)    , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output
    real(dp) :: q1
    call square_read_q(point, q1)
    call output % set_real([q1 * q1 - 4.0_dp])
  end subroutine square_value

  subroutine square_partial_action(this, point, request, directions, output)
    class(toy_square_form)    , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output
    real(dp), allocatable :: v1(:), v2(:)
    real(dp) :: q1
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       call square_read_q(point, q1)
       if (perturbs_q(request, 1)) then
          call directions(1) % values % get_real(v1)
          call output % set_real([2.0_dp * q1 * v1(1)])
       else
          call output % set_real([0.0_dp])
       end if
    case (2)
       if (perturbs_q(request, 1) .and. perturbs_q(request, 2)) then
          call directions(1) % values % get_real(v1)
          call directions(2) % values % get_real(v2)
          call output % set_real([2.0_dp * v1(1) * v2(1)])
       else
          call output % set_real([0.0_dp])
       end if
    end select
  end subroutine square_partial_action

  !===================================================================!
  ! Does slot k of the request perturb the state component q?
  !===================================================================!

  pure function perturbs_q(request, slot) result(does)
    type(gti_partial_request), intent(in) :: request
    integer                  , intent(in) :: slot
    logical :: does
    does = .false.
    if (request % argument_kind(slot) /= GTI_ARG_STATE) return
    does = request % state_component(slot) == GTI_STATE_Q
  end function perturbs_q

  !===================================================================!
  ! The constant form: R = [1, 1, 1] whatever the point says. Its
  ! Jacobian is exactly zero - the singular-pivot specimen.
  !===================================================================!

  pure function constant_name(this) result(name)
    class(toy_constant_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy constant'
  end function constant_name

  pure function constant_input_signature(this) result(signature)
    class(toy_constant_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE]
  end function constant_input_signature

  pure function constant_output_signature(this) result(signature)
    class(toy_constant_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [3, 1]
  end function constant_output_signature

  pure function constant_max_degree(this) result(degree)
    class(toy_constant_form), intent(in) :: this
    integer :: degree
    degree = 2
  end function constant_max_degree

  subroutine constant_value(this, point, output)
    class(toy_constant_form)  , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output
    ! the constant reads nothing of its point
    associate(unread => point)
    end associate
    call output % set_real([1.0_dp, 1.0_dp, 1.0_dp])
  end subroutine constant_value

  subroutine constant_partial_action(this, point, request, directions, output)
    class(toy_constant_form)  , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output
    call this % require_supported(request, directions)
    if (request % order == 0) then
       call this % value(point, output)
    else
       call output % set_real([0.0_dp, 0.0_dp, 0.0_dp])
    end if
  end subroutine constant_partial_action

  !===================================================================!
  ! The wide form: a lawful [3, 2] output that is not a vector.
  !===================================================================!

  pure function wide_name(this) result(name)
    class(toy_wide_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy wide'
  end function wide_name

  pure function wide_input_signature(this) result(signature)
    class(toy_wide_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE]
  end function wide_input_signature

  pure function wide_output_signature(this) result(signature)
    class(toy_wide_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [3, 2]
  end function wide_output_signature

  pure function wide_max_degree(this) result(degree)
    class(toy_wide_form), intent(in) :: this
    integer :: degree
    degree = 2
  end function wide_max_degree

  subroutine wide_value(this, point, output)
    class(toy_wide_form)      , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output
    ! the wide form reads nothing of its point
    associate(unread => point)
    end associate
    call output % set_real(spread(1.0_dp, dim=1, ncopies=6), ncomp=2)
  end subroutine wide_value

  subroutine wide_partial_action(this, point, request, directions, output)
    class(toy_wide_form)      , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output
    call this % require_supported(request, directions)
    if (request % order == 0) then
       call this % value(point, output)
    else
       call output % set_real(spread(0.0_dp, dim=1, ncopies=6), ncomp=2)
    end if
  end subroutine wide_partial_action

  !===================================================================!
  ! The short form: two equations, however many unknowns arrive.
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
    degree = 2
  end function short_max_degree

  subroutine short_value(this, point, output)
    class(toy_short_form)     , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output
    ! the short form reads nothing of its point
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
       call output % set_real([0.0_dp, 0.0_dp])
    end if
  end subroutine short_partial_action

end module gti_toy_forms
