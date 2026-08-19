!=====================================================================!
! THE TANGENT TOYS: exact differentiable residuals whose design
! partial is elementwise, so a design direction carries one value
! per equation.
!
!      R_i = q_i + qdot_i + xi        D_q R[v] = v
!                                     D_qdot R[v] = v
!                                     D_xi R[v] = v
!
!      S   = q(1)^2 + xi - 4          the nonlinear scalar
!      C_i = xi                       design-only: J_q = 0, the
!                                     singular specimen
!      W   = the [3,2] wide output    not a vector
!      H   = two equations of xi      shorter than the unknown
!
! qddot is absent everywhere; all partials are exact, and no
! finite difference exists in this suite.
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
  public :: toy_time_residual_form
  public :: toy_square_design_form
  public :: toy_design_only_form
  public :: toy_wide_form
  public :: toy_short_form

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

  type, extends(gti_differentiable_form) :: toy_square_design_form
   contains
     procedure :: name             => square_name
     procedure :: input_signature  => square_input_signature
     procedure :: output_signature => square_output_signature
     procedure :: max_degree       => square_max_degree
     procedure :: value            => square_value
     procedure :: partial_action   => square_partial_action
  end type toy_square_design_form

  type, extends(gti_differentiable_form) :: toy_design_only_form
   contains
     procedure :: name             => design_only_name
     procedure :: input_signature  => design_only_input_signature
     procedure :: output_signature => design_only_output_signature
     procedure :: max_degree       => design_only_max_degree
     procedure :: value            => design_only_value
     procedure :: partial_action   => design_only_partial_action
  end type toy_design_only_form

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
  ! Shared readers: q as required, qdot as zero when absent, and
  ! the one design value.
  !===================================================================!

  subroutine read_q_qdot(point, qv, qdotv)
    type(gti_evaluation_point), intent(in)  :: point
    real(dp), allocatable     , intent(out) :: qv(:), qdotv(:)
    if (.not. point % state % has_component(GTI_STATE_Q)) then
       error stop 'gti_toy_forms: the toys need the state q'
    end if
    call point % state % component(1 + GTI_STATE_Q) % value % get_real_vector(qv)
    if (point % state % has_component(GTI_STATE_QDOT)) then
       call point % state % component(1 + GTI_STATE_QDOT) % value % get_real_vector(qdotv)
    else
       qdotv = spread(0.0_dp, dim=1, ncopies=size(qv))
    end if
  end subroutine read_q_qdot

  subroutine read_xi(point, xi)
    type(gti_evaluation_point), intent(in)  :: point
    real(dp)                  , intent(out) :: xi
    real(dp), allocatable :: xv(:)
    if (.not. point % design % has_entries()) then
       error stop 'gti_toy_forms: the toys need a design entry'
    end if
    call point % design % component(1) % value % get_real_vector(xv)
    if (size(xv) < 1) then
       error stop 'gti_toy_forms: the toys need one design value'
    end if
    xi = xv(1)
  end subroutine read_xi

  pure function perturbs_q(request, slot) result(does)
    type(gti_partial_request), intent(in) :: request
    integer                  , intent(in) :: slot
    logical :: does
    does = .false.
    if (request % argument_kind(slot) /= GTI_ARG_STATE) return
    does = request % state_component(slot) == GTI_STATE_Q
  end function perturbs_q

  !===================================================================!
  ! R_i = q_i + qdot_i + xi, design partial elementwise.
  !===================================================================!

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
    real(dp), allocatable :: qv(:), qdotv(:)
    real(dp) :: xi
    call read_q_qdot(point, qv, qdotv)
    call read_xi(point, xi)
    if (size(qv) /= this % nequations) then
       error stop 'gti_toy_forms: the state must fill the residual domain exactly'
    end if
    call output % set_real(qv + qdotv + xi)
  end subroutine residual_value

  subroutine residual_partial_action(this, point, request, directions, output)
    class(toy_time_residual_form), intent(in)    :: this
    type(gti_evaluation_point)   , intent(in)    :: point
    type(gti_partial_request)    , intent(in)    :: request
    type(gti_direction_bundle)   , intent(in)    :: directions(:)
    type(gti_value_buffer)       , intent(inout) :: output
    real(dp), allocatable :: qv(:), qdotv(:), v1(:)
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       call read_q_qdot(point, qv, qdotv)
       if (request % argument_kind(1) == GTI_ARG_STATE .or. &
            & request % argument_kind(1) == GTI_ARG_DESIGN) then
          ! every first partial of this linear form is the identity:
          ! state directions and the elementwise design direction alike
          call directions(1) % values % get_real(v1)
          call output % set_real(v1)
       else
          call output % set_real(spread(0.0_dp, dim=1, ncopies=size(qv)))
       end if
    case (2)
       call read_q_qdot(point, qv, qdotv)
       call output % set_real(spread(0.0_dp, dim=1, ncopies=size(qv)))
    end select
  end subroutine residual_partial_action

  !===================================================================!
  ! S = q(1)^2 + xi - 4.
  !===================================================================!

  pure function square_name(this) result(name)
    class(toy_square_design_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy square design'
  end function square_name

  pure function square_input_signature(this) result(signature)
    class(toy_square_design_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE, GTI_ARG_DESIGN]
  end function square_input_signature

  pure function square_output_signature(this) result(signature)
    class(toy_square_design_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function square_output_signature

  pure function square_max_degree(this) result(degree)
    class(toy_square_design_form), intent(in) :: this
    integer :: degree
    degree = 2
  end function square_max_degree

  subroutine square_value(this, point, output)
    class(toy_square_design_form), intent(in)    :: this
    type(gti_evaluation_point)   , intent(in)    :: point
    type(gti_value_buffer)       , intent(inout) :: output
    real(dp), allocatable :: qv(:), qdotv(:)
    real(dp) :: xi
    call read_q_qdot(point, qv, qdotv)
    call read_xi(point, xi)
    call output % set_real([qv(1) * qv(1) + xi - 4.0_dp])
  end subroutine square_value

  subroutine square_partial_action(this, point, request, directions, output)
    class(toy_square_design_form), intent(in)    :: this
    type(gti_evaluation_point)   , intent(in)    :: point
    type(gti_partial_request)    , intent(in)    :: request
    type(gti_direction_bundle)   , intent(in)    :: directions(:)
    type(gti_value_buffer)       , intent(inout) :: output
    real(dp), allocatable :: qv(:), qdotv(:), v1(:), v2(:)
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       call read_q_qdot(point, qv, qdotv)
       if (perturbs_q(request, 1)) then
          call directions(1) % values % get_real(v1)
          call output % set_real([2.0_dp * qv(1) * v1(1)])
       else if (request % argument_kind(1) == GTI_ARG_DESIGN) then
          call directions(1) % values % get_real(v1)
          call output % set_real([v1(1)])
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
  ! C_i = xi: a residual with no q in it at all - J_q = 0, while
  ! R_xi lives. The singular-pivot specimen.
  !===================================================================!

  pure function design_only_name(this) result(name)
    class(toy_design_only_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy design only'
  end function design_only_name

  pure function design_only_input_signature(this) result(signature)
    class(toy_design_only_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_DESIGN]
  end function design_only_input_signature

  pure function design_only_output_signature(this) result(signature)
    class(toy_design_only_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [3, 1]
  end function design_only_output_signature

  pure function design_only_max_degree(this) result(degree)
    class(toy_design_only_form), intent(in) :: this
    integer :: degree
    degree = 2
  end function design_only_max_degree

  subroutine design_only_value(this, point, output)
    class(toy_design_only_form), intent(in)    :: this
    type(gti_evaluation_point) , intent(in)    :: point
    type(gti_value_buffer)     , intent(inout) :: output
    real(dp) :: xi
    call read_xi(point, xi)
    call output % set_real([xi, xi, xi])
  end subroutine design_only_value

  subroutine design_only_partial_action(this, point, request, directions, output)
    class(toy_design_only_form), intent(in)    :: this
    type(gti_evaluation_point) , intent(in)    :: point
    type(gti_partial_request)  , intent(in)    :: request
    type(gti_direction_bundle) , intent(in)    :: directions(:)
    type(gti_value_buffer)     , intent(inout) :: output
    real(dp), allocatable :: v1(:)
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       if (request % argument_kind(1) == GTI_ARG_DESIGN) then
          call directions(1) % values % get_real(v1)
          call output % set_real(v1)
       else
          call output % set_real([0.0_dp, 0.0_dp, 0.0_dp])
       end if
    case (2)
       call output % set_real([0.0_dp, 0.0_dp, 0.0_dp])
    end select
  end subroutine design_only_partial_action

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
    signature = [GTI_ARG_DESIGN]
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
       call output % set_real(spread(0.5_dp, dim=1, ncopies=6), ncomp=2)
    end if
  end subroutine wide_partial_action

  !===================================================================!
  ! The short form: two equations of xi, whatever the unknown.
  !===================================================================!

  pure function short_name(this) result(name)
    class(toy_short_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy short'
  end function short_name

  pure function short_input_signature(this) result(signature)
    class(toy_short_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_DESIGN]
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
    real(dp) :: xi
    call read_xi(point, xi)
    call output % set_real([xi, xi])
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
