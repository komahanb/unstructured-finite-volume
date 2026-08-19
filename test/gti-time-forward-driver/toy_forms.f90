!=====================================================================!
! THE MARCH TOYS: the exact linear time residual over a vector
! design,
!
!      R_i = q_i + qdot_i + xi_i       D_q R[v] = v
!                                      D_qdot R[v] = v
!                                      D_xi R[v] = v
!
! and the rootless scalar
!
!      P = q(1)^2 + 1                  D_q P[v] = 2 q v
!
! whose Newton must fail lawfully - there is no real root to find.
! All partials exact; no finite difference exists in this suite.
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
  public :: toy_square_plus_form

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

  type, extends(gti_differentiable_form) :: toy_square_plus_form
   contains
     procedure :: name             => square_plus_name
     procedure :: input_signature  => square_plus_input_signature
     procedure :: output_signature => square_plus_output_signature
     procedure :: max_degree       => square_plus_max_degree
     procedure :: value            => square_plus_value
     procedure :: partial_action   => square_plus_partial_action
  end type toy_square_plus_form

contains

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

  subroutine read_xiv(point, xv)
    type(gti_evaluation_point), intent(in)  :: point
    real(dp), allocatable     , intent(out) :: xv(:)
    if (.not. point % design % has_entries()) then
       error stop 'gti_toy_forms: the toys need a design entry'
    end if
    call point % design % component(1) % value % get_real_vector(xv)
  end subroutine read_xiv

  !===================================================================!
  ! R_i = q_i + qdot_i + xi_i.
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
    real(dp), allocatable :: qv(:), qdotv(:), xv(:)
    call read_q_qdot(point, qv, qdotv)
    call read_xiv(point, xv)
    if (size(qv) /= this % nequations .or. size(xv) /= size(qv)) then
       error stop 'gti_toy_forms: the state and design must fill the residual domain'
    end if
    call output % set_real(qv + qdotv + xv)
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
          ! every first partial of this linear form is the identity
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
  ! P = q(1)^2 + 1: strictly positive, so Newton has nothing to
  ! find and must say so lawfully.
  !===================================================================!

  pure function square_plus_name(this) result(name)
    class(toy_square_plus_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy square plus one'
  end function square_plus_name

  pure function square_plus_input_signature(this) result(signature)
    class(toy_square_plus_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE]
  end function square_plus_input_signature

  pure function square_plus_output_signature(this) result(signature)
    class(toy_square_plus_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function square_plus_output_signature

  pure function square_plus_max_degree(this) result(degree)
    class(toy_square_plus_form), intent(in) :: this
    integer :: degree
    degree = 2
  end function square_plus_max_degree

  subroutine square_plus_value(this, point, output)
    class(toy_square_plus_form), intent(in)    :: this
    type(gti_evaluation_point) , intent(in)    :: point
    type(gti_value_buffer)     , intent(inout) :: output
    real(dp), allocatable :: qv(:), qdotv(:)
    call read_q_qdot(point, qv, qdotv)
    call output % set_real([qv(1) * qv(1) + 1.0_dp])
  end subroutine square_plus_value

  subroutine square_plus_partial_action(this, point, request, directions, output)
    class(toy_square_plus_form), intent(in)    :: this
    type(gti_evaluation_point) , intent(in)    :: point
    type(gti_partial_request)  , intent(in)    :: request
    type(gti_direction_bundle) , intent(in)    :: directions(:)
    type(gti_value_buffer)     , intent(inout) :: output
    real(dp), allocatable :: qv(:), qdotv(:), v1(:)
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       call read_q_qdot(point, qv, qdotv)
       if (request % argument_kind(1) == GTI_ARG_STATE) then
          if (request % state_component(1) == GTI_STATE_Q) then
             call directions(1) % values % get_real(v1)
             call output % set_real([2.0_dp * qv(1) * v1(1)])
          else
             call output % set_real([0.0_dp])
          end if
       else
          call output % set_real([0.0_dp])
       end if
    case (2)
       call output % set_real([0.0_dp])
    end select
  end subroutine square_plus_partial_action

end module gti_toy_forms
