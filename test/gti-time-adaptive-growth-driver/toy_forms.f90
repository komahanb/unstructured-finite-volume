!=====================================================================!
! THE ADJOINT TOYS: exact differentiable residuals and scalar
! functionals over a VECTOR design xi, one design value per
! equation.
!
!      R_i = q_i + qdot_i + xi_i      D_q R[v] = v
!                                     D_qdot R[v] = v
!                                     D_xi R[v] = v
!
!      F   = sum q_i + sum xi_i       D_q F[v] = sum v
!                                     D_xi F[v] = sum v
!                                     D_qdot F = 0
!
!      S   = q(1)^2 + xi(1) - 4       the nonlinear scalar residual
!      G   = q(1)                     the scalar functional with no
!                                     design in it at all
!      C_i = xi_i                     design-only: J_q = 0, the
!                                     singular specimen
!      short / scalar-wide / wide     the refusal shapes
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
       & GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_ARG_TIME, &
       & GTI_STATE_Q, GTI_STATE_QDOT

  implicit none

  private
  public :: toy_time_residual_form
  public :: toy_sum_functional
  public :: toy_sum_time_functional
  public :: toy_square_design_form
  public :: toy_q_functional
  public :: toy_design_only_form
  public :: toy_short_form
  public :: toy_scalar_wide_form
  public :: toy_wide_form
  public :: toy_square_plus_form

  !===================================================================!
  ! P = q(1)^2 + 1: strictly positive, so Newton has nothing to
  ! find and must fail lawfully - the non-convergence specimen.
  !===================================================================!

  type, extends(gti_differentiable_form) :: toy_square_plus_form
   contains
     procedure :: name             => square_plus_name
     procedure :: input_signature  => square_plus_input_signature
     procedure :: output_signature => square_plus_output_signature
     procedure :: max_degree       => square_plus_max_degree
     procedure :: value            => square_plus_value
     procedure :: partial_action   => square_plus_partial_action
  end type toy_square_plus_form

  !===================================================================!
  ! F = sum q + sum xi + t: the sum functional with a live clock,
  ! the witness that the term's evaluation time reaches the form.
  !===================================================================!

  type, extends(gti_differentiable_form) :: toy_sum_time_functional
   contains
     procedure :: name             => sum_time_name
     procedure :: input_signature  => sum_time_input_signature
     procedure :: output_signature => sum_time_output_signature
     procedure :: max_degree       => sum_time_max_degree
     procedure :: value            => sum_time_value
     procedure :: partial_action   => sum_time_partial_action
  end type toy_sum_time_functional

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

  type, extends(gti_differentiable_form) :: toy_sum_functional
   contains
     procedure :: name             => sum_name
     procedure :: input_signature  => sum_input_signature
     procedure :: output_signature => sum_output_signature
     procedure :: max_degree       => sum_max_degree
     procedure :: value            => sum_value
     procedure :: partial_action   => sum_partial_action
  end type toy_sum_functional

  type, extends(gti_differentiable_form) :: toy_square_design_form
   contains
     procedure :: name             => square_name
     procedure :: input_signature  => square_input_signature
     procedure :: output_signature => square_output_signature
     procedure :: max_degree       => square_max_degree
     procedure :: value            => square_value
     procedure :: partial_action   => square_partial_action
  end type toy_square_design_form

  type, extends(gti_differentiable_form) :: toy_q_functional
   contains
     procedure :: name             => qf_name
     procedure :: input_signature  => qf_input_signature
     procedure :: output_signature => qf_output_signature
     procedure :: max_degree       => qf_max_degree
     procedure :: value            => qf_value
     procedure :: partial_action   => qf_partial_action
  end type toy_q_functional

  type, extends(gti_differentiable_form) :: toy_design_only_form
   contains
     procedure :: name             => design_only_name
     procedure :: input_signature  => design_only_input_signature
     procedure :: output_signature => design_only_output_signature
     procedure :: max_degree       => design_only_max_degree
     procedure :: value            => design_only_value
     procedure :: partial_action   => design_only_partial_action
  end type toy_design_only_form

  type, extends(gti_differentiable_form) :: toy_short_form
   contains
     procedure :: name             => short_name
     procedure :: input_signature  => short_input_signature
     procedure :: output_signature => short_output_signature
     procedure :: max_degree       => short_max_degree
     procedure :: value            => short_value
     procedure :: partial_action   => short_partial_action
  end type toy_short_form

  type, extends(gti_differentiable_form) :: toy_scalar_wide_form
   contains
     procedure :: name             => scalar_wide_name
     procedure :: input_signature  => scalar_wide_input_signature
     procedure :: output_signature => scalar_wide_output_signature
     procedure :: max_degree       => scalar_wide_max_degree
     procedure :: value            => scalar_wide_value
     procedure :: partial_action   => scalar_wide_partial_action
  end type toy_scalar_wide_form

  type, extends(gti_differentiable_form) :: toy_wide_form
   contains
     procedure :: name             => wide_name
     procedure :: input_signature  => wide_input_signature
     procedure :: output_signature => wide_output_signature
     procedure :: max_degree       => wide_max_degree
     procedure :: value            => wide_value
     procedure :: partial_action   => wide_partial_action
  end type toy_wide_form

contains

  !===================================================================!
  ! Shared readers.
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

  subroutine read_xiv(point, xv)
    type(gti_evaluation_point), intent(in)  :: point
    real(dp), allocatable     , intent(out) :: xv(:)
    if (.not. point % design % has_entries()) then
       error stop 'gti_toy_forms: the toys need a design entry'
    end if
    call point % design % component(1) % value % get_real_vector(xv)
    if (size(xv) < 1) then
       error stop 'gti_toy_forms: the toys need design values'
    end if
  end subroutine read_xiv

  pure function perturbs_q(request, slot) result(does)
    type(gti_partial_request), intent(in) :: request
    integer                  , intent(in) :: slot
    logical :: does
    does = .false.
    if (request % argument_kind(slot) /= GTI_ARG_STATE) return
    does = request % state_component(slot) == GTI_STATE_Q
  end function perturbs_q

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
  ! F = sum q + sum xi, and F never reads qdot.
  !===================================================================!

  pure function sum_name(this) result(name)
    class(toy_sum_functional), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy sum functional'
  end function sum_name

  pure function sum_input_signature(this) result(signature)
    class(toy_sum_functional), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE, GTI_ARG_DESIGN]
  end function sum_input_signature

  pure function sum_output_signature(this) result(signature)
    class(toy_sum_functional), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function sum_output_signature

  pure function sum_max_degree(this) result(degree)
    class(toy_sum_functional), intent(in) :: this
    integer :: degree
    degree = 2
  end function sum_max_degree

  subroutine sum_value(this, point, output)
    class(toy_sum_functional) , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output
    real(dp), allocatable :: qv(:), qdotv(:), xv(:)
    call read_q_qdot(point, qv, qdotv)
    call read_xiv(point, xv)
    call output % set_real([sum(qv) + sum(xv)])
  end subroutine sum_value

  subroutine sum_partial_action(this, point, request, directions, output)
    class(toy_sum_functional) , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output
    real(dp), allocatable :: v1(:)
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       if (perturbs_q(request, 1) .or. &
            & request % argument_kind(1) == GTI_ARG_DESIGN) then
          call directions(1) % values % get_real(v1)
          call output % set_real([sum(v1)])
       else
          call output % set_real([0.0_dp])
       end if
    case (2)
       call output % set_real([0.0_dp])
    end select
  end subroutine sum_partial_action

  !===================================================================!
  ! S = q(1)^2 + xi(1) - 4, the nonlinear scalar residual.
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
    real(dp), allocatable :: qv(:), qdotv(:), xv(:)
    call read_q_qdot(point, qv, qdotv)
    call read_xiv(point, xv)
    call output % set_real([qv(1) * qv(1) + xv(1) - 4.0_dp])
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
  ! G = q(1): the scalar functional with no design in it.
  !===================================================================!

  pure function qf_name(this) result(name)
    class(toy_q_functional), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy q functional'
  end function qf_name

  pure function qf_input_signature(this) result(signature)
    class(toy_q_functional), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE]
  end function qf_input_signature

  pure function qf_output_signature(this) result(signature)
    class(toy_q_functional), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function qf_output_signature

  pure function qf_max_degree(this) result(degree)
    class(toy_q_functional), intent(in) :: this
    integer :: degree
    degree = 2
  end function qf_max_degree

  subroutine qf_value(this, point, output)
    class(toy_q_functional)   , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output
    real(dp), allocatable :: qv(:), qdotv(:)
    call read_q_qdot(point, qv, qdotv)
    call output % set_real([qv(1)])
  end subroutine qf_value

  subroutine qf_partial_action(this, point, request, directions, output)
    class(toy_q_functional)   , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output
    real(dp), allocatable :: v1(:)
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       if (perturbs_q(request, 1)) then
          call directions(1) % values % get_real(v1)
          call output % set_real([v1(1)])
       else
          call output % set_real([0.0_dp])
       end if
    case (2)
       call output % set_real([0.0_dp])
    end select
  end subroutine qf_partial_action

  !===================================================================!
  ! C_i = xi_i: no q anywhere, J_q = 0 while R_xi lives.
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
    real(dp), allocatable :: xv(:)
    call read_xiv(point, xv)
    call output % set_real(xv(1:3))
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
          call output % set_real(v1(1:3))
       else
          call output % set_real([0.0_dp, 0.0_dp, 0.0_dp])
       end if
    case (2)
       call output % set_real([0.0_dp, 0.0_dp, 0.0_dp])
    end select
  end subroutine design_only_partial_action

  !===================================================================!
  ! The three refusal shapes: short [2,1], scalar-wide [1,2], and
  ! wide [3,2]. They read nothing of their points.
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

  !===================================================================!
  ! F = sum q + sum xi + t.
  !===================================================================!

  pure function sum_time_name(this) result(name)
    class(toy_sum_time_functional), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy sum time functional'
  end function sum_time_name

  pure function sum_time_input_signature(this) result(signature)
    class(toy_sum_time_functional), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_ARG_TIME]
  end function sum_time_input_signature

  pure function sum_time_output_signature(this) result(signature)
    class(toy_sum_time_functional), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function sum_time_output_signature

  pure function sum_time_max_degree(this) result(degree)
    class(toy_sum_time_functional), intent(in) :: this
    integer :: degree
    degree = 2
  end function sum_time_max_degree

  subroutine sum_time_value(this, point, output)
    class(toy_sum_time_functional), intent(in)    :: this
    type(gti_evaluation_point)    , intent(in)    :: point
    type(gti_value_buffer)        , intent(inout) :: output
    real(dp), allocatable :: qv(:), qdotv(:), xv(:)
    call read_q_qdot(point, qv, qdotv)
    call read_xiv(point, xv)
    call output % set_real([sum(qv) + sum(xv) + point % time])
  end subroutine sum_time_value

  subroutine sum_time_partial_action(this, point, request, directions, output)
    class(toy_sum_time_functional), intent(in)    :: this
    type(gti_evaluation_point)    , intent(in)    :: point
    type(gti_partial_request)     , intent(in)    :: request
    type(gti_direction_bundle)    , intent(in)    :: directions(:)
    type(gti_value_buffer)        , intent(inout) :: output
    real(dp), allocatable :: v1(:)
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
    case (1)
       if (perturbs_q(request, 1) .or. &
            & request % argument_kind(1) == GTI_ARG_DESIGN) then
          call directions(1) % values % get_real(v1)
          call output % set_real([sum(v1)])
       else if (request % argument_kind(1) == GTI_ARG_TIME) then
          call directions(1) % values % get_real(v1)
          call output % set_real([v1(1)])
       else
          call output % set_real([0.0_dp])
       end if
    case (2)
       call output % set_real([0.0_dp])
    end select
  end subroutine sum_time_partial_action

  pure function scalar_wide_name(this) result(name)
    class(toy_scalar_wide_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy scalar wide'
  end function scalar_wide_name

  pure function scalar_wide_input_signature(this) result(signature)
    class(toy_scalar_wide_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE]
  end function scalar_wide_input_signature

  pure function scalar_wide_output_signature(this) result(signature)
    class(toy_scalar_wide_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 2]
  end function scalar_wide_output_signature

  pure function scalar_wide_max_degree(this) result(degree)
    class(toy_scalar_wide_form), intent(in) :: this
    integer :: degree
    degree = 2
  end function scalar_wide_max_degree

  subroutine scalar_wide_value(this, point, output)
    class(toy_scalar_wide_form), intent(in)    :: this
    type(gti_evaluation_point) , intent(in)    :: point
    type(gti_value_buffer)     , intent(inout) :: output
    associate(unread => point)
    end associate
    call output % set_real([1.0_dp, 1.0_dp], ncomp=2)
  end subroutine scalar_wide_value

  subroutine scalar_wide_partial_action(this, point, request, directions, output)
    class(toy_scalar_wide_form), intent(in)    :: this
    type(gti_evaluation_point) , intent(in)    :: point
    type(gti_partial_request)  , intent(in)    :: request
    type(gti_direction_bundle) , intent(in)    :: directions(:)
    type(gti_value_buffer)     , intent(inout) :: output
    call this % require_supported(request, directions)
    if (request % order == 0) then
       call this % value(point, output)
    else
       call output % set_real([0.5_dp, 0.5_dp], ncomp=2)
    end if
  end subroutine scalar_wide_partial_action

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
  ! P = q(1)^2 + 1.
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
       if (perturbs_q(request, 1)) then
          call directions(1) % values % get_real(v1)
          call output % set_real([2.0_dp * qv(1) * v1(1)])
       else
          call output % set_real([0.0_dp])
       end if
    case (2)
       call output % set_real([0.0_dp])
    end select
  end subroutine square_plus_partial_action

end module gti_toy_forms
