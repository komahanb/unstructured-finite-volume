!=====================================================================!
! THE TOY PAIR: one residual form and one functional form over the
! same abstract interface, with exact partial actions to degree 2.
!
!      R_i(q, xi) = q_i^2 + xi          i = 1..n
!      F(q, xi)   = 1/2 sum q_i^2 + xi
!
! The exact multilinear actions the pair serves:
!
!      D R [v_q]        = 2 q . v_q        (elementwise)
!      D R [w_xi]       = w_xi 1
!      D^2 R [v1, v2]   = 2 v1 . v2        both slots q, else 0
!
!      D F [v_q]        = sum q_k v_k
!      D F [w_xi]       = w_xi
!      D^2 F [v1, v2]   = sum v1_k v2_k    both slots q, else 0
!
! Neither form reads qdot, qddot, or time, so those directions
! contract to zero - a direction names the perturbed argument, and
! an argument the form ignores has a zero partial, not an error.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_toy_forms

  use iso_fortran_env  , only : dp => REAL64
  use gti_value_buffers, only : gti_value_buffer
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_interface, only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, &
       & GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_STATE_Q

  implicit none

  private
  public :: toy_residual_form
  public :: toy_functional_form

  !===================================================================!
  ! R : R^n x R -> R^n. The one configuration a residual form
  ! declares is its codomain dimension n, visible through
  ! output_signature; evaluation state arrives only through the
  ! point.
  !===================================================================!

  type, extends(gti_differentiable_form) :: toy_residual_form

     integer :: nequations = 0

   contains

     procedure :: name             => residual_name
     procedure :: input_signature  => residual_input_signature
     procedure :: output_signature => residual_output_signature
     procedure :: max_degree       => residual_max_degree
     procedure :: value            => residual_value
     procedure :: partial_action   => residual_partial_action

  end type toy_residual_form

  !===================================================================!
  ! F : R^n x R -> R. Scalar codomain, no configuration at all.
  !===================================================================!

  type, extends(gti_differentiable_form) :: toy_functional_form

   contains

     procedure :: name             => functional_name
     procedure :: input_signature  => functional_input_signature
     procedure :: output_signature => functional_output_signature
     procedure :: max_degree       => functional_max_degree
     procedure :: value            => functional_value
     procedure :: partial_action   => functional_partial_action

  end type toy_functional_form

contains

  !===================================================================!
  ! What both toys read of a point: the state q and the one design
  ! value xi. Absence is refused loudly - a toy needs its
  ! arguments before it can be a function of them.
  !===================================================================!

  subroutine read_point(point, qv, xi)

    type(gti_evaluation_point), intent(in)  :: point
    real(dp), allocatable     , intent(out) :: qv(:)
    real(dp)                  , intent(out) :: xi

    real(dp), allocatable :: xv(:)

    if (.not. point % state % has_component(GTI_STATE_Q)) then
       error stop 'gti_toy_forms: the toy forms need the state q'
    end if
    if (.not. point % design % has_entries()) then
       error stop 'gti_toy_forms: the toy forms need a design entry'
    end if

    call point % state % component(1 + GTI_STATE_Q) % value % get_real_vector(qv)
    call point % design % component(1) % value % get_real_vector(xv)

    if (size(xv) < 1) then
       error stop 'gti_toy_forms: the toy forms need one design value'
    end if
    xi = xv(1)

  end subroutine read_point

  !===================================================================!
  ! One direction's values, as a plain vector.
  !===================================================================!

  subroutine direction_values(directions, slot, v)

    type(gti_direction_bundle), intent(in)  :: directions(:)
    integer                   , intent(in)  :: slot
    real(dp), allocatable     , intent(out) :: v(:)

    call directions(slot) % values % get_real(v)

  end subroutine direction_values

  !===================================================================!
  ! Does slot k of the request perturb the state component q? The
  ! state-component ledger is consulted only after the slot proves
  ! to be a state slot - a design-only request lawfully carries no
  ! ledger at all.
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
  ! The residual form R.
  !===================================================================!

  pure function residual_name(this) result(name)

    class(toy_residual_form), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'toy residual'

  end function residual_name

  pure function residual_input_signature(this) result(signature)

    class(toy_residual_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = [GTI_ARG_STATE, GTI_ARG_DESIGN]

  end function residual_input_signature

  pure function residual_output_signature(this) result(signature)

    class(toy_residual_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = [this % nequations, 1]

  end function residual_output_signature

  pure function residual_max_degree(this) result(degree)

    class(toy_residual_form), intent(in) :: this
    integer :: degree

    degree = 2

  end function residual_max_degree

  subroutine residual_value(this, point, output)

    class(toy_residual_form)  , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output

    real(dp), allocatable :: qv(:)
    real(dp) :: xi

    call read_point(point, qv, xi)
    if (size(qv) /= this % nequations) then
       error stop 'gti_toy_forms: the state must fill the residual domain exactly'
    end if

    call output % set_real(qv * qv + xi)

  end subroutine residual_value

  subroutine residual_partial_action(this, point, request, directions, output)

    class(toy_residual_form)  , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output

    real(dp), allocatable :: qv(:), v1(:), v2(:)
    real(dp) :: xi

    call this % require_supported(request, directions)

    select case (request % order)

    case (0)

       call this % value(point, output)

    case (1)

       call read_point(point, qv, xi)

       if (perturbs_q(request, 1)) then
          call direction_values(directions, 1, v1)
          call output % set_real(2.0_dp * qv * v1)
       else if (request % argument_kind(1) == GTI_ARG_DESIGN) then
          call direction_values(directions, 1, v1)
          call output % set_real(spread(v1(1), dim=1, ncopies=size(qv)))
       else
          call output % set_real(spread(0.0_dp, dim=1, ncopies=size(qv)))
       end if

    case (2)

       call read_point(point, qv, xi)

       if (perturbs_q(request, 1) .and. perturbs_q(request, 2)) then
          call direction_values(directions, 1, v1)
          call direction_values(directions, 2, v2)
          call output % set_real(2.0_dp * v1 * v2)
       else
          call output % set_real(spread(0.0_dp, dim=1, ncopies=size(qv)))
       end if

    end select

  end subroutine residual_partial_action

  !===================================================================!
  ! The functional form F.
  !===================================================================!

  pure function functional_name(this) result(name)

    class(toy_functional_form), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'toy functional'

  end function functional_name

  pure function functional_input_signature(this) result(signature)

    class(toy_functional_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = [GTI_ARG_STATE, GTI_ARG_DESIGN]

  end function functional_input_signature

  pure function functional_output_signature(this) result(signature)

    class(toy_functional_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = [1, 1]

  end function functional_output_signature

  pure function functional_max_degree(this) result(degree)

    class(toy_functional_form), intent(in) :: this
    integer :: degree

    degree = 2

  end function functional_max_degree

  subroutine functional_value(this, point, output)

    class(toy_functional_form), intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output

    real(dp), allocatable :: qv(:)
    real(dp) :: xi

    call read_point(point, qv, xi)

    call output % set_real([0.5_dp * sum(qv * qv) + xi])

  end subroutine functional_value

  subroutine functional_partial_action(this, point, request, directions, output)

    class(toy_functional_form), intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output

    real(dp), allocatable :: qv(:), v1(:), v2(:)
    real(dp) :: xi

    call this % require_supported(request, directions)

    select case (request % order)

    case (0)

       call this % value(point, output)

    case (1)

       call read_point(point, qv, xi)

       if (perturbs_q(request, 1)) then
          call direction_values(directions, 1, v1)
          call output % set_real([sum(qv * v1)])
       else if (request % argument_kind(1) == GTI_ARG_DESIGN) then
          call direction_values(directions, 1, v1)
          call output % set_real([v1(1)])
       else
          call output % set_real([0.0_dp])
       end if

    case (2)

       if (perturbs_q(request, 1) .and. perturbs_q(request, 2)) then
          call direction_values(directions, 1, v1)
          call direction_values(directions, 2, v2)
          call output % set_real([sum(v1 * v2)])
       else
          call output % set_real([0.0_dp])
       end if

    end select

  end subroutine functional_partial_action

end module gti_toy_forms
