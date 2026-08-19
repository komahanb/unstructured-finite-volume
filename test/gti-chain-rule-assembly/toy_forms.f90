!=====================================================================!
! THE TOY TRIPLE: a residual form, a functional form, and a mixed
! scalar form over the same abstract interface, with exact partial
! actions to degree 2.
!
!      R_i(q, xi) = q_i^2 + xi          i = 1..n
!      F(q, xi)   = 1/2 sum q_i^2 + xi
!      M(q, xi)   = xi sum q_i
!
! The exact multilinear actions the triple serves:
!
!      D R [v_q]        = 2 q . v_q        (elementwise)
!      D R [w_xi]       = w_xi 1
!      D^2 R [v1, v2]   = 2 v1 . v2        both slots q, else 0
!
!      D F [v_q]        = sum q_k v_k
!      D F [w_xi]       = w_xi
!      D^2 F [v1, v2]   = sum v1_k v2_k    both slots q, else 0
!
!      D M [v_q]        = xi sum v_k
!      D M [w_xi]       = w_xi sum q_k
!      D^2 M [v_q, w_xi] = w_xi sum v_k    one q slot and one xi
!                                          slot, in either order;
!                                          (q,q) and (xi,xi) vanish
!
! M is the specimen with a live mixed second partial: a total
! second derivative through it must assemble BOTH ordered cross
! terms, and their sum is the factor of two.
!
! No form reads qdot, qddot, or time, so those directions contract
! to zero - a direction names the perturbed argument, and an
! argument the form ignores has a zero partial, not an error.
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
  public :: toy_mixed_form
  public :: toy_polynomial_form
  public :: toy_liar4_form
  public :: reset_toy_counters
  public :: value_calls, partial_order_1_calls, partial_order_2_calls
  public :: partial_order_3_calls, partial_order_4_calls

  !===================================================================!
  ! The call counters of the polynomial toy: how often each order
  ! of its calculus was asked. Test instrumentation only - the
  ! assembler never sees them.
  !===================================================================!

  integer :: value_calls = 0
  integer :: partial_order_1_calls = 0
  integer :: partial_order_2_calls = 0
  integer :: partial_order_3_calls = 0
  integer :: partial_order_4_calls = 0

  !===================================================================!
  ! Phi = q^4 + q^3 xi + q^2 xi^2 + q xi^3 + xi^4, scalar q and
  ! scalar xi, with EXACT partial actions to order 4: a request
  ! with a of its slots on q and b on xi answers the mixed partial
  ! M(a,b) times the product of the direction values. The degree-4
  ! chain rule has every term alive here.
  !===================================================================!

  type, extends(gti_differentiable_form) :: toy_polynomial_form
   contains
     procedure :: name             => poly_name
     procedure :: input_signature  => poly_input_signature
     procedure :: output_signature => poly_output_signature
     procedure :: max_degree       => poly_max_degree
     procedure :: value            => poly_value
     procedure :: partial_action   => poly_partial_action
  end type toy_polynomial_form

  !===================================================================!
  ! The order-4 liar: lawful to order 3, two entries at order 4 -
  ! the evaluator's output-shape law must catch it there.
  !===================================================================!

  type, extends(gti_differentiable_form) :: toy_liar4_form
   contains
     procedure :: name             => liar4_name
     procedure :: input_signature  => liar4_input_signature
     procedure :: output_signature => liar4_output_signature
     procedure :: max_degree       => liar4_max_degree
     procedure :: value            => liar4_value
     procedure :: partial_action   => liar4_partial_action
  end type toy_liar4_form

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

  !===================================================================!
  ! M : R^n x R -> R. Scalar codomain, and the one toy whose mixed
  ! second partial is alive.
  !===================================================================!

  type, extends(gti_differentiable_form) :: toy_mixed_form

   contains

     procedure :: name             => mixed_name
     procedure :: input_signature  => mixed_input_signature
     procedure :: output_signature => mixed_output_signature
     procedure :: max_degree       => mixed_max_degree
     procedure :: value            => mixed_value
     procedure :: partial_action   => mixed_partial_action

  end type toy_mixed_form

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
  ! Does slot k of the request perturb the design xi?
  !===================================================================!

  pure function perturbs_xi(request, slot) result(does)

    type(gti_partial_request), intent(in) :: request
    integer                  , intent(in) :: slot
    logical :: does

    does = request % argument_kind(slot) == GTI_ARG_DESIGN

  end function perturbs_xi

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

  !===================================================================!
  ! The mixed form M.
  !===================================================================!

  pure function mixed_name(this) result(name)

    class(toy_mixed_form), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'toy mixed'

  end function mixed_name

  pure function mixed_input_signature(this) result(signature)

    class(toy_mixed_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = [GTI_ARG_STATE, GTI_ARG_DESIGN]

  end function mixed_input_signature

  pure function mixed_output_signature(this) result(signature)

    class(toy_mixed_form), intent(in) :: this
    integer, allocatable :: signature(:)

    signature = [1, 1]

  end function mixed_output_signature

  pure function mixed_max_degree(this) result(degree)

    class(toy_mixed_form), intent(in) :: this
    integer :: degree

    degree = 2

  end function mixed_max_degree

  subroutine mixed_value(this, point, output)

    class(toy_mixed_form)     , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output

    real(dp), allocatable :: qv(:)
    real(dp) :: xi

    call read_point(point, qv, xi)

    call output % set_real([xi * sum(qv)])

  end subroutine mixed_value

  subroutine mixed_partial_action(this, point, request, directions, output)

    class(toy_mixed_form)     , intent(in)    :: this
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
          call output % set_real([xi * sum(v1)])
       else if (perturbs_xi(request, 1)) then
          call direction_values(directions, 1, v1)
          call output % set_real([v1(1) * sum(qv)])
       else
          call output % set_real([0.0_dp])
       end if

    case (2)

       !------------------------------------------------------------!
       ! The live mixed partial: one q slot and one xi slot, in
       ! either order, answer w sum v; (q,q) and (xi,xi) vanish.
       !------------------------------------------------------------!

       if (perturbs_q(request, 1) .and. perturbs_xi(request, 2)) then
          call direction_values(directions, 1, v1)
          call direction_values(directions, 2, v2)
          call output % set_real([v2(1) * sum(v1)])
       else if (perturbs_xi(request, 1) .and. perturbs_q(request, 2)) then
          call direction_values(directions, 1, v1)
          call direction_values(directions, 2, v2)
          call output % set_real([v1(1) * sum(v2)])
       else
          call output % set_real([0.0_dp])
       end if

    end select

  end subroutine mixed_partial_action

  !===================================================================!
  ! The polynomial toy.
  !===================================================================!

  subroutine reset_toy_counters()
    value_calls = 0
    partial_order_1_calls = 0
    partial_order_2_calls = 0
    partial_order_3_calls = 0
    partial_order_4_calls = 0
  end subroutine reset_toy_counters

  pure function poly_name(this) result(name)
    class(toy_polynomial_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy polynomial'
  end function poly_name

  pure function poly_input_signature(this) result(signature)
    class(toy_polynomial_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE, GTI_ARG_DESIGN]
  end function poly_input_signature

  pure function poly_output_signature(this) result(signature)
    class(toy_polynomial_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function poly_output_signature

  pure function poly_max_degree(this) result(degree)
    class(toy_polynomial_form), intent(in) :: this
    integer :: degree
    degree = 4
  end function poly_max_degree

  subroutine poly_read(point, q1, x1)
    type(gti_evaluation_point), intent(in)  :: point
    real(dp)                  , intent(out) :: q1, x1
    real(dp), allocatable :: v(:)
    if (.not. point % state % has_component(GTI_STATE_Q)) then
       error stop 'gti_toy_forms: the polynomial needs the state q'
    end if
    call point % state % component(1 + GTI_STATE_Q) % value % get_real_vector(v)
    q1 = v(1)
    if (.not. point % design % has_entries()) then
       error stop 'gti_toy_forms: the polynomial needs a design entry'
    end if
    call point % design % component(1) % value % get_real_vector(v)
    x1 = v(1)
  end subroutine poly_read

  !-------------------------------------------------------------------!
  ! The mixed partial M(a, b) of Phi at (q, xi): a derivatives in
  ! q, b in xi.
  !-------------------------------------------------------------------!

  pure function poly_mixed(a, b, q, x) result(m)
    integer , intent(in) :: a, b
    real(dp), intent(in) :: q, x
    real(dp) :: m
    select case (10 * a + b)
    case (10);  m = 4.0_dp*q**3 + 3.0_dp*q*q*x + 2.0_dp*q*x*x + x**3
    case (1);   m = q**3 + 2.0_dp*q*q*x + 3.0_dp*q*x*x + 4.0_dp*x**3
    case (20);  m = 12.0_dp*q*q + 6.0_dp*q*x + 2.0_dp*x*x
    case (11);  m = 3.0_dp*q*q + 4.0_dp*q*x + 3.0_dp*x*x
    case (2);   m = 2.0_dp*q*q + 6.0_dp*q*x + 12.0_dp*x*x
    case (30);  m = 24.0_dp*q + 6.0_dp*x
    case (21);  m = 6.0_dp*q + 4.0_dp*x
    case (12);  m = 4.0_dp*q + 6.0_dp*x
    case (3);   m = 6.0_dp*q + 24.0_dp*x
    case (40);  m = 24.0_dp
    case (31);  m = 6.0_dp
    case (22);  m = 4.0_dp
    case (13);  m = 6.0_dp
    case (4);   m = 24.0_dp
    case default; m = 0.0_dp
    end select
  end function poly_mixed

  subroutine poly_value(this, point, output)
    class(toy_polynomial_form), intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output
    real(dp) :: q1, x1
    value_calls = value_calls + 1
    call poly_read(point, q1, x1)
    call output % set_real([q1**4 + q1**3*x1 + q1*q1*x1*x1 + q1*x1**3 + x1**4])
  end subroutine poly_value

  subroutine poly_partial_action(this, point, request, directions, output)
    class(toy_polynomial_form), intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output
    real(dp), allocatable :: v(:)
    real(dp) :: q1, x1, product
    integer :: j, a, b
    logical :: dead
    call this % require_supported(request, directions)
    select case (request % order)
    case (0)
       call this % value(point, output)
       return
    case (1)
       partial_order_1_calls = partial_order_1_calls + 1
    case (2)
       partial_order_2_calls = partial_order_2_calls + 1
    case (3)
       partial_order_3_calls = partial_order_3_calls + 1
    case (4)
       partial_order_4_calls = partial_order_4_calls + 1
    end select
    call poly_read(point, q1, x1)
    a = 0
    b = 0
    dead = .false.
    product = 1.0_dp
    do j = 1, request % order
       if (request % argument_kind(j) == GTI_ARG_STATE) then
          if (request % state_component(j) == GTI_STATE_Q) then
             a = a + 1
          else
             dead = .true.
          end if
       else if (request % argument_kind(j) == GTI_ARG_DESIGN) then
          b = b + 1
       else
          dead = .true.
       end if
       call directions(j) % values % get_real(v)
       product = product * v(1)
    end do
    if (dead) then
       call output % set_real([0.0_dp])
    else
       call output % set_real([poly_mixed(a, b, q1, x1) * product])
    end if
  end subroutine poly_partial_action

  !===================================================================!
  ! The order-4 liar.
  !===================================================================!

  pure function liar4_name(this) result(name)
    class(toy_liar4_form), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'toy order-4 liar'
  end function liar4_name

  pure function liar4_input_signature(this) result(signature)
    class(toy_liar4_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [GTI_ARG_STATE]
  end function liar4_input_signature

  pure function liar4_output_signature(this) result(signature)
    class(toy_liar4_form), intent(in) :: this
    integer, allocatable :: signature(:)
    signature = [1, 1]
  end function liar4_output_signature

  pure function liar4_max_degree(this) result(degree)
    class(toy_liar4_form), intent(in) :: this
    integer :: degree
    degree = 4
  end function liar4_max_degree

  subroutine liar4_value(this, point, output)
    class(toy_liar4_form)     , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_value_buffer)    , intent(inout) :: output
    ! the liar reads nothing of its point
    associate(unread => point)
    end associate
    call output % set_real([0.0_dp])
  end subroutine liar4_value

  subroutine liar4_partial_action(this, point, request, directions, output)
    class(toy_liar4_form)     , intent(in)    :: this
    type(gti_evaluation_point), intent(in)    :: point
    type(gti_partial_request) , intent(in)    :: request
    type(gti_direction_bundle), intent(in)    :: directions(:)
    type(gti_value_buffer)    , intent(inout) :: output
    call this % require_supported(request, directions)
    if (request % order == 4) then
       ! the lie: two entries against a declared scalar
       call output % set_real([0.0_dp, 0.0_dp])
    else if (request % order == 0) then
       call this % value(point, output)
    else
       call output % set_real([0.0_dp])
    end if
  end subroutine liar4_partial_action

end module gti_toy_forms
