!=====================================================================!
! GTI HIGHER-ORDER CHAIN-RULE ASSEMBLY (LOCAL SEAT)
!
! One form, one point, and per-argument derivative directions; the
! assembler produces the total derivative of Phi along the path
! x(s) of its arguments, through the public form contract alone:
!
!      degree 0 :  Phi
!      degree 1 :  D Phi [x^(1)]
!      degree 2 :  D Phi [x^(2)] + D^2 Phi [x^(1), x^(1)]
!
! expanded over the argument channels:
!
!      D Phi [x^(k)]           =  sum_c   D Phi        [x_c^(k)]
!      D^2 Phi [x^(1), x^(1)]  =  sum_a,b D^2 Phi [x_a^(1), x_b^(1)]
!
! with (a, b) running over ORDERED pairs: (q, xi) and (xi, q) are
! both assembled, which is exactly the factor of two a symmetric
! mixed partial carries in the total second derivative.
!
! The input is a bundle of channels, one per argument of Phi:
!
!      channel % derivative(1) = x^(1)   the first path derivative
!      channel % derivative(2) = x^(2)   the second path derivative
!
! each entry a self-describing direction (argument kind, state
! component, values). A seat with no values is an absent
! derivative, read as zero: the transport term it would have fed
! is simply not assembled. Two channels naming the same argument
! are refused - a duplicated channel would double-count chain-rule
! terms, so distinctness is a law of the input, not a convention.
!
! Every partial action goes through the form evaluator, so each
! term is held to the form's declared output shape before it is
! accumulated; accumulation itself re-checks shape against the
! running output, because the declaration is re-read per call and
! nothing else pins it across calls.
!
! The assembler carries nothing: no form, no graph, no scheme, no
! solver, no map. Local assembly only - no time graph, no Newton,
! no traversal lives here, and this phase serves degrees 0, 1, 2.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_chain_rule_assemblies

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, GTI_ARG_STATE
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_evaluators  , only : gti_form_evaluator

  implicit none

  private
  public :: gti_chain_rule_assembler
  public :: gti_chain_channel
  public :: gti_chain_input

  !===================================================================!
  ! One argument's path derivatives: seat k holds x^(k). A seat
  ! with no values is an absent derivative - absence is a lawful
  ! answer, and the assembly reads it as zero.
  !===================================================================!

  type :: gti_chain_channel

     type(gti_direction_bundle), allocatable :: derivative(:)

   contains

     procedure :: has_degree => channel_has_degree

  end type gti_chain_channel

  !===================================================================!
  ! The channels of one assembly, one per argument of Phi. Their
  ! argument identities must be distinct.
  !===================================================================!

  type :: gti_chain_input

     type(gti_chain_channel), allocatable :: channel(:)

   contains

     procedure :: size => input_size

  end type gti_chain_input

  !===================================================================!
  ! The stateless assembler verb. The types keep their public
  ! singular names; Fortran denies a type its host module's name,
  ! so the module speaks in the plural.
  !===================================================================!

  type :: gti_chain_rule_assembler

   contains

     procedure :: assemble

  end type gti_chain_rule_assembler

contains

  !===================================================================!
  ! Does the channel carry x^(order)? True when the seat exists
  ! and holds values.
  !===================================================================!

  pure function channel_has_degree(this, order) result(has)

    class(gti_chain_channel), intent(in) :: this
    integer                 , intent(in) :: order
    logical :: has

    has = .false.
    if (.not. allocated(this % derivative)) return
    if (order < 1 .or. order > size(this % derivative)) return

    has = this % derivative(order) % values % nentries > 0

  end function channel_has_degree

  !===================================================================!
  ! How many argument channels the input carries.
  !===================================================================!

  pure function input_size(this) result(nchannels)

    class(gti_chain_input), intent(in) :: this
    integer :: nchannels

    if (allocated(this % channel)) then
       nchannels = size(this % channel)
    else
       nchannels = 0
    end if

  end function input_size

  !===================================================================!
  ! The one public verb: assemble the total derivative of the
  ! given degree. The degree is admitted, the channels are proven
  ! distinct, and only then does any form call run.
  !===================================================================!

  subroutine assemble(this, form, point, degree, input, output)

    class(gti_chain_rule_assembler), intent(in)    :: this
    class(gti_differentiable_form) , intent(in)    :: form
    type(gti_evaluation_point)     , intent(in)    :: point
    integer                        , intent(in)    :: degree
    type(gti_chain_input)          , intent(in)    :: input
    type(gti_value_buffer)         , intent(inout) :: output

    type(gti_form_evaluator) :: evaluator
    type(gti_value_buffer)   :: term
    integer :: a, b

    if (degree < 0) then
       error stop 'gti_chain_rule_assembler: derivative degree is nonnegative'
    end if

    if (degree > 2) then
       error stop 'gti_chain_rule_assembler: degree above phase-2 support is refused'
    end if

    call require_distinct_channels(input)

    !----------------------------------------------------------------!
    ! Degree 0 is the value itself, held to its shape by the
    ! evaluator.
    !----------------------------------------------------------------!

    if (degree == 0) then
       call evaluator % value(form, point, output)
       return
    end if

    call zero_output(form, output)

    !----------------------------------------------------------------!
    ! Transport terms: D Phi [x^(degree)] - the first path
    ! derivative feeds degree 1, the second feeds degree 2. An
    ! absent derivative contributes nothing.
    !----------------------------------------------------------------!

    do a = 1, input % size()
       if (.not. input % channel(a) % has_degree(degree)) cycle
       call first_order_term(form, point, &
            & input % channel(a) % derivative(degree), term)
       call accumulate(output, term)
    end do

    !----------------------------------------------------------------!
    ! Curvature terms at degree 2: D^2 Phi [x_a^(1), x_b^(1)] over
    ! ordered pairs. (a, b) and (b, a) are both assembled - that
    ! is the factor of two on symmetric mixed partials, by
    ! construction rather than by a coefficient.
    !----------------------------------------------------------------!

    if (degree == 2) then
       do a = 1, input % size()
          if (.not. input % channel(a) % has_degree(1)) cycle
          do b = 1, input % size()
             if (.not. input % channel(b) % has_degree(1)) cycle
             call second_order_term(form, point, &
                  & input % channel(a) % derivative(1), &
                  & input % channel(b) % derivative(1), term)
             call accumulate(output, term)
          end do
       end do
    end if

  end subroutine assemble

  !===================================================================!
  ! One first-order partial action, D Phi [d], through the
  ! evaluator: the request is written from the direction itself,
  ! so the two always tell one story.
  !===================================================================!

  subroutine first_order_term(form, point, direction, term)

    class(gti_differentiable_form), intent(in)    :: form
    type(gti_evaluation_point)    , intent(in)    :: point
    type(gti_direction_bundle)    , intent(in)    :: direction
    type(gti_value_buffer)        , intent(inout) :: term

    type(gti_form_evaluator)  :: evaluator
    type(gti_partial_request) :: request

    request = gti_partial_request(order=1, &
         & argument_kind  =[direction % argument_kind], &
         & state_component=[direction % state_component])

    call evaluator % partial_action(form, point, request, [direction], term)

  end subroutine first_order_term

  !===================================================================!
  ! One second-order partial action, D^2 Phi [da, db], through the
  ! evaluator.
  !===================================================================!

  subroutine second_order_term(form, point, da, db, term)

    class(gti_differentiable_form), intent(in)    :: form
    type(gti_evaluation_point)    , intent(in)    :: point
    type(gti_direction_bundle)    , intent(in)    :: da, db
    type(gti_value_buffer)        , intent(inout) :: term

    type(gti_form_evaluator)  :: evaluator
    type(gti_partial_request) :: request

    request = gti_partial_request(order=2, &
         & argument_kind  =[da % argument_kind  , db % argument_kind], &
         & state_component=[da % state_component, db % state_component])

    call evaluator % partial_action(form, point, request, [da, db], term)

  end subroutine second_order_term

  !===================================================================!
  ! The distinctness law: two channels naming the same argument -
  ! same kind, and for state the same component - would double-
  ! count chain-rule terms, and are refused. A channel with no
  ! present derivative names nothing and collides with nothing.
  !===================================================================!

  pure subroutine require_distinct_channels(input)

    type(gti_chain_input), intent(in) :: input

    integer :: i, j

    do i = 1, input % size()
       do j = i + 1, input % size()
          if (same_argument(input % channel(i), input % channel(j))) then
             error stop 'gti_chain_rule_assembler: duplicate argument channel is refused'
          end if
       end do
    end do

  end subroutine require_distinct_channels

  pure function same_argument(a, b) result(same)

    type(gti_chain_channel), intent(in) :: a, b
    logical :: same

    logical :: found_a, found_b
    integer :: kind_a, kind_b, comp_a, comp_b

    call argument_identity(a, found_a, kind_a, comp_a)
    call argument_identity(b, found_b, kind_b, comp_b)

    same = .false.
    if (.not. (found_a .and. found_b)) return
    if (kind_a /= kind_b) return

    if (kind_a == GTI_ARG_STATE) then
       same = comp_a == comp_b
    else
       same = .true.
    end if

  end function same_argument

  !===================================================================!
  ! A channel's argument identity, read from its lowest present
  ! derivative: which argument this channel's path perturbs.
  !===================================================================!

  pure subroutine argument_identity(c, found, kind, component)

    type(gti_chain_channel), intent(in)  :: c
    logical                , intent(out) :: found
    integer                , intent(out) :: kind, component

    integer :: k

    found     = .false.
    kind      = 0
    component = 0

    if (.not. allocated(c % derivative)) return

    do k = 1, size(c % derivative)
       if (c % has_degree(k)) then
          found     = .true.
          kind      = c % derivative(k) % argument_kind
          component = c % derivative(k) % state_component
          return
       end if
    end do

  end subroutine argument_identity

  !===================================================================!
  ! Zero of the form's codomain: the accumulation starts from the
  ! declared output shape.
  !===================================================================!

  pure subroutine zero_output(form, output)

    class(gti_differentiable_form), intent(in)    :: form
    type(gti_value_buffer)        , intent(inout) :: output

    integer, allocatable :: signature(:)

    signature = form % output_signature()

    if (size(signature) /= 2) then
       error stop 'gti_chain_rule_assembler: output_signature has shape [nentries,ncomp]'
    end if

    call output % set_real(spread(0.0_dp, dim=1, &
         & ncopies=signature(1) * signature(2)), ncomp=signature(2))

  end subroutine zero_output

  !===================================================================!
  ! Add one term into the running output, entrywise. Shape is
  ! re-proven at every addition - the declaration is re-read per
  ! call, and nothing else pins it across calls.
  !===================================================================!

  pure subroutine accumulate(output, term)

    type(gti_value_buffer), intent(inout) :: output
    type(gti_value_buffer), intent(in)    :: term

    integer :: nout, nterm

    nout = 0
    if (allocated(output % rvals)) nout = size(output % rvals)

    nterm = 0
    if (allocated(term % rvals)) nterm = size(term % rvals)

    if (output % nentries /= term % nentries .or. &
         & output % ncomp /= term % ncomp .or. nout /= nterm) then
       error stop 'gti_chain_rule_assembler: output accumulation shape mismatch'
    end if

    if (nout == 0) return

    output % rvals = output % rvals + term % rvals

  end subroutine accumulate

end module gti_chain_rule_assemblies
