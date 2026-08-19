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
!      degree 3 :  D Phi [x^(3)] + 3 D^2 Phi [x^(1), x^(2)]
!                                +   D^3 Phi [x^(1), x^(1), x^(1)]
!      degree 4 :  D Phi [x^(4)] + 4 D^2 Phi [x^(1), x^(3)]
!                                + 3 D^2 Phi [x^(2), x^(2)]
!                                + 6 D^3 Phi [x^(1), x^(1), x^(2)]
!                                +   D^4 Phi [x^(1), ..., x^(1)]
!
! Each symbolic term is one PATTERN - a coefficient and a tuple of
! slot orders - and every slot expands over all channels providing
! that derivative seat, in ordered tuples:
!
!      D^2 Phi [x^(1), x^(1)]  =  sum_a,b D^2 Phi [x_a^(1), x_b^(1)]
!
! so (q, xi) and (xi, q) are both assembled, which is exactly the
! factor of two a symmetric mixed partial carries; the pattern
! coefficient multiplies OUTSIDE the form's partial action and
! carries the multinomial count of the composition - never any
! derivative logic of the user's.
!
! The input is a bundle of channels, one per argument of Phi:
!
!      channel % derivative(k) = x^(k)   the k-th path derivative
!
! each entry a self-describing direction (argument kind, state
! component, values). A seat with no values is an absent
! derivative, read as zero: the terms it would have fed are simply
! not assembled, and a seat above the requested degree is ignored,
! never refused.
!
! Three input laws hold before any form call. Every present seat
! names its argument from the GTI_ARG vocabulary. A channel is one
! argument's path, not a bag of unrelated directions: every
! present derivative in it must name the same argument. And two
! channels naming the same argument are refused - a duplicated
! channel would double-count chain-rule terms. Vocabulary,
! consistency, and distinctness are laws of the input, not
! conventions.
!
! Every partial action goes through the form evaluator, so each
! term is held to the form's declared output shape before it is
! accumulated; accumulation itself re-checks shape against the
! running output, because the declaration is re-read per call and
! nothing else pins it across calls.
!
! The assembler carries nothing: no form, no graph, no scheme, no
! solver, no map. Local assembly only - no time graph, no Newton,
! no traversal lives here, no tensor is materialized, and this
! seat serves degrees 0 through 4.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_chain_rule_assemblies

  use iso_fortran_env      , only : dp => REAL64
  use gti_value_buffers    , only : gti_value_buffer
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle, &
       & GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_ARG_TIME, GTI_ARG_GEOM
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

  !===================================================================!
  ! One symbolic term of the composition: a multinomial coefficient
  ! and the derivative order each slot carries. Private: patterns
  ! are how the assembler thinks, not part of its contract.
  !===================================================================!

  type :: chain_term_pattern

     integer :: coefficient = 1
     integer, allocatable :: slot_degree(:)

  end type chain_term_pattern

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
  ! consistent and distinct, and only then does any form call run.
  !===================================================================!

  subroutine assemble(this, form, point, degree, input, output)

    class(gti_chain_rule_assembler), intent(in)    :: this
    class(gti_differentiable_form) , intent(in)    :: form
    type(gti_evaluation_point)     , intent(in)    :: point
    integer                        , intent(in)    :: degree
    type(gti_chain_input)          , intent(in)    :: input
    type(gti_value_buffer)         , intent(inout) :: output

    type(gti_form_evaluator) :: evaluator
    type(chain_term_pattern), allocatable :: patterns(:)
    integer :: p

    if (degree < 0 .or. degree > 4) then
       error stop 'gti_chain_rule_assembler: degree is supported'
    end if

    call require_consistent_channels(input)
    call require_distinct_channels(input)

    !----------------------------------------------------------------!
    ! Degree 0 is the value itself, held to its shape by the
    ! evaluator.
    !----------------------------------------------------------------!

    if (degree == 0) then
       call evaluator % value(form, point, output)
       return
    end if

    !----------------------------------------------------------------!
    ! Every higher degree is its pattern table: zero the output,
    ! then accumulate one expanded pattern at a time.
    !----------------------------------------------------------------!

    call zero_output(form, output)

    call get_patterns(degree, patterns)

    do p = 1, size(patterns)
       call assemble_pattern(form, point, patterns(p), input, output)
    end do

  end subroutine assemble

  !===================================================================!
  ! The pattern table of the composition, one row per symbolic
  ! term, the coefficient carrying the multinomial count. Private:
  ! how the assembler thinks, not what it promises.
  !===================================================================!

  subroutine get_patterns(degree, patterns)

    integer                              , intent(in)  :: degree
    type(chain_term_pattern), allocatable, intent(out) :: patterns(:)

    select case (degree)
    case (1)
       allocate(patterns(1))
       patterns(1) % coefficient = 1;  patterns(1) % slot_degree = [1]
    case (2)
       allocate(patterns(2))
       patterns(1) % coefficient = 1;  patterns(1) % slot_degree = [2]
       patterns(2) % coefficient = 1;  patterns(2) % slot_degree = [1, 1]
    case (3)
       allocate(patterns(3))
       patterns(1) % coefficient = 1;  patterns(1) % slot_degree = [3]
       patterns(2) % coefficient = 3;  patterns(2) % slot_degree = [1, 2]
       patterns(3) % coefficient = 1;  patterns(3) % slot_degree = [1, 1, 1]
    case default
       allocate(patterns(5))
       patterns(1) % coefficient = 1;  patterns(1) % slot_degree = [4]
       patterns(2) % coefficient = 4;  patterns(2) % slot_degree = [1, 3]
       patterns(3) % coefficient = 3;  patterns(3) % slot_degree = [2, 2]
       patterns(4) % coefficient = 6;  patterns(4) % slot_degree = [1, 1, 2]
       patterns(5) % coefficient = 1;  patterns(5) % slot_degree = [1, 1, 1, 1]
    end select

  end subroutine get_patterns

  !===================================================================!
  ! Expand one pattern over the channels: every slot ranges over
  ! every channel providing that derivative seat, in ordered
  ! tuples - (q, xi) and (xi, q) both - and each admitted tuple
  ! becomes one partial action, scaled by the pattern coefficient.
  !===================================================================!

  subroutine assemble_pattern(form, point, pattern, input, output)

    class(gti_differentiable_form), intent(in)    :: form
    type(gti_evaluation_point)    , intent(in)    :: point
    type(chain_term_pattern)      , intent(in)    :: pattern
    type(gti_chain_input)         , intent(in)    :: input
    type(gti_value_buffer)        , intent(inout) :: output

    integer :: chosen(size(pattern % slot_degree))
    integer :: k, j, nchannels
    logical :: admitted

    k         = size(pattern % slot_degree)
    nchannels = input % size()
    if (nchannels == 0) return

    chosen = 1

    do

       admitted = .true.
       do j = 1, k
          if (.not. input % channel(chosen(j)) % &
               & has_degree(pattern % slot_degree(j))) then
             admitted = .false.
             exit
          end if
       end do

       if (admitted) then
          call emit_term(form, point, pattern, input, chosen, output)
       end if

       ! the odometer: advance the last slot, carrying leftwards
       j = k
       do
          chosen(j) = chosen(j) + 1
          if (chosen(j) <= nchannels) exit
          chosen(j) = 1
          j = j - 1
          if (j == 0) return
       end do

    end do

  end subroutine assemble_pattern

  !===================================================================!
  ! One admitted tuple, one partial action: the request is written
  ! from the chosen seats themselves, so the two always tell one
  ! story, and the coefficient scales the accumulated term outside
  ! the form's calculus.
  !===================================================================!

  subroutine emit_term(form, point, pattern, input, chosen, output)

    class(gti_differentiable_form), intent(in)    :: form
    type(gti_evaluation_point)    , intent(in)    :: point
    type(chain_term_pattern)      , intent(in)    :: pattern
    type(gti_chain_input)         , intent(in)    :: input
    integer                       , intent(in)    :: chosen(:)
    type(gti_value_buffer)        , intent(inout) :: output

    type(gti_form_evaluator)   :: evaluator
    type(gti_partial_request)  :: request
    type(gti_direction_bundle) :: directions(size(chosen))
    type(gti_value_buffer)     :: term

    integer :: kinds(size(chosen)), components(size(chosen))
    integer :: j, k

    k = size(chosen)

    do j = 1, k
       directions(j) = input % channel(chosen(j)) % &
            & derivative(pattern % slot_degree(j))
       kinds(j)      = directions(j) % argument_kind
       components(j) = directions(j) % state_component
    end do

    request = gti_partial_request(order=k, argument_kind=kinds, &
         & state_component=components)

    call evaluator % partial_action(form, point, request, directions, term)

    call accumulate_scaled(output, term, real(pattern % coefficient, dp))

  end subroutine emit_term

  !===================================================================!
  ! The consistency law: a channel is one argument's path, not a
  ! bag of unrelated directions. Every present derivative in a
  ! channel must name the same argument - same kind, and for
  ! state the same component - and a channel that mixes arguments
  ! is refused before any form call.
  !===================================================================!

  pure subroutine require_consistent_channels(input)

    type(gti_chain_input), intent(in) :: input

    logical :: found
    integer :: c, k, kind, component

    do c = 1, input % size()

       call argument_identity(input % channel(c), found, kind, component)
       if (.not. found) cycle

       do k = 1, size(input % channel(c) % derivative)

          if (.not. input % channel(c) % has_degree(k)) cycle

          select case (input % channel(c) % derivative(k) % argument_kind)
          case (GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_ARG_TIME, GTI_ARG_GEOM)
             continue
          case default
             error stop 'gti_chain_rule_assembler: derivative seat names its argument'
          end select

          if (input % channel(c) % derivative(k) % argument_kind /= kind) then
             error stop 'gti_chain_rule_assembler: inconsistent argument channel is refused'
          end if

          if (kind == GTI_ARG_STATE) then
             if (input % channel(c) % derivative(k) % state_component /= component) then
                error stop 'gti_chain_rule_assembler: inconsistent argument channel is refused'
             end if
          end if

       end do

    end do

  end subroutine require_consistent_channels

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
  ! derivative: which argument this channel's path perturbs. The
  ! consistency law makes every present derivative equivalent, so
  ! the lowest is not a choice but a representative.
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

    if (signature(1) < 0) then
       error stop 'gti_chain_rule_assembler: output nentries is nonnegative'
    end if

    if (signature(2) < 1) then
       error stop 'gti_chain_rule_assembler: output ncomp is at least one'
    end if

    call output % set_real(spread(0.0_dp, dim=1, &
         & ncopies=signature(1) * signature(2)), ncomp=signature(2))

  end subroutine zero_output

  !===================================================================!
  ! Add one scaled term into the running output, entrywise. Shape
  ! is re-proven at every addition - the declaration is re-read
  ! per call, and nothing else pins it across calls. The
  ! coefficient is the pattern's multinomial count, applied
  ! outside the form's calculus.
  !===================================================================!

  pure subroutine accumulate_scaled(output, term, coefficient)

    type(gti_value_buffer), intent(inout) :: output
    type(gti_value_buffer), intent(in)    :: term
    real(dp)              , intent(in)    :: coefficient

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

    output % rvals = output % rvals + coefficient * term % rvals

  end subroutine accumulate_scaled

end module gti_chain_rule_assemblies
