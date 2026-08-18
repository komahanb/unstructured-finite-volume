!=====================================================================!
! GTI FORM EVALUATOR (PHASE 1 DRIVER SEAT)
!
! The first layer above the public form contract: a driver that
! calls R and F through class(gti_differentiable_form) and holds
! every answer to the declared output law
!
!      output_signature() = [nentries, ncomp],  nentries >= 0,
!                                               ncomp    >= 1,
!      output % nentries  = nentries,
!      output % ncomp     = ncomp,
!      stored scalars     = nentries * ncomp,
!
! checked after value and after partial_action alike - a partial
! action is output-shaped too, an order-k contraction of a form
! lands in the same codomain as the form itself.
!
! The evaluator carries nothing: no form, no graph, no scheme, no
! solver, no map. It is two verbs and one law -
!
!      value            form % value, then the shape law
!      partial_action   form % require_supported, the form's
!                       action, then the same shape law
!
! so the driver-side validation promised by the architecture
! (GTI.md 16: the driver validates before solving) exists from the
! first phase, while time graphs, Newton, traversals, and
! chain-rule assembly do not exist here at all.
!
! A form that lies - a malformed signature, a filled shape that
! contradicts the declaration, storage that contradicts the shape -
! dies loudly at this seat, never downstream in arithmetic.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_form_evaluators

  use gti_value_buffers    , only : gti_value_buffer
  use gti_evaluation_points, only : gti_evaluation_point
  use gti_form_interface   , only : gti_differentiable_form, &
       & gti_partial_request, gti_direction_bundle

  implicit none

  private
  public :: gti_form_evaluator

  !===================================================================!
  ! The stateless driver verb-pair. The type keeps the public
  ! singular name; Fortran denies a type its host module's name,
  ! so the module speaks in the plural.
  !===================================================================!

  type :: gti_form_evaluator

   contains

     procedure :: value          => evaluator_value
     procedure :: partial_action => evaluator_partial_action

  end type gti_form_evaluator

contains

  !===================================================================!
  ! Phi at a point, held to the declared output shape.
  !===================================================================!

  subroutine evaluator_value(this, form, point, output)

    class(gti_form_evaluator)     , intent(in)    :: this
    class(gti_differentiable_form), intent(in)    :: form
    type(gti_evaluation_point)    , intent(in)    :: point
    type(gti_value_buffer)        , intent(inout) :: output

    call form % value(point, output)

    call check_output_signature(form, output)

  end subroutine evaluator_value

  !===================================================================!
  ! D^order Phi contracted against the directions: the request is
  ! admitted first - the driver validates, never assumes the form
  ! did - and the answer is held to the same output shape as the
  ! value itself.
  !===================================================================!

  subroutine evaluator_partial_action(this, form, point, request, &
       & directions, output)

    class(gti_form_evaluator)     , intent(in)    :: this
    class(gti_differentiable_form), intent(in)    :: form
    type(gti_evaluation_point)    , intent(in)    :: point
    type(gti_partial_request)     , intent(in)    :: request
    type(gti_direction_bundle)    , intent(in)    :: directions(:)
    type(gti_value_buffer)        , intent(inout) :: output

    call form % require_supported(request, directions)

    call form % partial_action(point, request, directions, output)

    call check_output_signature(form, output)

  end subroutine evaluator_partial_action

  !===================================================================!
  ! The one shape law, asked after both verbs. Five refusals, in
  ! order: the declaration is two entries, its entries are lawful,
  ! the filled buffer matches the declaration, and the storage
  ! matches the shape. Validation computes nothing, so it is pure;
  ! refusal is its only side.
  !===================================================================!

  pure subroutine check_output_signature(form, output)

    class(gti_differentiable_form), intent(in) :: form
    type(gti_value_buffer)        , intent(in) :: output

    integer, allocatable :: signature(:)
    integer :: nstored

    signature = form % output_signature()

    if (size(signature) /= 2) then
       error stop 'gti_form_evaluator: output_signature has shape [nentries,ncomp]'
    end if

    if (signature(1) < 0) then
       error stop 'gti_form_evaluator: output nentries is nonnegative'
    end if

    if (signature(2) < 1) then
       error stop 'gti_form_evaluator: output ncomp is at least one'
    end if

    if (output % nentries /= signature(1) .or. output % ncomp /= signature(2)) then
       error stop 'gti_form_evaluator: form output does not match output_signature'
    end if

    nstored = 0
    if (allocated(output % rvals)) nstored = size(output % rvals)
    if (nstored /= signature(1) * signature(2)) then
       error stop 'gti_form_evaluator: output storage does not match output shape'
    end if

  end subroutine check_output_signature

end module gti_form_evaluators
