!=====================================================================!
! GTI DIFFERENTIABLE FORM (PHASE 1) - THE PUBLIC CONTRACT
!
! A local differentiable form is
!
!      Phi(U, xi, t),      U = (q^(0), ..., q^(d)),
!
! and residual and functional forms differ by use, not by ontology:
!
!      R : state/design/time -> residual output
!      F : state/design/time -> functional output
!
! Both answer the same six questions, and this is the whole public
! surface a user implements:
!
!      name               which form
!      input_signature    the argument kinds Phi reads, in order,
!                         in the GTI_ARG_* vocabulary
!      output_signature   the output shape [nentries, ncomp] the
!                         value call fills
!      max_degree         the highest partial order served
!      value              Phi at a point
!      partial_action     D^k Phi [v_1, ..., v_k] at a point
!
! The partial action is a multilinear contraction, never a stored
! tensor:
!
!      order = 0 :  Phi(U, xi, t)
!      order = 1 :  D Phi [v1]
!      order = 2 :  D^2 Phi [v1, v2]
!      order = k :  D^k Phi [v1, ..., vk]
!
! Each direction names the argument it perturbs - a state component
! q, qdot, qddot, the design xi, time, or geometry - so mixed
! partials arrive through the same call without publishing tensor
! storage. The request is the ledger of the contraction:
!
!      order            how many directions are contracted
!      argument_kind(k) which argument slot k perturbs
!      state_component(k) which state order, read only where the
!                         slot is GTI_ARG_STATE
!
! require_supported is the phase-1 seat of max_degree validation:
! a request beyond the declared degree, a slot speaking outside
! the declared vocabulary, a request whose slot count disagrees
! with its order, or directions that tell a different story than
! their request, all die loudly before any arithmetic.
!
! What the layer above will do with these calls - time graphs,
! Newton wiring, traversals, higher-order chain-rule assembly - is
! not this module's concern, and no trace of it appears here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_form_interface

  use gti_value_buffers    , only : gti_value_buffer, GTI_VALUE_REAL
  use gti_state_bundles    , only : GTI_STATE_Q, GTI_STATE_QDOT, GTI_STATE_QDDOT
  use gti_evaluation_points, only : gti_evaluation_point

  implicit none

  private
  public :: gti_differentiable_form
  public :: gti_partial_request
  public :: gti_direction_bundle
  public :: GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_ARG_TIME, GTI_ARG_GEOM
  public :: GTI_STATE_Q, GTI_STATE_QDOT, GTI_STATE_QDDOT
  public :: GTI_VALUE_REAL

  !===================================================================!
  ! The argument kinds of Phi(U, xi, t | geometry). State and
  ! design are bundles, time is a real, geometry waits for a later
  ! phase; the vocabulary reserves its seat.
  !===================================================================!

  integer, parameter :: GTI_ARG_STATE  = 1
  integer, parameter :: GTI_ARG_DESIGN = 2
  integer, parameter :: GTI_ARG_TIME   = 3
  integer, parameter :: GTI_ARG_GEOM   = 4

  !===================================================================!
  ! The ledger of one partial contraction: which mixed partial,
  ! against how many directions. Both arrays run 1..order;
  ! state_component(k) is read only where argument_kind(k) is
  ! GTI_ARG_STATE, and a request with no state slot may leave it
  ! unallocated.
  !===================================================================!

  type :: gti_partial_request

     integer :: order = 0
     integer, allocatable :: argument_kind(:)
     integer, allocatable :: state_component(:)

  end type gti_partial_request

  !===================================================================!
  ! One direction: which argument it perturbs, and the values it
  ! perturbs it by. A direction is self-describing, and must tell
  ! the same story as the request slot it serves.
  !===================================================================!

  type :: gti_direction_bundle

     integer :: argument_kind   = GTI_ARG_STATE
     integer :: state_component = GTI_STATE_Q
     type(gti_value_buffer) :: values

  end type gti_direction_bundle

  !===================================================================!
  ! The abstract form. R and F both extend this and nothing else;
  ! they differ by use, not by ontology.
  !===================================================================!

  type, abstract :: gti_differentiable_form

   contains

     procedure(form_name_interface)            , deferred :: name
     procedure(form_input_signature_interface) , deferred :: input_signature
     procedure(form_output_signature_interface), deferred :: output_signature
     procedure(form_max_degree_interface)      , deferred :: max_degree
     procedure(form_value_interface)           , deferred :: value
     procedure(form_partial_action_interface)  , deferred :: partial_action

     procedure :: require_supported

  end type gti_differentiable_form

  abstract interface

     !----------------------------------------------------------------!
     ! Which form; an identity, not a description.
     !----------------------------------------------------------------!

     pure function form_name_interface(this) result(name)
       import :: gti_differentiable_form
       class(gti_differentiable_form), intent(in) :: this
       character(len=:), allocatable :: name
     end function form_name_interface

     !----------------------------------------------------------------!
     ! The ordered argument kinds Phi reads, in the GTI_ARG_*
     ! vocabulary. A form that ignores time omits GTI_ARG_TIME.
     !----------------------------------------------------------------!

     pure function form_input_signature_interface(this) result(signature)
       import :: gti_differentiable_form
       class(gti_differentiable_form), intent(in) :: this
       integer, allocatable :: signature(:)
     end function form_input_signature_interface

     !----------------------------------------------------------------!
     ! The output shape [nentries, ncomp]. The law a driver checks:
     ! after value(point, out), the buffer's shape equals this
     ! signature - residual-sized for R, scalar for F.
     !----------------------------------------------------------------!

     pure function form_output_signature_interface(this) result(signature)
       import :: gti_differentiable_form
       class(gti_differentiable_form), intent(in) :: this
       integer, allocatable :: signature(:)
     end function form_output_signature_interface

     !----------------------------------------------------------------!
     ! The highest partial order this form serves. Every request
     ! with 0 <= order <= max_degree must be answered; anything
     ! beyond is refused loudly.
     !----------------------------------------------------------------!

     pure function form_max_degree_interface(this) result(degree)
       import :: gti_differentiable_form
       class(gti_differentiable_form), intent(in) :: this
       integer :: degree
     end function form_max_degree_interface

     !----------------------------------------------------------------!
     ! Phi at a point. For R the output is residual data; for F it
     ! is scalar or small-vector functional data.
     !----------------------------------------------------------------!

     subroutine form_value_interface(this, point, output)
       import :: gti_differentiable_form, gti_evaluation_point, &
            & gti_value_buffer
       class(gti_differentiable_form), intent(in)    :: this
       type(gti_evaluation_point)    , intent(in)    :: point
       type(gti_value_buffer)        , intent(inout) :: output
     end subroutine form_value_interface

     !----------------------------------------------------------------!
     ! D^order Phi contracted against the directions, one direction
     ! per order. Order zero is Phi itself.
     !----------------------------------------------------------------!

     subroutine form_partial_action_interface(this, point, request, &
          & directions, output)
       import :: gti_differentiable_form, gti_evaluation_point, &
            & gti_partial_request, gti_direction_bundle, gti_value_buffer
       class(gti_differentiable_form), intent(in)    :: this
       type(gti_evaluation_point)    , intent(in)    :: point
       type(gti_partial_request)     , intent(in)    :: request
       type(gti_direction_bundle)    , intent(in)    :: directions(:)
       type(gti_value_buffer)        , intent(inout) :: output
     end subroutine form_partial_action_interface

  end interface

contains

  !===================================================================!
  ! The admission check every partial_action runs first. Five laws,
  ! refused loudly in this order:
  !
  !      0 <= order <= max_degree()
  !      one argument kind per order
  !      one direction per order
  !      every slot speaks the declared vocabulary - the argument
  !      kind is one of the four GTI_ARG_* values, and a state
  !      slot's component is one of the three GTI_STATE_* orders;
  !      the vocabulary is closed, an unknown word is not a request
  !      the request and its directions tell one story - kind
  !      agrees slot by slot, and where the slot is state, the
  !      state component is named and agrees too
  !
  ! Validation computes nothing, so it is pure; refusal is its
  ! only side.
  !===================================================================!

  pure subroutine require_supported(this, request, directions)

    class(gti_differentiable_form), intent(in) :: this
    type(gti_partial_request)     , intent(in) :: request
    type(gti_direction_bundle)    , intent(in) :: directions(:)

    integer :: nslots, k

    if (request % order < 0) then
       error stop 'gti_form_interface: a partial request has a nonnegative order'
    end if

    if (request % order > this % max_degree()) then
       error stop 'gti_form_interface: a partial request beyond the declared max_degree is refused'
    end if

    nslots = 0
    if (allocated(request % argument_kind)) nslots = size(request % argument_kind)
    if (nslots /= request % order) then
       error stop 'gti_form_interface: a partial request names one argument kind per order'
    end if

    if (size(directions) /= request % order) then
       error stop 'gti_form_interface: a partial request takes one direction per order'
    end if

    do k = 1, request % order

       select case (request % argument_kind(k))
       case (GTI_ARG_STATE, GTI_ARG_DESIGN, GTI_ARG_TIME, GTI_ARG_GEOM)
          continue
       case default
          error stop 'gti_form_interface: unknown argument kind'
       end select

       if (directions(k) % argument_kind /= request % argument_kind(k)) then
          error stop 'gti_form_interface: a request and its directions must tell one story'
       end if

       if (request % argument_kind(k) == GTI_ARG_STATE) then
          if (.not. allocated(request % state_component)) then
             error stop 'gti_form_interface: a state slot names its state component'
          end if
          if (size(request % state_component) /= request % order) then
             error stop 'gti_form_interface: a state slot names its state component'
          end if
          select case (request % state_component(k))
          case (GTI_STATE_Q, GTI_STATE_QDOT, GTI_STATE_QDDOT)
             continue
          case default
             error stop 'gti_form_interface: unknown state component'
          end select
          if (directions(k) % state_component /= request % state_component(k)) then
             error stop 'gti_form_interface: a request and its directions must tell one story'
          end if
       end if

    end do

  end subroutine require_supported

end module gti_form_interface
