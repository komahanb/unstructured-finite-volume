!=====================================================================!
! GTI STATE BUNDLE (PHASE 1)
!
! The state of a degree-d governing equation is the tuple
!
!      U = (q^(0), q^(1), ..., q^(d)),
!
! one graph field per differential order. The bundle is that tuple.
! Seat k holds q^(k-1), so the Fortran storage is one-based while
! the mathematical order stays zero-based:
!
!      component(1) = q          order GTI_STATE_Q     = 0
!      component(2) = qdot       order GTI_STATE_QDOT  = 1
!      component(3) = qddot      order GTI_STATE_QDDOT = 2
!
!      seat = order + 1
!
! A seat is not an occupant, and the bundle answers the two
! questions separately:
!
!      differential_degree()    the shape: seats - 1, and -1 when
!                               the bundle carries no state at all
!      has_component(order)     the occupancy: the seat exists AND
!                               holds a field
!
! Fields ride in slots because graph_field is abstract: a slot
! holds one concrete field polymorphically and adds no law of its
! own. The same slot carries design entries (gti_design_bundle) -
! it is the shared carrier, not a state concept.
!
! The bundle holds level-5 fields and nothing else: no time graph,
! no scheme, no solver, no mesh.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_state_bundles

  use graph_field_calculus, only : graph_field

  implicit none

  private
  public :: gti_field_slot
  public :: gti_state_bundle
  public :: GTI_STATE_Q, GTI_STATE_QDOT, GTI_STATE_QDDOT

  !===================================================================!
  ! The zero-based order vocabulary: q^(order).
  !===================================================================!

  integer, parameter :: GTI_STATE_Q     = 0
  integer, parameter :: GTI_STATE_QDOT  = 1
  integer, parameter :: GTI_STATE_QDDOT = 2

  !===================================================================!
  ! One slot: one concrete field, held polymorphically. An empty
  ! slot is a seat with no occupant.
  !===================================================================!

  type :: gti_field_slot

     class(graph_field), allocatable :: value

  end type gti_field_slot

  !===================================================================!
  ! The state tuple U. The type keeps the public singular name;
  ! Fortran denies a type its host module's name, so the module
  ! speaks in the plural.
  !===================================================================!

  type :: gti_state_bundle

     type(gti_field_slot), allocatable :: component(:)

   contains

     procedure :: differential_degree => state_differential_degree
     procedure :: has_component       => state_has_component

  end type gti_state_bundle

contains

  !===================================================================!
  ! The shape: a bundle of d+1 seats represents a degree-d state
  ! tuple, so degree = seats - 1. No seats at all is degree -1:
  ! not even q has a place to sit.
  !===================================================================!

  pure function state_differential_degree(this) result(degree)

    class(gti_state_bundle), intent(in) :: this
    integer :: degree

    if (allocated(this % component)) then
       degree = size(this % component) - 1
    else
       degree = -1
    end if

  end function state_differential_degree

  !===================================================================!
  ! The occupancy: q^(order) is present when its seat exists and
  ! holds a field. A seat outside the shape, or an empty seat,
  ! answers false - never an error, absence is a lawful answer.
  !===================================================================!

  pure function state_has_component(this, state_component) result(has)

    class(gti_state_bundle), intent(in) :: this
    integer                , intent(in) :: state_component
    logical :: has

    integer :: seat

    seat = state_component + 1

    has = .false.
    if (.not. allocated(this % component)) return
    if (seat < 1 .or. seat > size(this % component)) return

    has = allocated(this % component(seat) % value)

  end function state_has_component

end module gti_state_bundles
