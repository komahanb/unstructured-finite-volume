!=====================================================================!
! GTI VALUE CARRIER (PHASE 1)
!
! A value buffer is a finite real-valued rectangle
!
!      w : {1..nentries} x {1..ncomp} -> R
!
! stored interleaved, the shape law every value store in this
! library obeys:
!
!      stored scalars = nentries * ncomp
!      position       = (entry - 1) * ncomp + component
!
! It is the phase-1 seat of two things: the output a form fills,
! and the direction values a partial action contracts against. The
! buffer carries values and shape, and nothing else - no domain, no
! graph, no ownership. A form's output is not a field until a
! caller says on which domain it lives; a later PR may adapt this
! buffer to the graph field layer, which is why the value-kind axis
! exists while only GTI_VALUE_REAL inhabits it today.
!
! The laws:
!
!      set_real   replaces values and shape together, and refuses
!                 a vector that does not fill entries times
!                 components exactly
!      get_real   answers the values, or a zero-length array when
!                 none are held - the zero length is the signal,
!                 never a conversion, never an inference
!      clear      forgets the values and returns the shape to empty
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_value_buffers

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private
  public :: gti_value_buffer
  public :: GTI_VALUE_REAL

  !===================================================================!
  ! The one value kind of phase 1. The axis outlives the singleton.
  !===================================================================!

  integer, parameter :: GTI_VALUE_REAL = 1

  !===================================================================!
  ! One buffer: a kind, a shape, and one live store. The type keeps
  ! the public singular name; Fortran denies a type its host
  ! module's name, so the module speaks in the plural.
  !===================================================================!

  type :: gti_value_buffer

     integer :: value_kind = GTI_VALUE_REAL
     integer :: nentries = 0
     integer :: ncomp = 1

     real(dp), allocatable :: rvals(:)

   contains

     procedure :: clear    => buffer_clear
     procedure :: set_real => buffer_set_real
     procedure :: get_real => buffer_get_real

  end type gti_value_buffer

contains

  !===================================================================!
  ! Forget the values and return the shape to empty. The kind
  ! remains real: there is nothing else for it to become.
  !===================================================================!

  pure subroutine buffer_clear(this)

    class(gti_value_buffer), intent(inout) :: this

    if (allocated(this % rvals)) deallocate(this % rvals)
    this % nentries   = 0
    this % ncomp      = 1
    this % value_kind = GTI_VALUE_REAL

  end subroutine buffer_clear

  !===================================================================!
  ! Replace the values and the shape together. The vector must fill
  ! entries times components exactly - the shape law is refused
  ! loudly, never repaired silently.
  !===================================================================!

  pure subroutine buffer_set_real(this, values, ncomp)

    class(gti_value_buffer), intent(inout) :: this
    real(dp)               , intent(in)    :: values(:)
    integer , optional     , intent(in)    :: ncomp

    integer :: width

    width = 1
    if (present(ncomp)) width = ncomp

    if (width < 1) then
       error stop 'gti_value_buffer: a component count is at least one'
    end if

    if (mod(size(values), width) /= 0) then
       error stop 'gti_value_buffer: a value vector must fill entries times components exactly'
    end if

    this % rvals      = values
    this % ncomp      = width
    this % nentries   = size(values) / width
    this % value_kind = GTI_VALUE_REAL

  end subroutine buffer_set_real

  !===================================================================!
  ! Answer the values, or a zero-length array when none are held.
  !===================================================================!

  pure subroutine buffer_get_real(this, values)

    class(gti_value_buffer), intent(in)  :: this
    real(dp), allocatable  , intent(out) :: values(:)

    if (this % value_kind == GTI_VALUE_REAL .and. allocated(this % rvals)) then
       values = this % rvals
    else
       allocate(values(0))
    end if

  end subroutine buffer_get_real

end module gti_value_buffers
