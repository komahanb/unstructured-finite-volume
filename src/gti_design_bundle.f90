!=====================================================================!
! GTI DESIGN BUNDLE (PHASE 1)
!
! Design is the tuple
!
!      xi = (xi_1, ..., xi_m),
!
! one graph field per design component - a scalar parameter rides
! as a one-entry field, a design field rides whole. The bundle is
! that tuple, in the same slot carrier the state bundle uses: a
! slot holds one concrete field polymorphically and adds no law.
!
! A seat is not an entry, and the bundle answers the two questions
! separately:
!
!      size()          the shape: how many design components the
!                      bundle is shaped for
!      has_entries()   the occupancy: at least one seat holds a
!                      field
!
! so a bundle of declared seats with nobody sitting reports its
! size and still has no entries.
!
! The bundle holds design fields and nothing else: no optimizer,
! no gradient, no partition map.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_design_bundles

  use gti_state_bundles, only : gti_field_slot

  implicit none

  private
  public :: gti_design_bundle

  !===================================================================!
  ! The design tuple xi. The type keeps the public singular name;
  ! Fortran denies a type its host module's name, so the module
  ! speaks in the plural.
  !===================================================================!

  type :: gti_design_bundle

     type(gti_field_slot), allocatable :: component(:)

   contains

     procedure :: size        => design_size
     procedure :: has_entries => design_has_entries

  end type gti_design_bundle

contains

  !===================================================================!
  ! The shape: how many design components the bundle is shaped for.
  ! No seats at all is size zero - the undesigned problem is lawful.
  !===================================================================!

  pure function design_size(this) result(nslots)

    class(gti_design_bundle), intent(in) :: this
    integer :: nslots

    if (allocated(this % component)) then
       nslots = size(this % component)
    else
       nslots = 0
    end if

  end function design_size

  !===================================================================!
  ! The occupancy: a design entry exists when some seat holds a
  ! field. Declared-but-empty seats are not entries.
  !===================================================================!

  pure function design_has_entries(this) result(has)

    class(gti_design_bundle), intent(in) :: this
    logical :: has

    integer :: seat

    has = .false.
    if (.not. allocated(this % component)) return

    do seat = 1, size(this % component)
       if (allocated(this % component(seat) % value)) then
          has = .true.
          return
       end if
    end do

  end function design_has_entries

end module gti_design_bundles
