!=====================================================================!
! This module holds the global verbosity for library diagnostics. One
! knob decides how much the mesh loader, the mesh and the assembler
! print:
!
!   0 : quiet (default) - errors and warnings only
!   1 : progress        - "reading / forming / evaluating" milestones
!   2 : debug           - full per-entity dumps (cells, faces, ...)
!
! Set it once at start-up - for example, from the driver or a config
! file.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module module_verbosity

  implicit none

  private
  public :: verbosity, set_verbosity

  ! The knob is read everywhere and written only through
  ! set_verbosity.
  integer, protected, save :: verbosity = 0

contains

  !===================================================================!
  ! Set the global verbosity level.
  !===================================================================!

  impure subroutine set_verbosity(level)

    integer, intent(in) :: level

    verbosity = level

  end subroutine set_verbosity

end module module_verbosity
