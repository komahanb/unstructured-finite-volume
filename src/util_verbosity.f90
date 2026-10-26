!=====================================================================!
! Module that stores the global verbosity of the library's
! diagnostics. One parameter decides how much the mesh loader, the
! mesh and the assembler print:
!
!   0 : quiet (default) - errors and warnings only
!   1 : progress        - "reading / forming / evaluating" milestones
!   2 : debug           - full per-entity dumps (cells, faces, ...)
!
!   quiet       progress                    debug
!   |-----------|---------------------------|------------------->
!   errors      + milestones                + every entity
!
! Set it once at start-up, for example from the driver or a
! configuration file.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module util_verbosity

  implicit none

  private

  public :: verbosity, set_verbosity

  ! Read everywhere, written only through set_verbosity
  integer, protected, save :: verbosity = 0

contains

  !===================================================================!
  ! Sets the global verbosity level
  !===================================================================!

  impure subroutine set_verbosity(level)

    integer, intent(in) :: level

    verbosity = level

  end subroutine set_verbosity

end module util_verbosity
