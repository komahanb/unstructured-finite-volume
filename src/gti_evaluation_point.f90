!=====================================================================!
! GTI EVALUATION POINT (PHASE 1)
!
! A form sees a local point, never a trajectory:
!
!      point = ( U, xi, t | window, step, stage )
!
! The state tuple U and the design tuple xi ride as bundles, time
! is a real, and the three ids name where on a time graph the
! point stands - a graph that phase 1 does not build, which is why
! they default to zero and carry no law here.
!
! What the form must not see, the point does not own:
!
!      global sparse structure      time graph storage
!      optimizer                    adjoint vector
!      partition map                solver internals
!      mesh                         FV assembler
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module gti_evaluation_points

  use iso_fortran_env   , only : dp => REAL64
  use gti_state_bundles , only : gti_state_bundle
  use gti_design_bundles, only : gti_design_bundle

  implicit none

  private
  public :: gti_evaluation_point

  !===================================================================!
  ! The local argument of a form. The type keeps the public
  ! singular name; Fortran denies a type its host module's name,
  ! so the module speaks in the plural.
  !===================================================================!

  type :: gti_evaluation_point

     type(gti_state_bundle)  :: state
     type(gti_design_bundle) :: design
     real(dp)                :: time = 0.0_dp
     integer                 :: window_id = 0
     integer                 :: step_id = 0
     integer                 :: stage_id = 0

  end type gti_evaluation_point

end module gti_evaluation_points
