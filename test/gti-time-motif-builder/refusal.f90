!=====================================================================!
! The refusals that must die at the builder seat, one per
! invocation:
!
!      bdf0      BDF order 0, below the supported pair
!      bdf3      BDF order 3, above the supported pair
!      bdfh0     BDF with a zero step
!      bdfhneg   BDF with a negative step
!      dirkh0    DIRK with a zero step
!      dirkg0    DIRK with a zero diagonal
!      abm0      ABM order 0
!      abm3      ABM order 3
!      abmh0     ABM with a zero step
!
! Every case must error stop before any row is minted; a case
! that returns is a failure of the suite.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program refusal

  use iso_fortran_env       , only : dp => REAL64
  use gti_time_local_schemes , only : gti_time_motif
  use gti_time_motif_builders, only : gti_time_motif_builder

  implicit none

  type(gti_time_motif_builder) :: builder
  type(gti_time_motif)         :: motif
  character(len=32) :: which

  call get_command_argument(1, which)

  select case (trim(which))

  case ('bdf0')

     call builder % bdf_uniform(0, 0.5_dp, motif)

  case ('bdf3')

     call builder % bdf_uniform(3, 0.5_dp, motif)

  case ('bdfh0')

     call builder % bdf_uniform(1, 0.0_dp, motif)

  case ('bdfhneg')

     call builder % bdf_uniform(1, -0.5_dp, motif)

  case ('dirkh0')

     call builder % dirk_stage(0.5_dp, 0.0_dp, motif)

  case ('dirkg0')

     call builder % dirk_stage(0.0_dp, 2.0_dp, motif)

  case ('abm0')

     call builder % abm_corrector(0, 2.0_dp, motif)

  case ('abm3')

     call builder % abm_corrector(3, 2.0_dp, motif)

  case ('abmh0')

     call builder % abm_corrector(1, 0.0_dp, motif)

  case default

     error stop 'refusal: unknown case'

  end select

  write(*,*) 'refusal case survived: ', trim(which)

end program refusal
