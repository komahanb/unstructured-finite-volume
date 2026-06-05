#include "scalar.fpp"

!=====================================================================!
! Flux operator F(q, grad q): the vector under the divergence in a
! conservation law  dq/dt + div F(q, grad q) = S. Extracted from
! interface_physics so each concrete law (diffusion, advection, ...)
! extends `flux` from its own interface file, matching the
! one-interface-per-file layout the rest of the framework uses.
!
! value -> F(3,nv); state jacobians dflux_dq (3,nv,nv) and dflux_dgradq
! (3,nv,3,nv); the trailing variable index is the block-coupling seam
! (diagonal today). flux extends the common pointwise `physics` base and
! shares its point_state carrier (both still in interface_physics).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module interface_flux

  use iso_fortran_env  , only : dp => REAL64
  use interface_physics, only : physics, point_state

  implicit none

  private
  public :: flux

  type, extends(physics), abstract :: flux
   contains
     procedure(flux_value_interface)   , deferred :: value
     procedure(flux_dq_interface)      , deferred :: dflux_dq
     procedure(flux_dgradq_interface)  , deferred :: dflux_dgradq
     procedure :: dflux_ddesign => flux_ddesign_zero          ! (3,nv) for design var k

     ! face-normal projections of the flux's own jacobian (the diffusivity
     ! and advection speed a finite-volume face with normal nf presents)
     procedure :: normal_diffusivity
     procedure :: normal_speed
  end type flux

  abstract interface

     pure function flux_value_interface(this, st) result(F)
       import :: flux, point_state
       class(flux)      , intent(in) :: this
       type(point_state), intent(in) :: st
       type(scalar)                  :: F(3, this % num_components)
     end function flux_value_interface

     pure function flux_dq_interface(this, st) result(dF)
       import :: flux, point_state
       class(flux)      , intent(in) :: this
       type(point_state), intent(in) :: st
       type(scalar)                  :: dF(3, this % num_components, this % num_components)
     end function flux_dq_interface

     pure function flux_dgradq_interface(this, st) result(dF)
       import :: flux, point_state
       class(flux)      , intent(in) :: this
       type(point_state), intent(in) :: st
       type(scalar)                  :: dF(3, this % num_components, 3, this % num_components)
     end function flux_dgradq_interface

  end interface

contains

  !===================================================================!
  ! Default design partial: zero (overridden by laws that depend on a
  ! design variable, e.g. diffusion's conductivity).
  !===================================================================!

  pure function flux_ddesign_zero(this, st, k) result(dF)
    class(flux)      , intent(in) :: this
    type(point_state), intent(in) :: st
    integer          , intent(in) :: k
    type(scalar)                  :: dF(3, this % num_components)
    dF = 0.0_dp
  end function flux_ddesign_zero

  !===================================================================!
  ! Normal diffusivity: keff = -nf^T (dF/dgradq) nf, the conductivity the
  ! flux presents across a face with outward normal nf for variable ivar.
  !===================================================================!

  pure real(dp) function normal_diffusivity(this, st, nf, ivar)

    class(flux)      , intent(in) :: this
    type(point_state), intent(in) :: st
    real(dp)         , intent(in) :: nf(3)
    integer          , intent(in) :: ivar

    type(scalar) :: dFg(3, this % num_components, 3, this % num_components)

    dFg = this % dflux_dgradq(st)
    normal_diffusivity = -real(dot_product(nf, matmul(dFg(:,ivar,:,ivar), nf)), dp)

  end function normal_diffusivity

  !===================================================================!
  ! Normal speed: vn = nf^T (dF/dq), the advection speed the flux presents
  ! across a face with outward normal nf for variable ivar.
  !===================================================================!

  pure real(dp) function normal_speed(this, st, nf, ivar)

    class(flux)      , intent(in) :: this
    type(point_state), intent(in) :: st
    real(dp)         , intent(in) :: nf(3)
    integer          , intent(in) :: ivar

    type(scalar) :: dFq(3, this % num_components, this % num_components)

    dFq = this % dflux_dq(st)
    normal_speed = real(dot_product(nf, dFq(:,ivar,ivar)), dp)

  end function normal_speed

end module interface_flux
