#include "scalar.fpp"

!=====================================================================!
! The flux operator F(q, grad q) is the vector under the divergence in
! a conservation law  dq/dt + div F(q, grad q) = S. It is extracted
! from interface_physics so that each concrete law (diffusion,
! advection, ...) extends `flux` from its own interface file, matching
! the one-interface-per-file layout the rest of the framework uses.
!
! The value returns F(3,nv); the state jacobians are dflux_dq (3,nv,nv)
! and dflux_dgradq (3,nv,3,nv); the trailing variable index is the
! block-coupling seam (diagonal today). The flux type extends the
! common pointwise `physics` base and shares its point_state carrier
! (both still live in interface_physics).
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
     procedure :: dflux_ddesign => flux_ddesign_zero          ! The shape is (3,nv) for design variable k.

     ! These are the face-normal projections of the flux's own jacobian:
     ! the diffusivity and the advection speed a finite-volume face with
     ! normal nf presents.
     procedure :: normal_diffusivity
     procedure :: normal_speed
  end type flux

  abstract interface

     !================================================================!
     ! Return the flux value F(3,nv) at the point state.
     !================================================================!

     pure function flux_value_interface(this, st) result(F)
       import :: flux, point_state
       class(flux)      , intent(in) :: this
       type(point_state), intent(in) :: st
       type(scalar)                  :: F(3, this % num_vertices)
     end function flux_value_interface

     !================================================================!
     ! Return the state jacobian dF/dq (3,nv,nv) at the point state.
     !================================================================!

     pure function flux_dq_interface(this, st) result(dF)
       import :: flux, point_state
       class(flux)      , intent(in) :: this
       type(point_state), intent(in) :: st
       type(scalar)                  :: dF(3, this % num_vertices, this % num_vertices)
     end function flux_dq_interface

     !================================================================!
     ! Return the gradient jacobian dF/d(grad q) (3,nv,3,nv) at the
     ! point state.
     !================================================================!

     pure function flux_dgradq_interface(this, st) result(dF)
       import :: flux, point_state
       class(flux)      , intent(in) :: this
       type(point_state), intent(in) :: st
       type(scalar)                  :: dF(3, this % num_vertices, 3, this % num_vertices)
     end function flux_dgradq_interface

  end interface

contains

  !===================================================================!
  ! The default design partial is zero; laws that depend on a design
  ! variable (for example, diffusion's conductivity) override it.
  !===================================================================!

  pure function flux_ddesign_zero(this, st, k) result(dF)
    class(flux)      , intent(in) :: this
    type(point_state), intent(in) :: st
    integer          , intent(in) :: k
    type(scalar)                  :: dF(3, this % num_vertices)
    dF = 0.0_dp
  end function flux_ddesign_zero

  !===================================================================!
  ! The normal diffusivity  keff = -nf^T (dF/dgradq) nf  is the
  ! conductivity the flux presents across a face with outward normal nf
  ! for variable ivar.
  !===================================================================!

  pure real(dp) function normal_diffusivity(this, st, nf, ivar)

    class(flux)      , intent(in) :: this
    type(point_state), intent(in) :: st
    real(dp)         , intent(in) :: nf(3)
    integer          , intent(in) :: ivar

    type(scalar) :: dFg(3, this % num_vertices, 3, this % num_vertices)

    dFg = this % dflux_dgradq(st)
    normal_diffusivity = -real(dot_product(nf, matmul(dFg(:,ivar,:,ivar), nf)), dp)

  end function normal_diffusivity

  !===================================================================!
  ! The normal speed  vn = nf^T (dF/dq)  is the advection speed the flux
  ! presents across a face with outward normal nf for variable ivar.
  !===================================================================!

  pure real(dp) function normal_speed(this, st, nf, ivar)

    class(flux)      , intent(in) :: this
    type(point_state), intent(in) :: st
    real(dp)         , intent(in) :: nf(3)
    integer          , intent(in) :: ivar

    type(scalar) :: dFq(3, this % num_vertices, this % num_vertices)

    dFq = this % dflux_dq(st)
    normal_speed = real(dot_product(nf, dFq(:,ivar,ivar)), dp)

  end function normal_speed

end module interface_flux
