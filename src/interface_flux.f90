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

end module interface_flux
