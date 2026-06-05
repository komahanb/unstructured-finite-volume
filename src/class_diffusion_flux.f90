#include "scalar.fpp"

!=====================================================================!
! Diffusion as the first concrete law, expressed pointwise:
!
!     flux   F(:,j) = -K . grad q(:,j)        (Fourier / Fick)
!     source S(:,j) = s                        (constant volumetric)
!
! with a constant (here isotropic) conductivity tensor K. This is the
! law-agnostic restatement of the old class_diffusion: a finite-volume
! assembler dotting F with a face normal recovers exactly the normal
! diffusive flux keff = n^T K n, and dF/dgradq = -K gives that keff to
! the matrix-free jacobian. The single design variable is the isotropic
! conductivity kappa (K = kappa I), with dF/dkappa = -grad q = F/kappa.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_diffusion_flux

  use iso_fortran_env , only : dp => REAL64
  use interface_flux   , only : flux
  use interface_physics, only : source, point_state

  implicit none

  private
  public :: diffusion_flux, constant_source

  !-------------------------------------------------------------------!
  ! Diffusive flux  F = -K grad q
  !-------------------------------------------------------------------!

  type, extends(flux) :: diffusion_flux
     real(dp) :: kmat(3,3) = 0.0_dp   ! conductivity tensor
   contains
     procedure :: value          => diffusion_flux_value
     procedure :: dflux_dq       => diffusion_flux_dq
     procedure :: dflux_dgradq   => diffusion_flux_dgradq
     procedure :: dflux_ddesign  => diffusion_flux_ddesign
     procedure :: num_design_vars => diffusion_num_design_vars
     procedure :: set_design_vars => diffusion_set_design_vars
     procedure :: get_design_vars => diffusion_get_design_vars
  end type diffusion_flux

  interface diffusion_flux
     module procedure create_diffusion_flux
  end interface diffusion_flux

  !-------------------------------------------------------------------!
  ! Constant volumetric source  S = s
  !-------------------------------------------------------------------!

  type, extends(source) :: constant_source
     real(dp) :: s = 0.0_dp
   contains
     procedure :: value => constant_source_value
  end type constant_source

  interface constant_source
     module procedure create_constant_source
  end interface constant_source

contains

  !===================================================================!
  ! Constructors
  !===================================================================!

  pure type(diffusion_flux) function create_diffusion_flux(kappa, nvars) result(this)
    real(dp), intent(in)           :: kappa
    integer , intent(in), optional :: nvars
    integer :: i
    this % num_components = 1
    if (present(nvars)) this % num_components = nvars
    this % kmat = 0.0_dp
    do i = 1, 3
       this % kmat(i,i) = kappa
    end do
  end function create_diffusion_flux

  pure type(constant_source) function create_constant_source(s, nvars) result(this)
    real(dp), intent(in)           :: s
    integer , intent(in), optional :: nvars
    this % num_components = 1
    if (present(nvars)) this % num_components = nvars
    this % s = s
  end function create_constant_source

  !===================================================================!
  ! Flux value  F(:,j) = -K grad q(:,j)
  !===================================================================!

  pure function diffusion_flux_value(this, st) result(F)
    class(diffusion_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    type(scalar)                      :: F(3, this % num_components)
    integer :: j
    do j = 1, this % num_components
       F(:,j) = -matmul(this % kmat, st % gradq(:,j))
    end do
  end function diffusion_flux_value

  !===================================================================!
  ! dF/dq = 0 (the diffusive flux has no direct q dependence)
  !===================================================================!

  pure function diffusion_flux_dq(this, st) result(dF)
    class(diffusion_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    type(scalar)                      :: dF(3, this % num_components, this % num_components)
    dF = 0.0_dp
  end function diffusion_flux_dq

  !===================================================================!
  ! dF(:,i)/d(grad q(:,j)) = -K delta_ij  (variables decoupled)
  !===================================================================!

  pure function diffusion_flux_dgradq(this, st) result(dF)
    class(diffusion_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    type(scalar)                      :: dF(3, this % num_components, 3, this % num_components)
    integer :: i
    dF = 0.0_dp
    do i = 1, this % num_components
       dF(:, i, :, i) = -this % kmat
    end do
  end function diffusion_flux_dgradq

  !===================================================================!
  ! dF/dkappa = -grad q = F/kappa  (isotropic conductivity design var)
  !===================================================================!

  pure function diffusion_flux_ddesign(this, st, k) result(dF)
    class(diffusion_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    integer              , intent(in) :: k
    type(scalar)                      :: dF(3, this % num_components)
    integer :: j
    do j = 1, this % num_components
       dF(:,j) = -st % gradq(:,j)
    end do
  end function diffusion_flux_ddesign

  !===================================================================!
  ! One design variable: the isotropic conductivity kappa (K = kappa I)
  !===================================================================!

  pure integer function diffusion_num_design_vars(this) result(n)
    class(diffusion_flux), intent(in) :: this
    n = 1
  end function diffusion_num_design_vars

  subroutine diffusion_set_design_vars(this, x)
    class(diffusion_flux), intent(inout) :: this
    real(dp)             , intent(in)    :: x(:)
    integer :: i
    this % kmat = 0.0_dp
    do i = 1, 3
       this % kmat(i,i) = x(1)
    end do
  end subroutine diffusion_set_design_vars

  subroutine diffusion_get_design_vars(this, x)
    class(diffusion_flux), intent(in)  :: this
    real(dp)             , intent(out) :: x(:)
    x(1) = this % kmat(1,1)
  end subroutine diffusion_get_design_vars

  !===================================================================!
  ! Constant source value
  !===================================================================!

  pure function constant_source_value(this, st) result(S)
    class(constant_source), intent(in) :: this
    type(point_state)     , intent(in) :: st
    type(scalar)                       :: S(this % num_components)
    S = this % s
  end function constant_source_value

end module class_diffusion_flux
