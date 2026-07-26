#include "scalar.fpp"

!=====================================================================!
! Diffusion is the first concrete law, expressed pointwise:
!
!     flux   F(:,j) = -K . grad q(:,j)        (Fourier / Fick)
!     source S(:,j) = s                        (constant volumetric)
!
! Both use a constant (here isotropic) conductivity tensor K. This is the
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
  ! The diffusive flux is  F = -K grad q.
  !-------------------------------------------------------------------!

  type, extends(flux) :: diffusion_flux
     real(dp) :: kmat(3,3) = 0.0_dp   ! The conductivity tensor.
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
  ! The constant volumetric source is  S = s.
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
  ! Construct a diffusion flux with isotropic conductivity kappa, so
  ! that K = kappa I.
  !===================================================================!

  pure type(diffusion_flux) function create_diffusion_flux(kappa, nvars) result(this)
    real(dp), intent(in)           :: kappa
    integer , intent(in), optional :: nvars
    integer :: i
    ! One vertex per variable and no edges: today's laws couple nothing.
    if (present(nvars)) then
       call this % form_coupling(nvars)
    else
       call this % form_coupling(1)
    end if
    this % kmat = 0.0_dp
    do i = 1, 3
       this % kmat(i,i) = kappa
    end do
  end function create_diffusion_flux

  !===================================================================!
  ! Construct a constant volumetric source of strength s.
  !===================================================================!

  pure type(constant_source) function create_constant_source(s, nvars) result(this)
    real(dp), intent(in)           :: s
    integer , intent(in), optional :: nvars
    ! One vertex per variable and no edges: today's laws couple nothing.
    if (present(nvars)) then
       call this % form_coupling(nvars)
    else
       call this % form_coupling(1)
    end if
    this % s = s
  end function create_constant_source

  !===================================================================!
  ! The flux value is  F(:,j) = -K grad q(:,j).
  !===================================================================!

  pure function diffusion_flux_value(this, st) result(F)
    class(diffusion_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    type(scalar)                      :: F(3, this % num_vertices)
    integer :: j
    do j = 1, this % num_vertices
       F(:,j) = -matmul(this % kmat, st % gradq(:,j))
    end do
  end function diffusion_flux_value

  !===================================================================!
  ! The partial dF/dq is zero: the flux has no direct q dependence.
  !===================================================================!

  pure function diffusion_flux_dq(this, st) result(dF)
    class(diffusion_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    type(scalar)                      :: dF(3, this % num_vertices, this % num_vertices)
    dF = 0.0_dp
  end function diffusion_flux_dq

  !===================================================================!
  ! The gradient jacobian is  dF(:,i)/d(grad q(:,j)) = -K delta_ij;
  ! the variables are decoupled.
  !===================================================================!

  pure function diffusion_flux_dgradq(this, st) result(dF)
    class(diffusion_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    type(scalar)                      :: dF(3, this % num_vertices, 3, this % num_vertices)
    integer :: i
    dF = 0.0_dp
    do i = 1, this % num_vertices
       dF(:, i, :, i) = -this % kmat
    end do
  end function diffusion_flux_dgradq

  !===================================================================!
  ! The design partial is  dF/dkappa = -grad q = F/kappa  for the
  ! isotropic conductivity design variable.
  !===================================================================!

  pure function diffusion_flux_ddesign(this, st, k) result(dF)
    class(diffusion_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    integer              , intent(in) :: k
    type(scalar)                      :: dF(3, this % num_vertices)
    integer :: j
    do j = 1, this % num_vertices
       dF(:,j) = -st % gradq(:,j)
    end do
  end function diffusion_flux_ddesign

  !===================================================================!
  ! There is one design variable: the isotropic conductivity kappa,
  ! with K = kappa I.
  !===================================================================!

  pure integer function diffusion_num_design_vars(this) result(n)
    class(diffusion_flux), intent(in) :: this
    n = 1
  end function diffusion_num_design_vars

  !===================================================================!
  ! Set the isotropic conductivity from the design vector x.
  !===================================================================!

  pure subroutine diffusion_set_design_vars(this, x)
    class(diffusion_flux), intent(inout) :: this
    real(dp)             , intent(in)    :: x(:)
    integer :: i
    this % kmat = 0.0_dp
    do i = 1, 3
       this % kmat(i,i) = x(1)
    end do
  end subroutine diffusion_set_design_vars

  !===================================================================!
  ! Report the isotropic conductivity into the design vector x.
  !===================================================================!

  pure subroutine diffusion_get_design_vars(this, x)
    class(diffusion_flux), intent(in)  :: this
    real(dp)             , intent(out) :: x(:)
    x(1) = this % kmat(1,1)
  end subroutine diffusion_get_design_vars

  !===================================================================!
  ! The constant source returns its stored strength s.
  !===================================================================!

  pure function constant_source_value(this, st) result(S)
    class(constant_source), intent(in) :: this
    type(point_state)     , intent(in) :: st
    type(scalar)                       :: S(this % num_vertices)
    S = this % s
  end function constant_source_value

end module class_diffusion_flux
