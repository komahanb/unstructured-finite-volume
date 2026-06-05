#include "scalar.fpp"

!=====================================================================!
! Advection as a concrete law, expressed pointwise:
!
!     advection_flux            F(:,j) = v q(j)
!     advection_diffusion_flux  F(:,j) = v q(j) - K . grad q(:,j)
!
! with a constant velocity v (and, for the second, an isotropic
! conductivity K = kappa I). The advective term lives in dF/dq = v -
! the part the assembler picks up as the normal advection speed
! vn = v . n at a face (mirror of keff = n^T K n for diffusion). Unlike
! diffusion, the advective face flux uses the face VALUE of q, so the
! resulting operator is non-symmetric (central differencing gives a
! skew-symmetric contribution) - which is exactly why the nonsymmetric
! krylov solvers (class_gmres) exist.
!
! Pure steady advection is hyperbolic (ill-posed with all-dirichlet data),
! so advection_diffusion_flux is the well-posed operator to solve; the pure
! advection_flux is provided as the primitive (kappa = 0).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_advection_flux

  use iso_fortran_env  , only : dp => REAL64
  use interface_physics, only : flux, point_state

  implicit none

  private
  public :: advection_flux, advection_diffusion_flux

  !-------------------------------------------------------------------!
  ! Advective flux  F = v q
  !-------------------------------------------------------------------!

  type, extends(flux) :: advection_flux
     real(dp) :: v(3) = 0.0_dp        ! advection velocity
   contains
     procedure :: value        => advection_value
     procedure :: dflux_dq     => advection_dq
     procedure :: dflux_dgradq => advection_dgradq
  end type advection_flux

  interface advection_flux
     module procedure create_advection
  end interface advection_flux

  !-------------------------------------------------------------------!
  ! Advection + diffusion  F = v q - K grad q
  !-------------------------------------------------------------------!

  type, extends(flux) :: advection_diffusion_flux
     real(dp) :: v(3)      = 0.0_dp   ! advection velocity
     real(dp) :: kmat(3,3) = 0.0_dp   ! conductivity tensor (kappa I)
   contains
     procedure :: value        => advdiff_value
     procedure :: dflux_dq     => advdiff_dq
     procedure :: dflux_dgradq => advdiff_dgradq
  end type advection_diffusion_flux

  interface advection_diffusion_flux
     module procedure create_advdiff
  end interface advection_diffusion_flux

contains

  !===================================================================!
  ! Constructors
  !===================================================================!

  pure type(advection_flux) function create_advection(velocity, nvars) result(this)
    real(dp), intent(in)           :: velocity(3)
    integer , intent(in), optional :: nvars
    this % num_components = 1
    if (present(nvars)) this % num_components = nvars
    this % v = velocity
  end function create_advection

  pure type(advection_diffusion_flux) function create_advdiff(velocity, kappa, nvars) result(this)
    real(dp), intent(in)           :: velocity(3)
    real(dp), intent(in)           :: kappa
    integer , intent(in), optional :: nvars
    integer :: i
    this % num_components = 1
    if (present(nvars)) this % num_components = nvars
    this % v    = velocity
    this % kmat = 0.0_dp
    do i = 1, 3
       this % kmat(i,i) = kappa
    end do
  end function create_advdiff

  !===================================================================!
  ! advection_flux:  F(:,j) = v q(j);  dF/dq(:,i,i) = v;  dF/dgradq = 0
  !===================================================================!

  pure function advection_value(this, st) result(F)
    class(advection_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    type(scalar)                      :: F(3, this % num_components)
    integer :: j
    do j = 1, this % num_components
       F(:,j) = this % v * st % q(j)
    end do
  end function advection_value

  pure function advection_dq(this, st) result(dF)
    class(advection_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    type(scalar)                      :: dF(3, this % num_components, this % num_components)
    integer :: i
    dF = 0.0_dp
    do i = 1, this % num_components
       dF(:, i, i) = this % v
    end do
  end function advection_dq

  pure function advection_dgradq(this, st) result(dF)
    class(advection_flux), intent(in) :: this
    type(point_state)    , intent(in) :: st
    type(scalar)                      :: dF(3, this % num_components, 3, this % num_components)
    dF = 0.0_dp
  end function advection_dgradq

  !===================================================================!
  ! advection_diffusion_flux:  F(:,j) = v q(j) - K grad q(:,j)
  !   dF/dq(:,i,i) = v ;  dF/d(grad q)(:,i,:,i) = -K
  !===================================================================!

  pure function advdiff_value(this, st) result(F)
    class(advection_diffusion_flux), intent(in) :: this
    type(point_state)              , intent(in) :: st
    type(scalar)                                :: F(3, this % num_components)
    integer :: j
    do j = 1, this % num_components
       F(:,j) = this % v * st % q(j) - matmul(this % kmat, st % gradq(:,j))
    end do
  end function advdiff_value

  pure function advdiff_dq(this, st) result(dF)
    class(advection_diffusion_flux), intent(in) :: this
    type(point_state)              , intent(in) :: st
    type(scalar)                                :: dF(3, this % num_components, this % num_components)
    integer :: i
    dF = 0.0_dp
    do i = 1, this % num_components
       dF(:, i, i) = this % v
    end do
  end function advdiff_dq

  pure function advdiff_dgradq(this, st) result(dF)
    class(advection_diffusion_flux), intent(in) :: this
    type(point_state)              , intent(in) :: st
    type(scalar)                                :: dF(3, this % num_components, 3, this % num_components)
    integer :: i
    dF = 0.0_dp
    do i = 1, this % num_components
       dF(:, i, :, i) = -this % kmat
    end do
  end function advdiff_dgradq

end module class_advection_flux
