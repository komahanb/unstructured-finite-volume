!=====================================================================!
! Steady diffusion  -div(K grad phi) = s  with a constant (possibly
! anisotropic) tensor K and a constant source s. The first concrete
! equation - enough to recover the old laplace solver (K=I, s=0) and to
! drive anisotropic/poisson cases.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_diffusion

  use iso_fortran_env  , only : dp => REAL64
  use interface_equation, only : equation

  implicit none

  private
  public :: diffusion

  type, extends(equation) :: diffusion

     real(dp) :: kmat(3,3)   ! diffusion tensor
     real(dp) :: src         ! volumetric source

   contains

     procedure :: diffusion_tensor
     procedure :: source

  end type diffusion

  ! isotropic: diffusion(kappa [, source]) ; anisotropic: diffusion(K [, source])
  interface diffusion
     module procedure create_isotropic
     module procedure create_tensor
  end interface diffusion

contains

  pure type(diffusion) function create_isotropic(kappa, source, nvars) result(this)
    real(dp), intent(in)           :: kappa
    real(dp), intent(in), optional :: source
    integer , intent(in), optional :: nvars
    integer :: i
    this % num_variables = 1
    if (present(nvars)) this % num_variables = nvars
    this % kmat = 0.0_dp
    do i = 1, 3
       this % kmat(i,i) = kappa
    end do
    this % src = 0.0_dp
    if (present(source)) this % src = source
  end function create_isotropic

  pure type(diffusion) function create_tensor(K, source, nvars) result(this)
    real(dp), intent(in)           :: K(3,3)
    real(dp), intent(in), optional :: source
    integer , intent(in), optional :: nvars
    this % num_variables = 1
    if (present(nvars)) this % num_variables = nvars
    this % kmat = K
    this % src = 0.0_dp
    if (present(source)) this % src = source
  end function create_tensor

  pure function diffusion_tensor(this, x) result(K)
    class(diffusion), intent(in) :: this
    real(dp)        , intent(in) :: x(3)
    real(dp)                     :: K(3,3)
    K = this % kmat   ! constant in space (ignore x)
  end function diffusion_tensor

  pure function source(this, x) result(s)
    class(diffusion), intent(in) :: this
    real(dp)        , intent(in) :: x(3)
    real(dp)                     :: s
    s = this % src
  end function source

end module class_diffusion
