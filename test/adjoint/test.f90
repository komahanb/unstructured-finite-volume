#include "scalar.fpp"

!=====================================================================!
! Adjoint verification: steady poisson on a box (homogeneous dirichlet,
! unit source) with the isotropic conductivity kappa as the lone design
! variable and the state energy J = 1/2 ||u||^2 as the function of
! interest. The discrete adjoint gradient dJ/dkappa must match a central
! finite-difference of J(kappa) obtained by re-solving the forward
! problem.
!
! A nonzero exit (error stop) means the adjoint gradient is off.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_adjoint

  use iso_fortran_env   , only : dp => REAL64
  use class_gmsh_loader , only : gmsh_loader
  use class_mesh        , only : mesh
  use class_assembler   , only : assembler
  use class_newton_solver, only : newton
  use class_diffusion_flux  , only : diffusion_flux, constant_source
  use class_state_energy, only : state_energy

  implicit none

  class(gmsh_loader), allocatable :: gl
  class(mesh)       , allocatable :: grid
  class(assembler)  , allocatable :: fvm
  type(newton)                    :: nl
  type(state_energy)              :: func
  real(dp)    , allocatable :: g_adj(:), g_fd(:), fields(:,:)
  type(scalar), allocatable :: psi(:)
  real(dp)                  :: relerr
  integer                   :: nfail

  nfail = 0

  ! build the assembler on a small box
  allocate(gl  , source = gmsh_loader("box-36.msh"))
  allocate(grid, source = mesh(gl))
  allocate(fvm , source = assembler(grid))

  ! steady poisson: isotropic kappa = 2, unit volumetric source, all six
  ! faces held at zero (homogeneous dirichlet)
  call fvm % set_equation(diffusion_flux(2.0_dp), constant_source(1.0_dp))
  call fvm % set_dirichlet("front" , 0.0_dp)
  call fvm % set_dirichlet("back"  , 0.0_dp)
  call fvm % set_dirichlet("top"   , 0.0_dp)
  call fvm % set_dirichlet("bottom", 0.0_dp)
  call fvm % set_dirichlet("left"  , 0.0_dp)
  call fvm % set_dirichlet("right" , 0.0_dp)

  func = state_energy()

  ! adjoint gradient (also hand back the adjoint state) and the fd reference
  call nl % eval_func_grad   (fvm, func, g_adj, psi)
  call nl % eval_fd_func_grad(fvm, func, g_fd)

  ! export the state and adjoint state for post-processing (paraview .vtu
  ! and gmsh .msh)
  allocate(fields(fvm % num_state_vars, 2))
  fields(:,1) = real(fvm % S(:,1), dp)   ! state u
  fields(:,2) = real(psi         , dp)   ! adjoint state psi
  call fvm % write_solution_fields("steady_adjoint.vtu", fields, &
       & [character(len=7) :: "state", "adjoint"])
  call fvm % write_gmsh_series("box-36.msh", "steady_adjoint.msh", &
       & reshape(fields, [fvm % num_state_vars, 2, 1]), &
       & [character(len=7) :: "state", "adjoint"], [0.0_dp])

  relerr = abs(g_adj(1) - g_fd(1))/max(abs(g_fd(1)), tiny(1.0_dp))

  write(*,'(a)') " adjoint gradient check - dJ/dkappa, J = 1/2||u||^2"
  write(*,'(2x,a,es16.8)') "adjoint        : ", g_adj(1)
  write(*,'(2x,a,es16.8)') "finite-diff    : ", g_fd(1)
  write(*,'(2x,a,es16.8)') "relative error : ", relerr
  write(*,'(2x,a)')        "wrote steady_adjoint.vtu / .msh (state + adjoint)"

  if (relerr .gt. 1.0e-5_dp) nfail = nfail + 1

  write(*,*) "============================================="
  if (nfail .eq. 0) then
     write(*,*) "adjoint gradient matches finite difference"
  else
     write(*,*) "adjoint gradient check FAILED"
     error stop
  end if

end program test_adjoint
