!=====================================================================!
! Transient adjoint verification: march transient diffusion on a box
! (homogeneous dirichlet, unit source) with bdf order 2, take the
! isotropic conductivity kappa as the design variable and the time-
! integrated state energy  J = h sum_k 1/2 ||u_k||^2  as the function of
! interest. The transient discrete adjoint gradient dJ/dkappa (backward
! sweep over psi) must match a central finite difference of J(kappa)
! obtained by re-marching the forward problem.
!
! A nonzero exit (error stop) means the transient adjoint gradient is off.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program test_transient_adjoint

  use iso_fortran_env   , only : dp => REAL64
  use class_gmsh_loader , only : gmsh_loader
  use class_mesh        , only : mesh
  use class_assembler   , only : assembler
  use class_diffusion   , only : diffusion
  use class_state_energy, only : state_energy
  use class_bdf         , only : bdf
  use class_adjoint     , only : adjoint

  implicit none

  class(gmsh_loader), allocatable :: gl
  class(mesh)       , allocatable :: grid
  class(assembler)  , allocatable :: fvm
  type(bdf)                       :: ti
  type(state_energy)              :: func
  type(adjoint)                   :: adj
  real(dp), allocatable :: g_adj(:), g_fd(:)
  real(dp)              :: relerr
  integer               :: nfail

  nfail = 0

  ! build the assembler on a small box
  allocate(gl  , source = gmsh_loader("box-36.msh"))
  allocate(grid, source = mesh(gl))
  allocate(fvm , source = assembler(grid))

  ! transient diffusion: isotropic kappa = 2, unit source, faces at zero
  call fvm % set_equation(diffusion(2.0_dp, source = 1.0_dp))
  call fvm % set_dirichlet("front" , 0.0_dp)
  call fvm % set_dirichlet("back"  , 0.0_dp)
  call fvm % set_dirichlet("top"   , 0.0_dp)
  call fvm % set_dirichlet("bottom", 0.0_dp)
  call fvm % set_dirichlet("left"  , 0.0_dp)
  call fvm % set_dirichlet("right" , 0.0_dp)

  func = state_energy()

  ! march with bdf order 2 over t in [0,1], dt = 0.1 (the bandwidth ramps
  ! 1 -> 2, exercising the variable-bandwidth future-step coupling)
  ti  = bdf(fvm, 0.0_dp, 1.0_dp, 0.1_dp, 2)
  adj = adjoint(ti, func)

  ! transient adjoint gradient and the finite-difference reference
  call adj % eval_func_grad   (g_adj)
  call adj % eval_fd_func_grad(g_fd)

  ! export the state and adjoint-state trajectories: paraview (one .vtu
  ! per step) and gmsh (one .msh, two animated views over the time steps)
  call adj % write_solution("transient_adjoint")
  call adj % write_gmsh_solution("box-36.msh", "transient_adjoint.msh")

  relerr = abs(g_adj(1) - g_fd(1))/max(abs(g_fd(1)), tiny(1.0_dp))

  write(*,'(a)') " transient adjoint check - dJ/dkappa, J = h sum 1/2||u_k||^2"
  write(*,'(2x,a,es16.8)') "adjoint        : ", g_adj(1)
  write(*,'(2x,a,es16.8)') "finite-diff    : ", g_fd(1)
  write(*,'(2x,a,es16.8)') "relative error : ", relerr
  write(*,'(2x,a,i0,a)')   "wrote transient_adjoint_0000..", ti % num_steps-1, &
       & ".vtu + transient_adjoint.msh (state + adjoint)"

  if (relerr .gt. 1.0e-5_dp) nfail = nfail + 1

  write(*,*) "============================================="
  if (nfail .eq. 0) then
     write(*,*) "transient adjoint gradient matches finite difference"
  else
     write(*,*) "transient adjoint gradient check FAILED"
     error stop
  end if

end program test_transient_adjoint
