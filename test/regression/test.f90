!=====================================================================!
! Regression / verification suite. Each check is an analytic oracle for
! the finite volume solver, so it needs no stored golden numbers:
!
!   1. constant dirichlet on the sphere  -> phi = 5 everywhere
!   2. kappa scaling                     -> doubling kappa halves phi
!   3. transient (t -> inf)              -> the steady solution
!   4. operator split                    -> A = L + U + D
!   5. neumann insulated box             -> linear profile, exact mean
!   6. solver wrappers (cgnr, dist-cg)   -> match plain cg
!
! A nonzero exit (error stop) means a check failed.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program regression

  use iso_fortran_env          , only : dp => REAL64
  use class_gmsh_loader        , only : gmsh_loader
  use class_mesh               , only : mesh
  use class_assembler          , only : assembler
  use class_diffusion_flux  , only : diffusion_flux, constant_source
  use class_conjugate_gradient , only : conjugate_gradient
  use class_normal_cg          , only : normal_cg, CGNR_METHOD
  use class_distributed_cg     , only : distributed_cg_solver
  use class_bdf                , only : bdf

  implicit none

  integer  :: nfail
  real(dp) :: tol

  nfail = 0
  tol   = 100.0_dp*epsilon(1.0_dp)

  call check_sphere_constant(nfail)
  call check_kappa_scaling(nfail)
  call check_transient_steady(nfail)
  call check_operator_split(nfail)
  call check_neumann_insulated(nfail)
  call check_solver_wrappers(nfail)

  write(*,*) "============================================="
  if (nfail .eq. 0) then
     write(*,*) "all regression checks passed"
  else
     write(*,*) nfail, " regression checks FAILED"
     error stop
  end if

contains

  !===================================================================!
  ! Build a fresh assembler on a mesh file
  !===================================================================!

  subroutine make(meshfile, fvm)

    character(len=*)             , intent(in)  :: meshfile
    class(assembler), allocatable, intent(out) :: fvm

    class(gmsh_loader), allocatable :: gl
    class(mesh)       , allocatable :: grid

    allocate(gl, source = gmsh_loader(meshfile))
    allocate(grid, source = mesh(gl))
    allocate(fvm, source = assembler(grid))

  end subroutine make

  !===================================================================!
  ! Report the outcome of one check
  !===================================================================!

  subroutine report(name, ok, nfail)

    character(len=*), intent(in)    :: name
    logical         , intent(in)    :: ok
    integer         , intent(inout) :: nfail

    if (ok) then
       write(*,'(1x,a,a)') "PASS : ", name
    else
       write(*,'(1x,a,a)') "FAIL : ", name
       nfail = nfail + 1
    end if

  end subroutine report

  !===================================================================!
  ! Constant dirichlet (=5) on every boundary of the sphere must give a
  ! uniform field phi = 5 (laplace, no source).
  !===================================================================!

  subroutine check_sphere_constant(nfail)

    integer, intent(inout) :: nfail

    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg
    real(dp)                 , allocatable :: x(:)

    call make("../sphere.msh", fvm)
    call fvm % set_dirichlet("boundary", 5.0_dp)

    allocate(cg, source = conjugate_gradient(fvm, 500, tol, 0))
    call cg % solve(x)

    call report("sphere constant dirichlet -> phi = 5", &
         & maxval(abs(x - 5.0_dp)) .lt. 1.0e-10_dp, nfail)

  end subroutine check_sphere_constant

  !===================================================================!
  ! Doubling an isotropic kappa must halve the solution (linear poisson
  ! with fixed source and homogeneous dirichlet).
  !===================================================================!

  subroutine check_kappa_scaling(nfail)

    integer, intent(inout) :: nfail

    real(dp), allocatable :: x1(:), x2(:)

    call poisson("../box-36.msh", 1.0_dp, x1)
    call poisson("../box-36.msh", 2.0_dp, x2)

    call report("kappa scaling -> phi(2k) = phi(k)/2", &
         & maxval(abs(x2 - 0.5_dp*x1)) .lt. 1.0e-10_dp, nfail)

  end subroutine check_kappa_scaling

  ! Helper: homogeneous-dirichlet poisson with unit source and given kappa
  subroutine poisson(meshfile, kappa, x)

    character(len=*)     , intent(in)  :: meshfile
    real(dp)             , intent(in)  :: kappa
    real(dp), allocatable, intent(out) :: x(:)

    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg

    call make(meshfile, fvm)
    call fvm % set_equation(diffusion_flux(kappa), constant_source(1.0_dp))

    allocate(cg, source = conjugate_gradient(fvm, 500, tol, 0))
    call cg % solve(x)

  end subroutine poisson

  !===================================================================!
  ! Backward-euler marched to large time must reach the steady solution
  !===================================================================!

  subroutine check_transient_steady(nfail)

    integer, intent(inout) :: nfail

    class(assembler)         , allocatable :: fsteady, ftrans
    class(conjugate_gradient), allocatable :: cg
    type(bdf)                              :: ti
    real(dp), allocatable :: xs(:), xt(:)

    ! steady
    call make("../box-3.msh", fsteady)
    call box_bc(fsteady)
    allocate(cg, source = conjugate_gradient(fsteady, 500, tol, 0))
    call cg % solve(xs)

    ! transient from zero, marched far (bdf order 1 = backward euler)
    call make("../box-3.msh", ftrans)
    call box_bc(ftrans)
    ti = bdf(ftrans, 0.0_dp, 200.0_dp, 10.0_dp, max_order = 1)
    call ti % solve()
    xt = real(ti % U(ti % num_steps, :, 1), dp)

    call report("transient (t -> inf) -> steady", &
         & maxval(abs(xt - xs)) .lt. 1.0e-8_dp, nfail)

  end subroutine check_transient_steady

  !===================================================================!
  ! The normal_cg (cgnr) and distributed_cg linear_solver wrappers must
  ! solve the same dirichlet problem as plain CG.
  !===================================================================!

  subroutine check_solver_wrappers(nfail)

    integer, intent(inout) :: nfail

    class(assembler)           , allocatable :: f1, f2, f3
    class(conjugate_gradient)  , allocatable :: cg
    type(normal_cg)                          :: ncg
    type(distributed_cg_solver)              :: dcg
    real(dp), allocatable :: xref(:), xn(:), xd(:)

    ! reference: plain CG
    call make("../box-3.msh", f1); call box_bc(f1)
    allocate(cg, source = conjugate_gradient(f1, 5000, tol, 0))
    call cg % solve(xref)

    ! normal_cg (cgnr): CG on the normal equations (kappa^2 -> looser tol)
    call make("../box-3.msh", f2); call box_bc(f2)
    ncg = normal_cg(FVAssembler=f2, max_it=50000, max_tol=1.0e-10_dp, &
         & method=CGNR_METHOD, print_level=0)
    call ncg % solve(xn)

    ! distributed_cg wrapper (serial build -> plain CG)
    call make("../box-3.msh", f3); call box_bc(f3)
    dcg = distributed_cg_solver(FVAssembler=f3, max_it=5000, max_tol=tol, print_level=0)
    call dcg % solve(xd)

    call report("cgnr wrapper      -> matches cg", &
         & maxval(abs(xn - xref)) .lt. 1.0e-5_dp, nfail)
    call report("distributed_cg    -> matches cg", &
         & maxval(abs(xd - xref)) .lt. 1.0e-8_dp, nfail)

  end subroutine check_solver_wrappers

  !===================================================================!
  ! The assembled operator must split as A = L + U + D
  !===================================================================!

  subroutine check_operator_split(nfail)

    integer, intent(inout) :: nfail

    class(assembler), allocatable :: fvm
    real(dp), allocatable :: A(:,:), L(:,:), U(:,:), D(:,:)

    call make("../box-36.msh", fvm)

    call fvm % get_jacobian(A)
    call fvm % get_jacobian(L, filter = fvm % LOWER_TRIANGLE)
    call fvm % get_jacobian(U, filter = fvm % UPPER_TRIANGLE)
    call fvm % get_jacobian(D, filter = fvm % DIAGONAL)

    call report("operator split A = L + U + D", &
         & maxval(abs(A - L - U - D)) .lt. tiny(1.0_dp), nfail)

  end subroutine check_operator_split

  !===================================================================!
  ! Box with front=5, back=0 and the four sides insulated (neumann 0)
  ! gives a 1d conduction: the cell-average must be the midpoint 2.5 and
  ! the field must stay within the end values.
  !===================================================================!

  subroutine check_neumann_insulated(nfail)

    integer, intent(inout) :: nfail

    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg
    real(dp)                 , allocatable :: x(:)
    logical :: ok

    call make("../box-3.msh", fvm)
    call fvm % set_dirichlet("front", 5.0_dp)
    call fvm % set_dirichlet("back" , 0.0_dp)
    call fvm % set_neumann  ("top"   , 0.0_dp)
    call fvm % set_neumann  ("bottom", 0.0_dp)
    call fvm % set_neumann  ("left"  , 0.0_dp)
    call fvm % set_neumann  ("right" , 0.0_dp)

    allocate(cg, source = conjugate_gradient(fvm, 500, tol, 0))
    call cg % solve(x)

    ok =       (abs(sum(x)/size(x) - 2.5_dp) .lt. 1.0e-8_dp) &
         .and. (minval(x) .gt. 0.0_dp) &
         .and. (maxval(x) .lt. 5.0_dp)

    call report("neumann insulated box -> mean 2.5, bounded", ok, nfail)

  end subroutine check_neumann_insulated

  !===================================================================!
  ! The box-3 dirichlet boundary values used by several checks
  !===================================================================!

  subroutine box_bc(fvm)

    class(assembler), intent(inout) :: fvm

    call fvm % set_dirichlet("front" , 5.0_dp)
    call fvm % set_dirichlet("bottom", 10.0_dp)
    call fvm % set_dirichlet("right" , 15.0_dp)
    call fvm % set_dirichlet("top"   , 0.0_dp)
    call fvm % set_dirichlet("left"  , 0.0_dp)
    call fvm % set_dirichlet("back"  , 0.0_dp)

  end subroutine box_bc

end program regression
