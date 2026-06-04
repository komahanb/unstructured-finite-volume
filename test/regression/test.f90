!=====================================================================!
! Regression / verification suite. Each check is an analytic oracle for
! the finite volume solver, so it needs no stored golden numbers:
!
!   1. constant dirichlet on the sphere  -> phi = 5 everywhere
!   2. kappa scaling                     -> doubling kappa halves phi
!   3. transient (t -> inf)              -> the steady solution
!   4. operator split                    -> A = L + U + D
!   5. neumann insulated box             -> linear profile, exact mean
!   6. bdf-1 integrator                  -> matches legacy backward-euler
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
  use class_time_integrator    , only : time_integrator
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
  call check_bdf_matches_be(nfail)

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
    type(time_integrator)                  :: ti
    real(dp), allocatable :: xs(:), xt(:), phi0(:)

    ! steady
    call make("../box-3.msh", fsteady)
    call box_bc(fsteady)
    allocate(cg, source = conjugate_gradient(fsteady, 500, tol, 0))
    call cg % solve(xs)

    ! transient from zero, marched far
    call make("../box-3.msh", ftrans)
    call box_bc(ftrans)
    ti = time_integrator(ftrans, 0.0_dp, 200.0_dp, 10.0_dp, 500, tol)
    allocate(phi0(ftrans % num_state_vars))
    phi0 = 0.0_dp
    call ti % integrate(phi0, xt)

    call report("transient (t -> inf) -> steady", &
         & maxval(abs(xt - xs)) .lt. 1.0e-8_dp, nfail)

  end subroutine check_transient_steady

  !===================================================================!
  ! The order-1 bdf integrator (M*udot - A*u + b = 0, marched by newton)
  ! must reproduce the legacy backward-euler march step for step. Same
  ! mesh, bcs, dt and step count, both from a zero initial field.
  !===================================================================!

  subroutine check_bdf_matches_be(nfail)

    integer, intent(inout) :: nfail

    class(assembler)      , allocatable :: fbe, fbdf
    type(time_integrator)               :: ti
    type(bdf)                           :: bi
    real(dp), allocatable :: xbe(:), xbdf(:), phi0(:)

    ! legacy backward-euler reference
    call make("../box-3.msh", fbe)
    call box_bc(fbe)
    ti = time_integrator(fbe, 0.0_dp, 200.0_dp, 10.0_dp, 500, tol)
    allocate(phi0(fbe % num_state_vars)); phi0 = 0.0_dp
    call ti % integrate(phi0, xbe)

    ! same march through the general bdf integrator at order 1 (transient
    ! flag stays off: the integrator, not the assembler, owns the time term)
    call make("../box-3.msh", fbdf)
    call box_bc(fbdf)
    bi = bdf(fbdf, 0.0_dp, 200.0_dp, 10.0_dp, max_order = 1)
    call bi % solve()
    xbdf = real(bi % U(bi % num_steps, :, 1), dp)

    call report("bdf-1 integrator -> matches backward-euler", &
         & maxval(abs(xbdf - xbe)) .lt. 1.0e-8_dp, nfail)

  end subroutine check_bdf_matches_be

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
