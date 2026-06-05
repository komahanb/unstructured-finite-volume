#include "scalar.fpp"

!=====================================================================!
! Discrete adjoint for the steady system  R(u, x) = 0,  giving the total
! derivative of a function of interest  J = f(u, x)  with respect to the
! design variables x:
!
!     dJ/dx = df/dx + psi^T dR/dx,     (dR/du)^T psi = -(df/du).
!
! The state jacobian action and its transpose are matrix-free (from the
! assembler), so the adjoint solve is a conjugate gradient on the
! transpose operator - no assembled matrix. eval_func_grad computes the
! adjoint gradient; eval_fd_func_grad re-solves at perturbed designs to
! verify it by central finite differences.
!
! Steady only for now (no time term); the transient backward sweep over
! psi will extend this.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module adjoint

  use iso_fortran_env    , only : dp => REAL64
  use interface_assembler, only : assembler
  use interface_function , only : functional
  use class_newton_solver, only : newton

  implicit none

  private
  public :: eval_func_grad, eval_fd_func_grad

  ! steady linearization coefficients [alpha, beta] = [1, 0]: dR/du only
  real(dp), parameter :: steady_coeff(2) = [1.0_dp, 0.0_dp]

contains

  !===================================================================!
  ! Adjoint total derivative dJ/dx for the steady problem
  !===================================================================!

  subroutine eval_func_grad(system, func, dJdx, adjoint_state)

    class(assembler) , intent(inout)                       :: system
    class(functional), intent(in)                          :: func
    real(dp)         , intent(out), allocatable            :: dJdx(:)
    type(scalar)     , intent(out), allocatable, optional  :: adjoint_state(:)

    type(scalar), allocatable :: dfdu(:), psi(:), rhs(:)
    integer                   :: n, ndv

    n   = system % get_num_state_vars()
    ndv = system % get_num_design_vars()

    ! 1. forward steady solve  R = 0  ->  u  (left in system % S)
    call solve_forward(system)

    ! 2. adjoint right-hand side  -df/du
    allocate(dfdu(n)); dfdu = 0.0d0
    call func % add_dfdu(system, dfdu)

    allocate(rhs(n)); rhs = -dfdu

    ! 3. adjoint solve  (dR/du)^T psi = -df/du
    call cg_solve_transpose(system, steady_coeff, rhs, psi)

    ! 4. total derivative  dJ/dx = df/dx + psi^T dR/dx
    allocate(dJdx(ndv)); dJdx = 0.0_dp
    call func   % add_dfdx(system, dJdx)
    call system % add_design_residual_transpose_product(dJdx, psi)

    ! hand back the adjoint state for post-processing if requested
    if (present(adjoint_state)) adjoint_state = psi

    deallocate(dfdu, rhs, psi)

  end subroutine eval_func_grad

  !===================================================================!
  ! Verification gradient: central finite differences of J over each
  ! design variable, re-solving the forward problem at each perturbation
  !===================================================================!

  subroutine eval_fd_func_grad(system, func, dJdx)

    class(assembler) , intent(inout)            :: system
    class(functional), intent(in)               :: func
    real(dp)         , intent(out), allocatable :: dJdx(:)

    real(dp)    , allocatable :: x0(:), x(:)
    type(scalar)              :: jp, jm
    real(dp)                  :: delta
    integer                   :: i, ndv

    ndv = system % get_num_design_vars()
    allocate(dJdx(ndv), x0(ndv), x(ndv))

    call system % get_design_vars(x0)

    do i = 1, ndv

       delta = 1.0e-6_dp*max(1.0_dp, abs(x0(i)))

       x = x0; x(i) = x0(i) + delta
       call system % set_design_vars(x)
       call solve_forward(system)
       call func % eval(system, jp)

       x = x0; x(i) = x0(i) - delta
       call system % set_design_vars(x)
       call solve_forward(system)
       call func % eval(system, jm)

       dJdx(i) = real(jp - jm, dp)/(2.0_dp*delta)

    end do

    ! restore the baseline design and state
    call system % set_design_vars(x0)
    call solve_forward(system)

    deallocate(x0, x)

  end subroutine eval_fd_func_grad

  !===================================================================!
  ! Solve the steady residual R(u) = 0 with the matrix-free newton (one
  ! linear step for the linear fvm); the converged u is left in S(:,1).
  !===================================================================!

  subroutine solve_forward(system)

    class(assembler), intent(inout) :: system

    type(scalar), allocatable :: U(:,:)
    type(newton)              :: nlsolver
    integer                   :: n, norder

    n      = system % get_num_state_vars()
    norder = system % get_differential_order()

    allocate(U(n, norder + 1)); U = 0.0d0
    call nlsolver % solve(system, steady_coeff, U)

    deallocate(U)

  end subroutine solve_forward

  !===================================================================!
  ! Matrix-free conjugate gradient on the TRANSPOSE jacobian operator
  !   B v = [scalar(i) dR/dU(i)]^T v
  ! supplied by add_jacobian_vector_product_transpose. For the steady
  ! diffusion problem B = -A^T is symmetric positive definite.
  !===================================================================!

  subroutine cg_solve_transpose(system, coeff, b, x)

    class(assembler), intent(in)               :: system
    type(scalar)    , intent(in)               :: coeff(:)
    type(scalar)    , intent(in)               :: b(:)
    type(scalar)    , intent(out), allocatable :: x(:)

    type(scalar), allocatable :: r(:), p(:), Bp(:)
    type(scalar)              :: rs_old, rs_new, alpha, pBp
    integer                   :: nvars, it, max_it

    nvars  = system % get_num_state_vars()
    max_it = nvars + 100

    allocate(x(nvars), r(nvars), p(nvars), Bp(nvars))

    x = 0.0d0
    r = b
    p = r

    rs_old = dot_product(r, r)

    cg: do it = 1, max_it

       if (vector_norm(r) .le. 1.0d-14) exit cg

       Bp = 0.0d0
       call system % add_jacobian_vector_product_transpose(Bp, p, coeff)

       pBp = dot_product(p, Bp)
       if (abs(pBp) .le. tiny(1.0_dp)) exit cg

       alpha = rs_old/pBp

       x = x + alpha*p
       r = r - alpha*Bp

       rs_new = dot_product(r, r)

       p      = r + (rs_new/rs_old)*p
       rs_old = rs_new

    end do cg

    deallocate(r, p, Bp)

  end subroutine cg_solve_transpose

  !===================================================================!
  ! 2-norm of a state vector (real or complex-step safe)
  !===================================================================!

  pure real(dp) function vector_norm(v)

    type(scalar), intent(in) :: v(:)

    vector_norm = sqrt(sum(abs(v)**2))

  end function vector_norm

end module adjoint
