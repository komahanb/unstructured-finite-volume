#include "scalar.fpp"

!=====================================================================!
! Condensed Newton solver for one implicit step, R(U) = 0, where
! U = [q, qdot, qddot, ...] is the assembler's state. The linearised
! system is matrix-free:
!
!     J dq = -R,     J v = sum_n coeff(n) dR/dU(n) v
!
! supplied by the assembler's add_jacobian_vector_product, and solved
! with a conjugate-gradient inner iteration. The single increment dq then
! updates every derivative order via U(:,n) += coeff(n) dq.
!
! Concrete newton extends the nonlinear_solver interface; the BDF marcher
! (class_bdf) and the adjoint forward solve call its solve. Everything is
! vector-level (no element indexing) so flat arrays can become distributed
! vectors later.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_newton_solver

  use iso_fortran_env           , only : dp => REAL64
  use interface_assembler       , only : assembler
  use interface_nonlinear_solver, only : nonlinear_solver
  use class_differential_state  , only : differential_state
  use interface_function        , only : functional
  use class_conjugate_gradient  , only : conjugate_gradient
  use module_solve_mode         , only : FORWARD, REVERSE

  implicit none

  private
  public :: newton

  ! steady linearization coefficients [alpha, beta] = [1, 0]: dR/du only
  real(dp), parameter :: steady_coeff(2) = [1.0_dp, 0.0_dp]

  !-------------------------------------------------------------------!
  ! Concrete Newton solver. Stopping criteria are members (defaults
  ! reproduce the previous nonlinear_marching constants).
  !-------------------------------------------------------------------!

  type, extends(nonlinear_solver) :: newton

     real(dp) :: abs_tol = 1.0d-12
     real(dp) :: rel_tol = 1.0d-11
     ! max_it (default 25), max_tol, print_level inherited from nonlinear_solver

   contains

     procedure :: solve
     procedure, private :: prepare_inner_solver

     ! adjoint total derivative dJ/dx (+ fd verify); the steady forward
     ! solve is the inherited march at the steady linearization
     procedure :: eval_func_grad
     procedure :: eval_fd_func_grad

  end type newton

contains

  !===================================================================!
  ! Drive R(U) -> 0 for the newest state U = (nvars, order+1)
  !===================================================================!

  impure subroutine solve(this, system, coeff, U)

    class(newton)   , intent(inout) :: this
    class(assembler), intent(inout) :: system
    type(scalar)    , intent(in)    :: coeff(:)
    type(scalar)    , intent(inout) :: U(:,:)

    type(scalar), allocatable :: res(:)
    real(dp)    , allocatable :: dq(:)
    type(differential_state)  :: s
    real(dp)                  :: r0, rnorm
    integer                   :: iter

    ! point the held linear solver at this linearization (J = sum coeff dR/dU)
    call this % prepare_inner_solver(coeff)

    allocate(res(system % get_num_state_vars()))

    ! the state owns the condensed update rule (one correction moves
    ! every derivative order, weighted by the linearization)
    s = differential_state(size(U, 1), size(U, 2) - 1, coeff)
    s % U = U

    r0 = 0.0_dp

    newton_iters: do iter = 1, this % max_it

       ! Evaluate the residual at the current state
       system % S = s % U

       res = 0.0d0
       call system % add_residual(res)

       rnorm = vector_norm(res)
       if (iter .eq. 1) r0 = rnorm

       if (rnorm .le. this % abs_tol .or. rnorm .le. this % rel_tol*r0) exit newton_iters

       ! Matrix-free linear solve  J dq = -res, delegated to the held solver
       ! (rhs unset => it forms -residual itself at the current state).
       call this % linear_solver % solve(system, dq, FORWARD)

       call s % update(dq)

    end do newton_iters

    U = s % U

    deallocate(res)

  end subroutine solve

  !===================================================================!
  ! Allocate (once) and configure the held inner linear solver to apply
  ! the current Newton/BDF linearization J = sum_n coeff(n) dR/dU(n).
  !===================================================================!

  impure subroutine prepare_inner_solver(this, coeff)

    class(newton), intent(inout) :: this
    type(scalar) , intent(in)    :: coeff(:)

    if (.not. allocated(this % linear_solver)) &
         & allocate(this % linear_solver, source = conjugate_gradient(1, 1.0d-14, 0))

    select type (ls => this % linear_solver)
    type is (conjugate_gradient)
       if (allocated(ls % lin_coeff)) deallocate(ls % lin_coeff)
       if (allocated(ls % rhs))       deallocate(ls % rhs)
       allocate(ls % lin_coeff(size(coeff)))
       ls % lin_coeff = real(coeff, dp)
    end select

  end subroutine prepare_inner_solver

  !===================================================================!
  ! Discrete adjoint total derivative of a functional:
  !   forward     R(u) = 0
  !   adjoint     (dR/du)^T psi = -df/du           (reverse-mode solve)
  !   derivative  dJ/dx = df/dx + psi^T dR/dx
  ! The adjoint state is handed back if requested.
  !===================================================================!

  impure subroutine eval_func_grad(this, system, func, dJdx, adjoint_state)

    class(newton)    , intent(inout)                     :: this
    class(assembler) , intent(inout)                     :: system
    class(functional), intent(in)                        :: func
    real(dp)         , intent(out), allocatable          :: dJdx(:)
    type(scalar)     , intent(out), allocatable, optional :: adjoint_state(:)

    type(scalar), allocatable :: dfdu(:)
    real(dp)    , allocatable :: psi(:)
    integer                   :: n, ndv

    n   = system % get_num_state_vars()
    ndv = system % get_num_design_vars()

    ! 1. forward steady solve R = 0 -> u (left in system % S): the
    ! inherited march at the steady linearization
    call march_steady(this, system)

    ! 2. adjoint right-hand side  -df/du
    allocate(dfdu(n))
    dfdu = 0.0d0
    call func % add_dfdu(system, dfdu)

    ! 3. adjoint solve  (dR/du)^T psi = -df/du, via the held solver in reverse
    call this % prepare_inner_solver(steady_coeff)
    select type (ls => this % linear_solver)
    type is (conjugate_gradient)
       allocate(ls % rhs(n))
       ls % rhs = real(-dfdu, dp)
    end select
    call this % linear_solver % solve(system, psi, REVERSE)

    ! 4. total derivative  dJ/dx = df/dx + psi^T dR/dx
    allocate(dJdx(ndv))
    dJdx = 0.0_dp
    call func % add_dfdx(system, dJdx)
    call system % add_design_residual_transpose_product(dJdx, psi)

    if (present(adjoint_state)) adjoint_state = psi

    deallocate(dfdu)

  end subroutine eval_func_grad

  !===================================================================!
  ! Verification gradient: central finite differences of J over each
  ! design variable, re-solving the forward problem at each perturbation.
  !===================================================================!

  impure subroutine eval_fd_func_grad(this, system, func, dJdx)

    class(newton)    , intent(inout)            :: this
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
       call march_steady(this, system)
       call func % eval(system, jp)

       x = x0; x(i) = x0(i) - delta
       call system % set_design_vars(x)
       call march_steady(this, system)
       call func % eval(system, jm)

       dJdx(i) = real(jp - jm, dp)/(2.0_dp*delta)

    end do

    ! restore the baseline design and state
    call system % set_design_vars(x0)
    call march_steady(this, system)

    deallocate(x0, x)

  end subroutine eval_fd_func_grad

  !===================================================================!
  ! Steady forward solve through the inherited march (fresh zero state
  ! at the steady linearization; the solution is left in system % S)
  !===================================================================!

  impure subroutine march_steady(this, system)

    class(newton)   , intent(inout) :: this
    class(assembler), intent(inout) :: system

    type(differential_state) :: s
    integer                  :: norder

    norder = system % get_differential_order()
    s = differential_state(system % get_num_state_vars(), norder, &
         & steady_coeff(1 : norder + 1))

    call this % march(system, s)

  end subroutine march_steady

  !===================================================================!
  ! 2-norm of a state vector (real or complex-step safe)
  !===================================================================!

  pure real(dp) function vector_norm(v)

    type(scalar), intent(in) :: v(:)

    vector_norm = sqrt(sum(abs(v)**2))

  end function vector_norm

end module class_newton_solver
