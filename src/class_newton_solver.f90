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

  implicit none

  private
  public :: newton

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

  end type newton

contains

  !===================================================================!
  ! Drive R(U) -> 0 for the newest state U = (nvars, order+1)
  !===================================================================!

  subroutine solve(this, system, coeff, U)

    class(newton)   , intent(in)    :: this
    class(assembler), intent(inout) :: system
    type(scalar)    , intent(in)    :: coeff(:)
    type(scalar)    , intent(inout) :: U(:,:)

    type(scalar), allocatable :: res(:), dq(:)
    real(dp)                  :: r0, rnorm
    integer                   :: nvars, n, iter

    nvars = system % get_num_state_vars()

    allocate(res(nvars))
    allocate(dq(nvars))

    r0 = 0.0_dp

    newton_iters: do iter = 1, this % max_it

       ! Evaluate the residual at the current state
       system % S = U

       res = 0.0d0
       call system % add_residual(res)

       rnorm = vector_norm(res)
       if (iter .eq. 1) r0 = rnorm

       if (rnorm .le. this % abs_tol .or. rnorm .le. this % rel_tol*r0) exit newton_iters

       ! Matrix-free linear solve  J dq = -res
       call cg_solve(system, coeff, -res, dq)

       ! Condensed update of every derivative order
       do n = 1, size(coeff)
          U(:,n) = U(:,n) + coeff(n)*dq
       end do

    end do newton_iters

    deallocate(res)
    deallocate(dq)

  end subroutine solve

  !===================================================================!
  ! Matrix-free conjugate gradient: solve J x = b where the operator
  ! J v = sum_n coeff(n) dR/dU(n) v is the assembler's jacobian-vector
  ! product (the state has already been set on the system).
  !===================================================================!

  subroutine cg_solve(system, coeff, b, x)

    class(assembler), intent(in)               :: system
    type(scalar)    , intent(in)               :: coeff(:)
    type(scalar)    , intent(in)               :: b(:)
    type(scalar)    , intent(out), allocatable :: x(:)

    type(scalar), allocatable :: r(:), p(:), Jp(:)
    type(scalar)              :: rs_old, rs_new, alpha, pJp
    integer                   :: nvars, it, max_it

    nvars  = system % get_num_state_vars()
    max_it = nvars + 100

    allocate(x(nvars))
    allocate(r(nvars))
    allocate(p(nvars))
    allocate(Jp(nvars))

    x = 0.0d0
    r = b
    p = r

    rs_old = dot_product(r, r)

    cg: do it = 1, max_it

       if (vector_norm(r) .le. 1.0d-14) exit cg

       Jp = 0.0d0
       call system % add_jacobian_vector_product(Jp, p, coeff)

       pJp = dot_product(p, Jp)
       if (abs(pJp) .le. tiny(1.0_dp)) exit cg

       alpha = rs_old/pJp

       x = x + alpha*p
       r = r - alpha*Jp

       rs_new = dot_product(r, r)

       p      = r + (rs_new/rs_old)*p
       rs_old = rs_new

    end do cg

    deallocate(r)
    deallocate(p)
    deallocate(Jp)

  end subroutine cg_solve

  !===================================================================!
  ! 2-norm of a state vector (real or complex-step safe)
  !===================================================================!

  pure real(dp) function vector_norm(v)

    type(scalar), intent(in) :: v(:)

    vector_norm = sqrt(sum(abs(v)**2))

  end function vector_norm

end module class_newton_solver
