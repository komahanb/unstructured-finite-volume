#include "scalar.fpp"

!=====================================================================!
! Transient discrete adjoint. Given a forward integrator that marches
! R_k(u_k, udot_k, x) = 0 and a function of interest J = sum_k w_k f(u_k),
! it sweeps the adjoint variables psi backward in time and forms the
! total derivative dJ/dx.
!
! At step m the adjoint equation is
!
!   J_m^T psi_m = -w_m df/du(u_m)
!                 - sum_{j>=1} (A(p_{m+j}, j+1)/h) M^T psi_{m+j},
!
! where J_m = (A(p_m,1)/h) M - A is the forward step jacobian and the
! sum couples to future steps whose BDF stencil still reaches back to m
! (j <= p_{m+j}). Both J_m^T and M^T act matrix-free through the
! assembler, so each step is a conjugate-gradient solve on the transpose
! operator - no assembled matrix. The total derivative is
!
!   dJ/dx = sum_m w_m df/dx + sum_m psi_m^T dR_m/dx.
!
! eval_func_grad marches forward, sweeps back and accumulates dJ/dx;
! eval_fd_func_grad re-marches at perturbed designs to verify it.
!
! First order in time (u, udot) - the order the fvm needs; the order-
! general sweep (udot, uddot) is a later extension.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_adjoint

  use iso_fortran_env     , only : dp => REAL64
  use interface_integrator, only : integrator
  use interface_assembler , only : assembler
  use interface_function  , only : functional

  implicit none

  private
  public :: adjoint

  !-------------------------------------------------------------------!
  ! Transient adjoint solver
  !-------------------------------------------------------------------!

  type :: adjoint

     class(integrator), allocatable :: fwd       ! forward integrator (+ trajectory)
     class(functional), allocatable :: func      ! function of interest
     type(scalar)     , allocatable :: psi(:,:)  ! adjoint history (step, nvars)
     real(dp)         , allocatable :: weight(:) ! per-step quadrature weight

   contains

     procedure :: eval_func_grad
     procedure :: eval_fd_func_grad
     procedure :: march_backwards
     procedure :: compute_total_derivative
     procedure :: integrate_functional
     procedure :: write_solution
     procedure :: write_gmsh_solution

  end type adjoint

  interface adjoint
     module procedure create
  end interface adjoint

contains

  !===================================================================!
  ! Construct from a forward integrator and a function of interest. The
  ! time integral uses J = h * sum_{k=2}^{N} f(u_k) (the initial state is
  ! given, so it carries no quadrature weight).
  !===================================================================!

  type(adjoint) function create(fwd, func) result(this)

    class(integrator), intent(in) :: fwd
    class(functional), intent(in) :: func

    integer :: n

    allocate(this % fwd , source = fwd)
    allocate(this % func, source = func)

    n = this % fwd % num_steps
    allocate(this % weight(n))
    this % weight    = this % fwd % h
    this % weight(1) = 0.0_dp

  end function create

  !===================================================================!
  ! Adjoint total derivative: march forward, sweep backward, accumulate
  !===================================================================!

  subroutine eval_func_grad(this, dJdx)

    class(adjoint), intent(inout)            :: this
    real(dp)      , intent(out), allocatable :: dJdx(:)

    call this % fwd % solve()
    call this % march_backwards()
    call this % compute_total_derivative(dJdx)

  end subroutine eval_func_grad

  !===================================================================!
  ! Backward sweep filling psi(N..2). psi(1) (the given initial state)
  ! stays zero - it has no residual equation.
  !===================================================================!

  subroutine march_backwards(this)

    class(adjoint), intent(inout) :: this

    type(scalar), allocatable :: rhs(:), dfdu(:), Mpsi(:), psi_m(:)
    type(scalar), allocatable :: scoeff(:), lin_coeff(:)
    type(scalar)              :: mass_coeff(2)
    integer                   :: n, m, k, j, p_m, p_k, kmax, nvars
    real(dp)                  :: h

    n = this % fwd % num_steps
    h = this % fwd % h

    ! mass action selector  [alpha, beta] = [0, 1]  ->  M v
    mass_coeff = [0.0_dp, 1.0_dp]

    associate(system => this % fwd % system)

      nvars = system % get_num_state_vars()

      if (allocated(this % psi)) deallocate(this % psi)
      allocate(this % psi(n, nvars)); this % psi = 0.0d0

      allocate(rhs(nvars), dfdu(nvars), Mpsi(nvars), lin_coeff(2))

      backward: do m = n, 2, -1

         ! forward step jacobian coefficients at m: J_m = beta M - alpha A
         ! with [alpha, beta] = [1, A(p_m,1)/h]
         p_m = this % fwd % get_bandwidth(m)
         call this % fwd % get_stencil_coeff(p_m, h, scoeff)
         lin_coeff = [1.0d0, scoeff(1)]

         ! position the assembler at step m (state u_m, udot_m)
         system % S = this % fwd % U(m,:,:)

         ! functional contribution  -w_m df/du(u_m)
         dfdu = 0.0d0
         call this % func % add_dfdu(system, dfdu)
         rhs = -this % weight(m)*dfdu

         ! future-step coupling: step k = m+j contributes (A(p_k,j+1)/h) M psi_k
         ! whenever its stencil still reaches back to m (j <= p_k). Bandwidth
         ! is monotone, so it suffices to scan to m + bandwidth(N).
         kmax = min(n, m + this % fwd % get_bandwidth(n))
         couple: do k = m+1, kmax
            j   = k - m
            p_k = this % fwd % get_bandwidth(k)
            if (j .gt. p_k) cycle couple
            call this % fwd % get_stencil_coeff(p_k, h, scoeff)
            Mpsi = 0.0d0
            call system % add_jacobian_vector_product_transpose( &
                 & Mpsi, this % psi(k,:), mass_coeff)
            rhs = rhs - scoeff(j+1)*Mpsi
         end do couple

         ! solve  J_m^T psi_m = rhs
         call cg_solve_transpose(system, lin_coeff, rhs, psi_m)
         this % psi(m,:) = psi_m
         deallocate(psi_m)

      end do backward

      deallocate(rhs, dfdu, Mpsi, lin_coeff)
      if (allocated(scoeff)) deallocate(scoeff)

    end associate

  end subroutine march_backwards

  !===================================================================!
  ! Total derivative  dJ/dx = sum_m w_m df/dx + sum_m psi_m^T dR_m/dx
  !===================================================================!

  subroutine compute_total_derivative(this, dJdx)

    class(adjoint), intent(inout)            :: this
    real(dp)      , intent(out), allocatable :: dJdx(:)

    real(dp), allocatable :: dfdx_m(:)
    integer               :: n, m, ndv

    n = this % fwd % num_steps

    associate(system => this % fwd % system)

      ndv = system % get_num_design_vars()
      allocate(dJdx(ndv));   dJdx   = 0.0_dp
      allocate(dfdx_m(ndv))

      do m = 2, n

         ! the design partial uses the step's full state (u_m, udot_m)
         system % S = this % fwd % U(m,:,:)

         ! explicit functional design dependence (weighted; default zero)
         dfdx_m = 0.0_dp
         call this % func % add_dfdx(system, dfdx_m)
         dJdx = dJdx + this % weight(m)*dfdx_m

         ! adjoint residual design contribution  psi_m^T dR_m/dx
         call system % add_design_residual_transpose_product(dJdx, this % psi(m,:))

      end do

      deallocate(dfdx_m)

    end associate

  end subroutine compute_total_derivative

  !===================================================================!
  ! Verification gradient: central differences of the time-integrated
  ! functional, re-marching the forward problem at each perturbation
  !===================================================================!

  subroutine eval_fd_func_grad(this, dJdx)

    class(adjoint), intent(inout)            :: this
    real(dp)      , intent(out), allocatable :: dJdx(:)

    real(dp), allocatable :: x0(:), x(:)
    real(dp)              :: jp, jm, delta
    integer               :: i, ndv

    associate(system => this % fwd % system)

      ndv = system % get_num_design_vars()
      allocate(dJdx(ndv), x0(ndv), x(ndv))

      call system % get_design_vars(x0)

      do i = 1, ndv

         delta = 1.0e-6_dp*max(1.0_dp, abs(x0(i)))

         x = x0; x(i) = x0(i) + delta
         call system % set_design_vars(x)
         call this % fwd % solve()
         jp = this % integrate_functional()

         x = x0; x(i) = x0(i) - delta
         call system % set_design_vars(x)
         call this % fwd % solve()
         jm = this % integrate_functional()

         dJdx(i) = (jp - jm)/(2.0_dp*delta)

      end do

      ! restore the baseline design and trajectory
      call system % set_design_vars(x0)
      call this % fwd % solve()

      deallocate(x0, x)

    end associate

  end subroutine eval_fd_func_grad

  !===================================================================!
  ! Export the state and adjoint-state trajectories for post-processing:
  ! one file per step, "<basename>_NNNN.vtu", each carrying the cell
  ! fields "state" (u_k) and "adjoint" (psi_k). Numbered so paraview
  ! reads the sequence as time steps. (psi_1 = 0 - the initial state has
  ! no adjoint equation.)
  !===================================================================!

  subroutine write_solution(this, basename)

    class(adjoint), intent(inout) :: this
    character(len=*), intent(in)  :: basename

    real(dp), allocatable :: fields(:,:)
    character(len=256)    :: fname
    integer               :: m, n, nvars

    n = this % fwd % num_steps

    associate(system => this % fwd % system)

      nvars = system % get_num_state_vars()
      allocate(fields(nvars, 2))

      do m = 1, n
         fields(:,1) = real(this % fwd % U(m,:,1), dp)   ! state u_m
         fields(:,2) = real(this % psi(m,:)      , dp)   ! adjoint psi_m
         write(fname, '(a,a,i0.4,a)') trim(basename), "_", m-1, ".vtu"
         call system % write_solution_fields( &
              & trim(fname), fields, [character(len=7) :: "state", "adjoint"])
      end do

      deallocate(fields)

    end associate

  end subroutine write_solution

  !===================================================================!
  ! Export the state and adjoint-state trajectories to a single gmsh
  ! post-processing file: two animated views ("state", "adjoint") over
  ! the time steps. meshfile is the source mesh copied verbatim. (psi_1
  ! = 0 - the initial state has no adjoint equation.)
  !===================================================================!

  subroutine write_gmsh_solution(this, meshfile, filename)

    class(adjoint), intent(inout) :: this
    character(len=*), intent(in)  :: meshfile, filename

    real(dp), allocatable :: fields(:,:,:)   ! (ndof, 2, nstep)
    real(dp), allocatable :: times(:)
    integer               :: m, n, nvars

    n = this % fwd % num_steps

    associate(system => this % fwd % system)

      nvars = system % get_num_state_vars()
      allocate(fields(nvars, 2, n), times(n))

      do m = 1, n
         fields(:,1,m) = real(this % fwd % U(m,:,1), dp)   ! state u_m
         fields(:,2,m) = real(this % psi(m,:)      , dp)   ! adjoint psi_m
         times(m)      = this % fwd % time(m)
      end do

      call system % write_gmsh_series(meshfile, filename, fields, &
           & [character(len=7) :: "state", "adjoint"], times)

      deallocate(fields, times)

    end associate

  end subroutine write_gmsh_solution

  !===================================================================!
  ! Time-integrated functional  J = sum_k w_k f(u_k) over the trajectory
  !===================================================================!

  real(dp) function integrate_functional(this) result(jval)

    class(adjoint), intent(inout) :: this

    type(scalar) :: fm
    integer      :: n, m

    n    = this % fwd % num_steps
    jval = 0.0_dp

    associate(system => this % fwd % system)
      do m = 1, n
         system % S = this % fwd % U(m,:,:)
         call this % func % eval(system, fm)
         jval = jval + this % weight(m)*real(fm, dp)
      end do
    end associate

  end function integrate_functional

  !===================================================================!
  ! Matrix-free conjugate gradient on the TRANSPOSE jacobian operator
  !   B v = [scalar(i) dR/dU(i)]^T v
  ! supplied by add_jacobian_vector_product_transpose. For diffusion the
  ! step jacobian beta M - alpha A is symmetric positive definite.
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

end module class_adjoint
