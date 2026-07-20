#include "scalar.fpp"

!=====================================================================!
! Backward Difference Formula time integrator (orders 1-6).
!
! The time steps form a chain; a stencil that reads p states back
! raises that chain to its p-th power - and this integrator carries
! the resulting step dag (a chain object, edges by rule):
!
!                .-----------.-----------.
!                |           v           v
!                1 --> 2 --> 3 --> 4 --> 5        edge m --> k
!                      |           ^              whenever
!                      '-----------'              k - m <= power
!
! One dag, two directions: forward, the in-edges of a vertex deliver
! the past states that form udot (step); backward, the out-edges
! hand back weighted mass actions (march_backwards). The same edge
! weights A(p,j+1)/h ride both ways.
!
! Now zoom into any vertex, and the picture repeats at a smaller
! scale - this solver is a fractal of graphs:
!
!    the step dag           1 --> 2 --> (m) --> ... --> n
!                                        |
!    inside vertex m:                    v
!    one transpose solve       it1 --> it2 --> ... --> itj
!    = a solver chain of                 |
!    matvecs and dots                    v
!                                    o---o---o      a matvec: at
!    inside one matvec:              | \ | / |      every mesh
!    the mesh graph                  o---o---o      vertex, a dot
!                                        |          over its edges
!    inside one dot:                     v
!    the partition           (o o o)  +  (o o o)  +  (o o)
!                             part 1      part 2      part 3
!
! Four rungs, one grammar: weighted sums over edges and sums over
! vertex sets - inner products all the way down.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_bdf

  use iso_fortran_env         , only : dp => REAL64
  use interface_integrator    , only : integrator
  use interface_assembler     , only : assembler
  use interface_function      , only : functional
  use class_chain             , only : chain
  use class_newton_solver     , only : newton
  use class_conjugate_gradient, only : conjugate_gradient
  use module_solve_mode       , only : REVERSE

  implicit none

  private
  public :: bdf

  !-------------------------------------------------------------------!
  ! BDF integrator type
  !-------------------------------------------------------------------!

  type, extends(integrator) :: bdf

     integer                   :: max_order = 6
     type(scalar), allocatable :: A(:,:)        ! BDF coefficient table
     type(chain)               :: steps         ! the step dag: chain(n) at power max_order
     ! the adjoint trajectory psi(step, nvars) lives on the integrator base
     ! alongside the primal trajectory U

   contains

     procedure :: step
     procedure :: get_bandwidth
     procedure :: get_linearization_coeff
     procedure :: get_stencil_coeff

     ! transient discrete adjoint (backward sweep over the trajectory)
     procedure :: integrate_adjoint
     procedure :: integrate_adjoint_fd
     procedure :: write_adjoint_solution
     procedure :: write_adjoint_gmsh
     procedure, private :: march_backwards
     procedure, private :: compute_total_derivative
     procedure, private :: integrate_functional

  end type bdf

  interface bdf
     module procedure create
  end interface bdf

contains

  !===================================================================!
  ! Construct a BDF integrator of the requested accuracy order
  !===================================================================!

  pure type(bdf) function create(system, tinit, tfinal, h, max_order) result(this)

    class(assembler), intent(in) :: system
    real(dp)        , intent(in) :: tinit, tfinal, h
    integer         , intent(in) :: max_order

    call this % construct(system, tinit, tfinal, h, implicit = .true.)

    this % max_order = min(max_order, 6)

    ! Coefficient table - http://www.scholarpedia.org/article/Backward_differentiation_formulas
    allocate(this % A(this % max_order, this % max_order + 1))
    this % A = 0.0d0

    if (this % max_order .ge. 1) this % A(1,1:2) = [1.0d0, -1.0d0]
    if (this % max_order .ge. 2) this % A(2,1:3) = [3.0d0, -4.0d0, 1.0d0]/2.0d0
    if (this % max_order .ge. 3) this % A(3,1:4) = [11.0d0, -18.0d0, 9.0d0, -2.0d0]/6.0d0
    if (this % max_order .ge. 4) this % A(4,1:5) = [25.0d0, -48.0d0, 36.0d0, -16.0d0, 3.0d0]/12.0d0
    if (this % max_order .ge. 5) this % A(5,1:6) = [137.0d0, -300.0d0, 300.0d0, -200.0d0, 75.0d0, -12.0d0]/60.0d0
    if (this % max_order .ge. 6) this % A(6,1:7) = [147.0d0, -360.0d0, 450.0d0, -400.0d0, 225.0d0, -72.0d0, 10.0d0]/60.0d0

  end function create

  !===================================================================!
  ! Number of past steps the formula uses (ramps up to max_order)
  !===================================================================!

  pure integer function get_bandwidth(this, step_index) result(width)

    class(bdf), intent(in) :: this
    integer   , intent(in) :: step_index

    width = step_index - 1

    if (width .gt. this % max_order) width = this % max_order

  end function get_bandwidth

  !===================================================================!
  ! Linearization coefficients  coeff(n+1) = (A(p,1)/h)^n,  n = 0..order
  !===================================================================!

  pure subroutine get_linearization_coeff(this, p, h, coeff)

    class(bdf)  , intent(in)    :: this
    integer     , intent(in)    :: p
    real(dp)    , intent(in)    :: h
    type(scalar), intent(inout) :: coeff(:)

    integer :: n

    do n = 0, this % system % get_differential_order()
       coeff(n+1) = (this % A(p,1)/h)**n
    end do

  end subroutine get_linearization_coeff

  !===================================================================!
  ! First-derivative stencil at bandwidth p:  scoeff(i+1) = A(p,i+1)/h,
  ! i = 0..p - the same coefficients that form udot from past states in
  ! step(). The adjoint backward sweep uses the whole stencil.
  !===================================================================!

  pure subroutine get_stencil_coeff(this, p, h, scoeff)

    class(bdf)  , intent(in)               :: this
    integer     , intent(in)               :: p
    real(dp)    , intent(in)               :: h
    type(scalar), intent(out), allocatable :: scoeff(:)

    integer :: i

    allocate(scoeff(p+1))

    do i = 0, p
       scoeff(i+1) = this % A(p, i+1)/h
    end do

  end subroutine get_stencil_coeff

  !===================================================================!
  ! Advance one step: the in-edges of the new vertex deliver the past
  ! states, and udot is their weighted sum -
  !
  !    k-p  ...  k-2   k-1                                  A(p,j+1)
  !      \        |     |         udot_k = sum  w_j u_{k-j},  w_j = --------
  !       \       v     v              in-edges                        h
  !        '----> ( k )
  !
  ! then (implicit) newton drives R(u_k, udot_k) to zero at the vertex.
  !===================================================================!

  impure subroutine step(this, t, U, h, p, ierr)

    class(bdf)  , intent(inout) :: this
    real(dp)    , intent(inout) :: t(:)
    type(scalar), intent(inout) :: U(:,:,:)       ! (window, nvars, order+1)
    integer     , intent(in)    :: p
    real(dp)    , intent(in)    :: h
    integer     , intent(out)   :: ierr

    type(scalar), allocatable :: coeff(:)
    type(newton)              :: nlsolver
    integer                   :: kk, torder, n, i

    ierr   = 0
    kk     = size(U, dim=1)                       ! newest state at window end
    torder = this % system % get_differential_order()

    ! Advance the time
    t(kk) = t(kk-1) + h

    ! Predictor: carry the last value for the lowest order
    U(kk,:,1) = U(kk-1,:,1)

    ! Higher-order states from the BDF stencil
    do n = 1, torder

       U(kk,:,n+1) = 0.0d0

       do i = 0, p
          U(kk,:,n+1) = U(kk,:,n+1) + (this % A(p,i+1)/h)*U(kk-i,:,n)
       end do

    end do

    ! Implicit correction driving the residual to zero
    if (this % implicit) then

       allocate(coeff(torder+1))

       call this % get_linearization_coeff(p, h, coeff)
       call nlsolver % solve(this % system, coeff, U(kk,:,:))

       deallocate(coeff)

    end if

  end subroutine step

  !===================================================================!
  ! Transient discrete adjoint, J = h sum_{k>=2} f(u_k):
  !
  !    march forward  ============>  u_1 ... u_n     (integrate)
  !    sweep backward <============  psi_n ... psi_2 (march_backwards)
  !    gather dJ/dx over every vertex visited        (total derivative)
  !
  ! The step jacobian J_m^T and the mass M^T act matrix-free through
  ! the assembler, so each backward vertex is one transpose CG solve.
  !===================================================================!

  impure subroutine integrate_adjoint(this, func, dJdx)

    class(bdf)       , intent(inout)            :: this
    class(functional), intent(in)               :: func
    real(dp)         , intent(out), allocatable :: dJdx(:)

    call this % integrate()                         ! forward march -> this % U
    call this % march_backwards(func)           ! fill this % psi (N..2)
    call this % compute_total_derivative(func, dJdx)

  end subroutine integrate_adjoint

  !===================================================================!
  ! Backward sweep: visit the step dag in reverse dependency order,
  ! and at each vertex pull one weighted sum in over the out-edges -
  ! the same inner-product primitive the linear solvers run on:
  !
  !               ( m ) ----> k = m+1 ... m+p      each out-edge hands
  !                 ^          |                   back its mass action
  !                 |          v                   M^T psi_k, scaled by
  !                 |     A(p_k, k-m+1)            the edge weight
  !                 |    --------------- M^T psi_k
  !                 '---------- h
  !
  !    rhs_m = -w_m df/du(u_m) - sum ( weight * M^T psi_k )
  !                         out-edges
  !
  ! then one transpose solve at the vertex: J_m^T psi_m = rhs_m.
  ! psi(1) stays zero - the initial state has no residual equation.
  !===================================================================!

  impure subroutine march_backwards(this, func)

    class(bdf)       , intent(inout) :: this
    class(functional), intent(in)    :: func

    type(scalar), allocatable :: rhs(:), dfdu(:), Mpsi(:), scoeff(:)
    real(dp)    , allocatable :: psi_m(:)
    integer     , allocatable :: order(:), nbrs(:)
    type(scalar)              :: mass_coeff(2), lin_coeff(2)
    type(conjugate_gradient)  :: lsolver
    integer                   :: n, m, k, j, i, e, p_m, p_k, nvars
    real(dp)                  :: h, w

    n = this % num_steps
    h = this % h

    ! the step dag: the chain of n steps raised to the stencil depth.
    ! its edges ARE the couplings - edge m -> k exactly when step k's
    ! stencil reads step m - so no reach/cutoff bookkeeping survives
    ! here; the sweep just walks the structure.
    this % steps = chain(n, power = this % max_order)
    order = this % steps % dependency_order()

    ! mass action selector [alpha, beta] = [0, 1] -> M v
    mass_coeff = [0.0_dp, 1.0_dp]

    associate(system => this % system)

      nvars = system % get_num_state_vars()

      if (allocated(this % psi)) deallocate(this % psi)
      allocate(this % psi(n, nvars)); this % psi = 0.0d0

      allocate(rhs(nvars), dfdu(nvars), Mpsi(nvars))

      backward: do i = n, 1, -1

         m = order(i)
         if (m .eq. 1) cycle backward   ! the initial state: no residual eqn

         ! forward step jacobian coefficients at m: [alpha, beta] = [1, A(p_m,1)/h]
         p_m = this % get_bandwidth(m)
         call this % get_stencil_coeff(p_m, h, scoeff)
         lin_coeff = [1.0d0, scoeff(1)]

         ! position the assembler at step m (state u_m, udot_m)
         system % S = this % U(m,:,:)

         ! functional contribution -w_m df/du(u_m)
         w    = merge(h, 0.0_dp, m >= 2)
         dfdu = 0.0d0
         call func % add_dfdu(system, dfdu)
         rhs  = -w*dfdu

         ! every out-edge hands back its weighted mass action
         nbrs = this % steps % out_neighbours(m)
         couple: do e = 1, size(nbrs)
            k   = nbrs(e)
            j   = k - m
            p_k = this % get_bandwidth(k)
            call this % get_stencil_coeff(p_k, h, scoeff)
            Mpsi = 0.0d0
            call system % add_jacobian_vector_product_transpose(Mpsi, this % psi(k,:), mass_coeff)
            rhs = rhs - scoeff(j+1)*Mpsi
         end do couple

         ! solve J_m^T psi_m = rhs via the linear solver in reverse (transpose)
         if (allocated(lsolver % lin_coeff)) deallocate(lsolver % lin_coeff)
         if (allocated(lsolver % rhs))       deallocate(lsolver % rhs)
         allocate(lsolver % lin_coeff(2));      lsolver % lin_coeff = real(lin_coeff, dp)
         allocate(lsolver % rhs(nvars));        lsolver % rhs       = real(rhs, dp)
         call lsolver % solve(system, psi_m, REVERSE)
         this % psi(m,:) = psi_m
         deallocate(psi_m)

      end do backward

      deallocate(rhs, dfdu, Mpsi)
      if (allocated(scoeff)) deallocate(scoeff)

    end associate

  end subroutine march_backwards

  !===================================================================!
  ! Total derivative  dJ/dx = sum_m w_m df/dx + sum_m psi_m^T dR_m/dx
  !===================================================================!

  impure subroutine compute_total_derivative(this, func, dJdx)

    class(bdf)       , intent(inout)            :: this
    class(functional), intent(in)               :: func
    real(dp)         , intent(out), allocatable :: dJdx(:)

    real(dp), allocatable :: dfdx_m(:)
    integer               :: n, m, ndv
    real(dp)              :: w

    n = this % num_steps

    associate(system => this % system)

      ndv = system % get_num_design_vars()
      allocate(dJdx(ndv)); dJdx = 0.0_dp
      allocate(dfdx_m(ndv))

      do m = 2, n
         system % S = this % U(m,:,:)
         w = merge(this % h, 0.0_dp, m >= 2)

         dfdx_m = 0.0_dp
         call func % add_dfdx(system, dfdx_m)
         dJdx = dJdx + w*dfdx_m

         call system % add_design_residual_transpose_product(dJdx, this % psi(m,:))
      end do

      deallocate(dfdx_m)

    end associate

  end subroutine compute_total_derivative

  !===================================================================!
  ! Time-integrated functional  J = sum_k w_k f(u_k) over the trajectory
  !===================================================================!

  impure real(dp) function integrate_functional(this, func) result(jval)

    class(bdf)       , intent(inout) :: this
    class(functional), intent(in)    :: func

    type(scalar) :: fm
    integer      :: n, m
    real(dp)     :: w

    n    = this % num_steps
    jval = 0.0_dp

    associate(system => this % system)
      do m = 1, n
         system % S = this % U(m,:,:)
         w = merge(this % h, 0.0_dp, m >= 2)
         call func % eval(system, fm)
         jval = jval + w*real(fm, dp)
      end do
    end associate

  end function integrate_functional

  !===================================================================!
  ! Verification gradient: central differences of the time-integrated
  ! functional, re-marching the forward problem at each perturbation.
  !===================================================================!

  impure subroutine integrate_adjoint_fd(this, func, dJdx)

    class(bdf)       , intent(inout)            :: this
    class(functional), intent(in)               :: func
    real(dp)         , intent(out), allocatable :: dJdx(:)

    real(dp), allocatable :: x0(:), x(:)
    real(dp)              :: jp, jm, delta
    integer               :: i, ndv

    associate(system => this % system)

      ndv = system % get_num_design_vars()
      allocate(dJdx(ndv), x0(ndv), x(ndv))

      call system % get_design_vars(x0)

      do i = 1, ndv

         delta = 1.0e-6_dp*max(1.0_dp, abs(x0(i)))

         x = x0; x(i) = x0(i) + delta
         call system % set_design_vars(x)
         call this % integrate()
         jp = this % integrate_functional(func)

         x = x0; x(i) = x0(i) - delta
         call system % set_design_vars(x)
         call this % integrate()
         jm = this % integrate_functional(func)

         dJdx(i) = (jp - jm)/(2.0_dp*delta)

      end do

      ! restore the baseline design and trajectory
      call system % set_design_vars(x0)
      call this % integrate()

      deallocate(x0, x)

    end associate

  end subroutine integrate_adjoint_fd

  !===================================================================!
  ! Export the state + adjoint-state trajectories: one .vtu per step,
  ! "<basename>_NNNN.vtu", each with cell fields "state" (u_k) and
  ! "adjoint" (psi_k). (psi_1 = 0 - the initial state has no adjoint eqn.)
  !===================================================================!

  impure subroutine write_adjoint_solution(this, basename)

    class(bdf)      , intent(inout) :: this
    character(len=*), intent(in)    :: basename

    real(dp), allocatable :: fields(:,:)
    character(len=256)    :: fname
    integer               :: m, n, nvars

    n = this % num_steps

    associate(system => this % system)
      nvars = system % get_num_state_vars()
      allocate(fields(nvars, 2))
      do m = 1, n
         fields(:,1) = real(this % U(m,:,1), dp)
         fields(:,2) = real(this % psi(m,:), dp)
         write(fname, '(a,a,i0.4,a)') trim(basename), "_", m-1, ".vtu"
         call system % write_solution_fields(trim(fname), fields, &
              & [character(len=7) :: "state", "adjoint"])
      end do
      deallocate(fields)
    end associate

  end subroutine write_adjoint_solution

  !===================================================================!
  ! Export the state + adjoint-state trajectories to one gmsh post file
  ! (two animated views over the time steps).
  !===================================================================!

  impure subroutine write_adjoint_gmsh(this, meshfile, filename)

    class(bdf)      , intent(inout) :: this
    character(len=*), intent(in)    :: meshfile, filename

    real(dp), allocatable :: fields(:,:,:), times(:)
    integer               :: m, n, nvars

    n = this % num_steps

    associate(system => this % system)
      nvars = system % get_num_state_vars()
      allocate(fields(nvars, 2, n), times(n))
      do m = 1, n
         fields(:,1,m) = real(this % U(m,:,1), dp)
         fields(:,2,m) = real(this % psi(m,:), dp)
         times(m)      = this % time(m)
      end do
      call system % write_gmsh_series(meshfile, filename, fields, &
           & [character(len=7) :: "state", "adjoint"], times)
      deallocate(fields, times)
    end associate

  end subroutine write_adjoint_gmsh

end module class_bdf
