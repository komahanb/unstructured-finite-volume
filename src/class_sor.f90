!=====================================================================!
! Relaxed-triangle-sweep linear solver (traditionally: successive
! over-relaxation, SOR): supplies only the sweep (`iterate`) and
! inherits the residual-minimization iteration from linear_solver.
! omega is the pseudo-time step size of the outer iteration.
!
! Two sweeps, chosen by what the solver carries:
!
!    no graph:  the relaxed lower-triangle solve (D+wL)y = R is
!               itself done iteratively (products through the
!               triangle filters)
!
!    carrying a graph:  the inherited colored sweep - exact, color
!               by color, relaxed by omega. the graph may be the
!               system's mesh or the assembled operator itself
!               (mesh-free smoothing; see gauss_seidel).
!=====================================================================!

module class_sor

  use iso_fortran_env         , only : dp => REAL64
  use interface_graph         , only : graph, power_iteration
  use interface_linear_solver , only : linear_solver
  use interface_assembler     , only : assembler, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE

  implicit none

  ! Expose only the linear solver datatype
  private
  public :: sor

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!

  type, extends(linear_solver) :: sor

     real(dp) :: omega

     ! parameter-selection state: the previous parameter value and the
     ! rate it achieved, for reverting a change that worsens the rate
     real(dp) :: omega_saved = -1.0_dp
     real(dp) :: rate_saved  = huge(1.0_dp)

   contains

     ! the sweep consumed by the inherited outer iteration
     procedure :: iterate

     ! auto-tuning: static optimal omega from the measured convergence
     ! factor, dynamic rollback if a parameter change worsens the rate
     procedure :: tune

  end type sor

  !===================================================================!
  ! Interface for multiple constructors
  !===================================================================!

  interface sor
     module procedure construct
  end interface sor

contains

  !===================================================================!
  ! Constructor for linear solver
  !===================================================================!

  pure type(sor) function construct(omega, max_it, &
       & max_tol, print_level, g) result (this)

    type(real(dp)) , intent(in)           :: omega
    type(integer)  , intent(in)           :: max_it
    type(real(dp)) , intent(in)           :: max_tol
    type(integer)  , intent(in)           :: print_level
    class(graph)   , intent(in), optional :: g

    this % max_it      = max_it
    this % max_tol     = max_tol
    this % print_level = print_level
    this % omega       = omega
    this % res_file    = 'sor.res'

    if (present(g)) call this % carry(g)

  end function construct

  !===================================================================!
  ! Parameter selection. At entry (pass 0): measure the convergence
  ! factor of the unrelaxed sweep - the spectral radius of D^-1 (L+U),
  ! by power iteration composed from the operator parts - and set the
  ! classical optimum omega = 2 / (1 + sqrt(1 - rho^2)).
  ! During the iteration (pass > 1): if the rate measured after a
  ! parameter change is worse than before, revert to the previous
  ! value.
  !===================================================================!

  impure subroutine tune(this, system, pass, rate)

    class(sor)      , intent(inout) :: this
    class(assembler), intent(in)    :: system
    integer         , intent(in)    :: pass
    real(dp)        , intent(in)    :: rate

    real(dp), allocatable :: D(:), identity(:)
    real(dp) :: rho
    integer  :: n

    static: if (pass .eq. 0) then

       n = system % num_state_vars
       allocate(D(n), identity(n))

       ! the diagonal (the self-loop subgraph on ones)
       identity = 1.0_dp
       call system % get_jacobian_residual_product(D, identity, part = DIAGONAL)

       ! the growth rate of the unrelaxed sweep map v -> D^-1 (L+U) v,
       ! measured by the graph module's power kernel
       call power_iteration(sweep_map, n, 50, rho)
       rho = min(rho, 1.0_dp - epsilon(1.0_dp))

       this % omega_saved = this % omega
       this % rate_saved  = huge(1.0_dp)
       this % omega       = 2.0_dp/(1.0_dp + sqrt(max(0.0_dp, 1.0_dp - rho*rho)))

       if (this % print_level .gt. 0) then
          write(*,'(1x,a,es12.5,a,f8.5)') &
               & "sor tune: measured rho = ", rho, ", omega set to ", this % omega
       end if

    else static

       ! revert a change that worsened the measured rate
       if (rate .gt. this % rate_saved .and. &
            & abs(this % omega - this % omega_saved) .gt. tiny(1.0_dp)) then
          this % omega = this % omega_saved
          if (this % print_level .gt. 0) then
             write(*,'(1x,a,f8.5)') "sor tune: rate worsened, reverting omega to ", this % omega
          end if
       else
          this % rate_saved  = rate
          this % omega_saved = this % omega
       end if

    end if static

  contains

    ! the unrelaxed sweep as a linear action (D and system captured)
    subroutine sweep_map(v, w)
      real(dp), intent(in)  :: v(:)
      real(dp), intent(out) :: w(:)
      real(dp), allocatable :: wu(:)
      allocate(wu(size(w)))
      call system % get_jacobian_residual_product(w , v, part = LOWER_TRIANGLE)
      call system % get_jacobian_residual_product(wu, v, part = UPPER_TRIANGLE)
      w = (w + wu)/D
    end subroutine sweep_map

  end subroutine tune

  !===================================================================!
  ! The sweep: relaxed (omega) triangle sweep on the correction
  ! equation A dx = r from dx = 0. The relaxed lower-triangle solve
  ! (D+wL)y = R is itself iterative.
  !===================================================================!

  impure subroutine iterate(this, system, r, dx, iter)

    class(sor)          , intent(inout) :: this
    class(assembler)    , intent(in)  :: system
    real(dp)            , intent(in)  :: r(:)
    real(dp)            , intent(out) :: dx(:)
    integer             , intent(out) :: iter

    ! Locals
    real(dp) :: tol, bnorm
    real(dp) , allocatable :: Ux(:), D(:), Ddx(:), R2(:), xnew(:), identity(:)
    real(dp) , allocatable :: Ly(:), y(:), ynew(:)

    ! carrying the system's graph buys the exact relaxed sweep
    if (allocated(this % g)) then
       call this % colored_sweep(system, r, dx, iter, this % omega)
       return
    end if

    dx = 0.0_dp
    allocate(Ux, D, Ddx, R2, xnew, identity, mold = dx)
    allocate(y, ynew, Ly, mold = dx)

    ! Extract the diagonal entries (the self-loop subgraph on ones)
    identity = 1.0d0
    call system % get_jacobian_residual_product(D, identity, part = DIAGONAL)

    bnorm = sqrt(system % inner_product(r, r))

    ! Homogeneous case (nothing to do)
    if (bnorm .le. this % max_tol) then
       iter = 0
       return
    end if

    !-----------------------------------------------------------------!
    ! Apply the relaxed sweep until tolerance is achieved
    !-----------------------------------------------------------------!

    iter = 1; tol  = huge(1.0d0)
    do while ((tol .gt. this % max_tol) .and. (iter .lt. this % max_it))

       ! Form Ux and Dx
       call system % get_jacobian_residual_product(Ux, dx, part = UPPER_TRIANGLE)
       call system % get_jacobian_residual_product(Ddx, dx, part = DIAGONAL)

       ! R = w(r-Ux_k)+(1-w)Dx_k
       R2 = this % omega * (r - Ux) + (1.0_dp-this % omega)*Ddx

       !--------------------------------------------------------------!
       ! Solve the linear system: By=R ; (D+wL)y=R ;  Dy=R-wLy
       !--------------------------------------------------------------!

       solve_lower_triangle: block

         real(dp) :: tol2
         integer :: iter2

         ! Initial guess is the current correction
         y = dx

         iter2 = 1; tol2 = huge(1.0d0)
         do while ((tol2 .gt. this % max_tol) .and. (iter2 .lt. this % max_it))

            call system % get_jacobian_residual_product(Ly, y, part = LOWER_TRIANGLE)

            ynew  = (R2 - this % omega*Ly)/D
            tol2  = sqrt(system % inner_product(y - ynew, y - ynew))

            if (this % print_level .gt. 2) then
               write(*,*) "inner (2)", iter2, tol2
            end if

            y = ynew
            iter2 = iter2 + 1

         end do

       end block solve_lower_triangle

       xnew = y
       tol  = sqrt(system % inner_product(dx - xnew, dx - xnew))

       if (this % print_level .gt. 1) then
          write(*,*) "inner (1)", iter, tol
       end if

       dx   = xnew
       iter = iter + 1

    end do

    deallocate(Ux, D, Ddx, R2, xnew, identity)
    deallocate(y, ynew, Ly)

  end subroutine iterate

end module class_sor
