#include "scalar.fpp"

!=====================================================================!
! The generic linear solver, plus the preconditioner contract it
! consumes. A linear solver takes an initial guess and improves it
! until the system's residual is eliminated. It uses matrix, vector
! and inner product only - it asks the system for the residual and the
! product, and never names a source or a correction scheme.
!
! The provided march (implemented by converge) is the residual-minimization
! iteration - one cycle of state -> residual -> correction -> state per
! pass, terminating on the true residual. The deferred iterate is
! the one thing a concrete solver must supply: one application of the
! approximate inverse, driving A dx = r from dx = 0. This mirrors
! interface_integrator, where the provided integrate marches the
! deferred step.
!
! Two optional operator slots accelerate the kernels: pre_conditioner
! (acts on residuals - the left preconditioner / smoother) and
! post_conditioner (acts on answers - the right preconditioner). The
! slots are storage, not outer-loop actions: each kernel applies its
! slot where its own algorithm requires (cg inside every iteration,
! gmres inside its arnoldi loop, undoing the map at the end). Under
! REVERSE the two must swap roles and apply their transposes; until
! that is implemented, converge rejects a preconditioned REVERSE solve
! with a clear error rather than silently mis-preconditioning.
!
! One march override exists in the family: cg's newton/bdf linearized
! path, where the linearization (coefficients and an external
! right-hand side) defines a different frozen operator. That state
! belongs on the system and moves there in the linearization commit;
! until then the override is the documented exception.
!
! linear_solver extends the common marcher base: march (the marcher's
! deferred contract) is implemented by converge, and solve is the
! family's provided wrapper over a plain solution vector.
!
! Author : Komahan Boopathy
!=====================================================================!

module interface_linear_solver

  use iso_fortran_env           , only : dp => REAL64
  use interface_graph           , only : graph, counting_sort
  use interface_marcher         , only : marcher
  use interface_assembler       , only : assembler, DIAGONAL
  use interface_state           , only : state
  use class_differential_state  , only : differential_state
  use module_solve_mode         , only : FORWARD, REVERSE, is_valid_mode

  implicit none

  private
  public :: linear_solver
  public :: preconditioner
  public :: MANUAL, AUTO

  ! Acceptance tolerance of the transpose-consistency gate. The verify
  ! returns a relative defect (scaled by max(|lhs|, |rhs|, 1)), so this
  ! is dimensionless and sits just above machine precision.
  real(dp), parameter :: transpose_defect_tol = 1.0d-12

  ! tuning modes: MANUAL takes the parameters as given (default); AUTO lets
  ! the solver select its parameters itself (opt-in, never default)
  integer, parameter :: MANUAL = 0
  integer, parameter :: AUTO   = 1

  !===================================================================!
  ! Abstract preconditioner: applies an approximate inverse z = M^-1 r -
  ! a cheap approximation of the residual -> correction map. A solver's
  ! kernel calls apply where its algorithm requires; concrete
  ! preconditioners (algebraic multigrid, jacobi, ...) extend this.
  !===================================================================!

  type, abstract :: preconditioner
   contains
     procedure(apply_interface), deferred :: apply
  end type preconditioner

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!

  type, abstract, extends(marcher) :: linear_solver

     ! max_tol, max_it and print_level are inherited from marcher

     ! iteration-history trace written when print_level == -1
     ! (machine-readable columns, not prose)
     character(len=:), allocatable :: res_file

     ! the two operator slots (storage; kernels apply them)
     class(preconditioner), allocatable :: pre_conditioner
     class(preconditioner), allocatable :: post_conditioner

     ! MANUAL (default) or AUTO - see tune below
     integer :: tuning = MANUAL

     ! inner iterations accumulated over the most recent march - the
     ! diagnostic lives on the object (the caller owns the state; no
     ! module-scope side channels)
     integer :: last_inner_iters = 0

     ! a carried graph and its coloring: when a stationary solver
     ! carries the system's graph, its sweep goes color by color
     ! (see colored_sweep). color_list holds the dofs grouped by
     ! color; color c's dofs are color_list(color_ptr(c) : ptr(c+1)-1)
     class(graph), allocatable :: g
     integer                   :: ncolors = 0
     integer     , allocatable :: color_ptr(:)
     integer     , allocatable :: color_list(:)

   contains

     ! march (the marcher contract): the residual-minimization iteration,
     ! generic over the state. converge is also bound by name so the one
     ! documented override can delegate to it; solve is the family
     ! wrapper over a plain solution vector.
     procedure :: march => converge
     procedure :: converge
     procedure :: solve

     ! deferred: the one-step correction (the analogue of the
     ! integrator's step)
     procedure(iterate_interface), deferred :: iterate

     ! passthrough apply helpers for the two slots
     procedure :: apply_pre_conditioner
     procedure :: apply_post_conditioner

     ! carry a graph (color it, group the dofs), and the exact
     ! stationary sweep the coloring buys
     procedure :: carry
     procedure :: colored_sweep

     ! parameter-selection hook, called by converge when tuning == AUTO:
     ! once at entry (pass 0 - measure the convergence factor, set the
     ! parameter) and once per pass with the measured rate (analogous to
     ! adaptive time-step control). three requirements: AUTO is opt-in,
     ! the residual rate is the progress measurement, and a parameter
     ! change that worsens the rate is reverted and discarded. default:
     ! no parameters to select - a solver with a selection rule
     ! overrides.
     procedure :: tune

     ! the solver's objective measurement: largest eigenvalue size of the
     ! operator by power iteration (products and norms only)
     procedure :: estimate_spectral_radius

     ! shared reporting, residual-based
     procedure :: residual_norm
     procedure :: monitor_step

  end type linear_solver

  !===================================================================!
  ! Deferred interfaces
  !===================================================================!

  abstract interface

     ! One correction step: drive A dx = r starting from dx = 0, where
     ! A is the system's operator and r the residual handed in by the
     ! outer iteration. iter reports the kernel's iterations for the
     ! monitor. this is intent(inout): the honest signature for a
     ! solver that keeps notes on itself.
     impure subroutine iterate_interface(this, system, r, dx, iter)
       import linear_solver
       import assembler
       import dp
       class(linear_solver), intent(inout) :: this
       class(assembler)    , intent(in)    :: system
       real(dp)            , intent(in)    :: r(:)
       real(dp)            , intent(out)   :: dx(:)
       integer             , intent(out)   :: iter
     end subroutine iterate_interface

     ! z = M^-1 r  (the approximate-inverse action)
     subroutine apply_interface(this, r, z)
       import preconditioner
       import dp
       class(preconditioner), intent(in)  :: this
       real(dp)             , intent(in)  :: r(:)
       real(dp)             , intent(out) :: z(:)
     end subroutine apply_interface

  end interface

contains

  !===================================================================!
  ! The provided solve: evaluate the residual, stop when it has dropped
  ! below the tolerance relative to its initial value, otherwise apply
  ! the correction step and update the solution. When the residual does
  ! not depend on the solution this reduces to a single solve plus a
  ! confirming pass; when it does (a non-orthogonal mesh), the
  ! iteration converges the true residual.
  !===================================================================!

  impure subroutine converge(this, system, s, mode)

    class(linear_solver)  , intent(inout)        :: this
    class(assembler)      , intent(inout)        :: system
    class(state)          , intent(inout)        :: s
    integer               , intent(in), optional :: mode  ! FORWARD (default) / REVERSE

    real(dp), allocatable :: r(:), dx(:), x(:)
    real(dp) :: rnorm, rnorm0, rnorm_prev, dxnorm
    integer  :: pass, inner

    ! Reject a preconditioned REVERSE solve: under the transpose the two
    ! slots swap roles and must apply their transposes, and that is not
    ! implemented yet. Failing loudly here prevents a silently
    ! mis-preconditioned adjoint solve.
    if (present(mode)) then
       if (mode .eq. REVERSE .and. &
            & (allocated(this % pre_conditioner) .or. &
            &  allocated(this % post_conditioner))) then
          error stop "linear_solver: preconditioned REVERSE solve refused - " // &
               & "the swap-and-transpose contract is not implemented"
       end if
    end if

    ! a wrong tag dies at the door with its name
    if (present(mode)) then
       if (.not. is_valid_mode(mode)) then
          write(*,'(1x,a,i0)') "converge: invalid mode tag ", mode
          error stop "converge: mode must be FORWARD or REVERSE"
       end if
    end if

    ! A REVERSE march through this shared iteration is valid only on a
    ! declared-symmetric operator: the loop below evaluates the forward
    ! residual and forward sweeps, which IS the adjoint solve exactly
    ! when J^T = J. A genuinely non-symmetric REVERSE solve needs the
    ! direction threaded through the residual and the sweep - deferred
    ! and tracked in the register; refusing here prevents certifying a
    ! claim the iteration would then ignore.
    if (present(mode)) then
       if (mode .eq. REVERSE .and. .not. system % operator_is_symmetric) then
          error stop "converge: a REVERSE march through the shared iteration requires " // &
               & "a declared-symmetric operator (the loop is direction-blind); a " // &
               & "genuine transpose march is a tracked deferral"
       end if
    end if

    ! Verify-before: a declared symmetry claim is checked against the
    ! analytic identity <w, J v> = <J^T w, v> before the FIRST march of
    ! any direction - forward marches consume the claim too (the
    ! normal-equation sweeps issue transpose products from inside a
    ! forward march). Once per system, a handful of products, machine
    ! precision, no truncation error. The comparison is written so a NaN
    ! defect also refuses (NaN fails every ordered comparison).
    if (system % operator_is_symmetric .and. .not. system % transpose_verified) then
       verify_transpose: block
         real(dp) :: defect
         defect = system % verify_transpose_consistency()
         if (.not. (defect .le. transpose_defect_tol)) then
            write(*,'(1x,a,es12.5)') &
                 & "converge: transpose-consistency defect ", defect
            error stop "converge: the system's transpose claim failed verification"
         end if
         system % transpose_verified = .true.
       end block verify_transpose
    end if

    ! parameter selection at entry (AUTO only)
    if (this % tuning .eq. AUTO) call this % tune(system, 0, 1.0_dp)

    allocate(x(system % num_state_vars))
    allocate(r, dx, mold = x)

    if (this % print_level .eq. -1) then
       open(13, file = history_file(this), action = 'write')
       write(13,'(a)') "# pass  inner  rnorm  rnorm_rel  dxnorm"
    end if

    rnorm0 = 1.0_dp
    dxnorm = 0.0_dp
    inner  = 0

    this % last_inner_iters = 0

    outer_iterations: do pass = 1, this % max_it + 1

       ! the residual at the state's current solution values
       x = current_solution(s)
       call system % get_residual(r, x)
       rnorm = sqrt(system % inner_product(r, r))

       if (pass .eq. 1) then
          rnorm0 = rnorm
          if (rnorm0 .le. tiny(1.0_dp)) then
             write(*,*) "warning: zero right-hand side; the solution is x = 0"
             exit outer_iterations
          end if
       end if

       if (this % print_level .eq. -1) then
          write(13,'(2(1x,i6),3(1x,es22.14))') pass, inner, rnorm, rnorm/rnorm0, dxnorm
       end if
       if (this % print_level .gt. 0) then
          call this % monitor_step(pass, inner, rnorm, rnorm0, dxnorm, &
               & sqrt(system % inner_product(x, x)))
       end if

       ! termination on the true residual
       if (rnorm .le. this % max_tol * rnorm0) exit outer_iterations

       ! parameter adaptation from the measured rate (AUTO only)
       if (this % tuning .eq. AUTO .and. pass .gt. 1) then
          call this % tune(system, pass, rnorm/rnorm_prev)
       end if

       ! stagnation: the sweep can no longer reduce the residual (its
       ! floor - typically machine precision). this is the criterion the
       ! old update-norm loop applied implicitly.
       if (pass .gt. 1) then
          if (rnorm .ge. 0.999_dp*rnorm_prev) then
             if (this % print_level .gt. 0) then
                write(*,'(1x,a,es12.5)') &
                     & "converge: stagnated at the smallest residual the sweep can attain, ||r||/||r0|| = ", rnorm/rnorm0
             end if
             exit outer_iterations
          end if
       end if
       rnorm_prev = rnorm

       if (pass .gt. this % max_it) then
          write(*,*) "warning: converge hit max_it without meeting tolerance; ||r||/||r0|| =", &
               & rnorm/rnorm0
          exit outer_iterations
       end if

       ! one correction step: A dx = r from dx = 0; the state applies it
       call this % iterate(system, r, dx, inner)
       this % last_inner_iters = this % last_inner_iters + inner

       call s % update(dx)
       dxnorm = sqrt(system % inner_product(dx, dx))

    end do outer_iterations

    if (this % print_level .eq. -1) close(13)

  end subroutine converge

  !===================================================================!
  ! The family wrapper: solve for a plain solution vector. Wraps x in
  ! an order-0 state (its update rule is u <- u + du), marches, and
  ! copies the values back out. Every existing call site keeps this
  ! signature.
  !===================================================================!

  impure subroutine solve(this, system, x, mode)

    class(linear_solver)  , intent(inout)        :: this
    class(assembler)      , intent(inout)        :: system
    real(dp), allocatable , intent(out)          :: x(:)
    integer               , intent(in), optional :: mode  ! FORWARD (default) / REVERSE

    type(differential_state) :: s

    s = differential_state(system % num_state_vars, 0)

    call this % march(system, s, mode)

    x = current_solution(s)

  end subroutine solve

  !===================================================================!
  ! The solution values of a state (its zeroth derivative order)
  !===================================================================!

  function current_solution(s) result(x)

    class(state), intent(in) :: s
    real(dp), allocatable    :: x(:)

    type(scalar), allocatable :: v(:,:)

    v = s % values()
    x = real(v(:,1), dp)

  end function current_solution

  !===================================================================!
  ! Apply the pre-operator z = M^-1 r, or pass through if none attached
  !===================================================================!

  impure subroutine apply_pre_conditioner(this, r, z)

    class(linear_solver), intent(in)  :: this
    real(dp)            , intent(in)  :: r(:)
    real(dp)            , intent(out) :: z(:)

    if (allocated(this % pre_conditioner)) then
       call this % pre_conditioner % apply(r, z)
    else
       z = r
    end if

  end subroutine apply_pre_conditioner

  !===================================================================!
  ! Apply the post-operator z = M^-1 r, or pass through if none attached
  !===================================================================!

  impure subroutine apply_post_conditioner(this, r, z)

    class(linear_solver), intent(in)  :: this
    real(dp)            , intent(in)  :: r(:)
    real(dp)            , intent(out) :: z(:)

    if (allocated(this % post_conditioner)) then
       call this % post_conditioner % apply(r, z)
    else
       z = r
    end if

  end subroutine apply_post_conditioner

  !===================================================================!
  ! Carry a graph: keep a copy, color it, and group the dofs color by
  ! color (the counting kernel groups, dofs_of expands). Carry after
  ! the system is fully configured - the dof layout is read from the
  ! graph as handed in.
  !===================================================================!

  pure subroutine carry(this, g)

    class(linear_solver), intent(inout) :: this
    class(graph)        , intent(in)    :: g

    integer, allocatable :: colors(:), vptr(:), vlist(:)
    integer :: c, v, nv

    allocate(this % g, source = g)

    nv     = g % num_vertices
    colors = g % coloring()
    this % ncolors = 0
    if (nv .gt. 0) this % ncolors = maxval(colors)

    call counting_sort(this % ncolors, colors, [(v, v = 1, nv)], vptr, vlist)

    allocate(this % color_ptr(this % ncolors + 1))
    this % color_ptr(1) = 1
    do c = 1, this % ncolors
       this % color_ptr(c+1) = this % color_ptr(c) &
            &                + (vptr(c+1) - vptr(c))*g % num_variables
    end do
    allocate(this % color_list(this % color_ptr(this % ncolors + 1) - 1))
    do c = 1, this % ncolors
       this % color_list(this % color_ptr(c) : this % color_ptr(c+1) - 1) = &
            & g % dofs_of(vlist(vptr(c) : vptr(c+1) - 1))
    end do

  end subroutine carry

  !===================================================================!
  ! The colored sweep: the carried graph's coloring chops the
  ! vertices into independent sets - no edge stays inside a color -
  ! so a whole color updates at once, exactly:
  !
  !    1 - 2 - 1 - 2        for each color c:
  !    |   |   |   |           dx(c) += omega * (r - A dx)(c) / D(c)
  !    2 - 1 - 2 - 1
  !    |   |   |   |        nothing inside c couples, so this IS the
  !    1 - 2 - 1 - 2        triangular sweep in the color ordering -
  !                         no inner iteration approximating it
  !
  ! omega = 1 is the gauss-seidel sweep; omega /= 1 relaxes it (sor).
  !===================================================================!

  impure subroutine colored_sweep(this, system, r, dx, iter, omega)

    class(linear_solver), intent(inout) :: this
    class(assembler)    , intent(in)    :: system
    real(dp)            , intent(in)    :: r(:)
    real(dp)            , intent(out)   :: dx(:)
    integer             , intent(out)   :: iter
    real(dp)            , intent(in)    :: omega

    real(dp), allocatable :: Adx(:), D(:), dxold(:), identity(:)
    real(dp) :: tol, bnorm
    integer  :: c

    dx = 0.0_dp
    allocate(Adx, D, dxold, identity, mold = dx)

    ! the diagonal (the self-loop weights)
    identity = 1.0_dp
    call system % get_jacobian_residual_product(D, identity, part = DIAGONAL)

    bnorm = sqrt(system % inner_product(r, r))
    if (bnorm .le. this % max_tol) then
       iter = 0
       return
    end if

    iter = 1; tol = huge(1.0_dp)
    do while (tol .gt. this % max_tol .and. iter .lt. this % max_it)

       dxold = dx
       do c = 1, this % ncolors
          associate(cd => this % color_list(this % color_ptr(c) : this % color_ptr(c+1) - 1))
            call system % get_jacobian_residual_product(Adx, dx)
            dx(cd) = dx(cd) + omega*(r(cd) - Adx(cd))/D(cd)
          end associate
       end do

       tol = sqrt(system % inner_product(dx - dxold, dx - dxold))
       if (this % print_level .gt. 1) write(*,*) "inner (colored)", iter, tol
       iter = iter + 1

    end do

  end subroutine colored_sweep

  !===================================================================!
  ! Norm of the system's residual at x
  !===================================================================!

  impure function residual_norm(this, sys, x) result(rnorm)

    class(linear_solver), intent(in) :: this
    class(assembler)    , intent(in) :: sys
    real(dp)            , intent(in) :: x(:)

    real(dp)              :: rnorm
    real(dp), allocatable :: r(:)

    allocate(r, mold = x)
    call sys % get_residual(r, x)
    rnorm = sqrt(sys % inner_product(r, r))

  end function residual_norm

  !===================================================================!
  ! One row of the convergence table: the residual (absolute and as a
  ! drop from the initial one), the correction size (absolute and
  ! relative to the answer). Prints the header on the first pass.
  !===================================================================!

  impure subroutine monitor_step(this, pass, inner, rnorm, rnorm0, dxnorm, xnorm)

    class(linear_solver), intent(in) :: this
    integer             , intent(in) :: pass
    integer             , intent(in) :: inner
    real(dp)            , intent(in) :: rnorm
    real(dp)            , intent(in) :: rnorm0
    real(dp)            , intent(in) :: dxnorm
    real(dp)            , intent(in) :: xnorm

    real(dp) :: r_drop, dx_rel

    r_drop = rnorm
    if (rnorm0 .gt. epsilon(1.0_dp)) r_drop = rnorm/rnorm0

    dx_rel = dxnorm
    if (xnorm .gt. epsilon(1.0_dp)) dx_rel = dxnorm/xnorm

    if (pass .eq. 1) then
       write(*,'(2x,a5,2x,a5,4(2x,a13))') &
            & "pass", "inner", "||r||", "||r||/||r0||", "||dx||", "||dx||/||x||"
       write(*,'(2x,a5,2x,a5,4(2x,a13))') &
            & "----", "-----", "-------------", "-------------", &
            & "-------------", "-------------"
    end if

    write(*,'(2x,i5,2x,i5,4(2x,es13.5))') pass, inner, rnorm, r_drop, dxnorm, dx_rel

  end subroutine monitor_step

  !===================================================================!
  ! Parameter-selection hook (see the binding comment). Default: this
  ! solver has no parameters to select.
  !===================================================================!

  impure subroutine tune(this, system, pass, rate)

    class(linear_solver), intent(inout) :: this
    class(assembler)    , intent(in)    :: system
    integer             , intent(in)    :: pass   ! 0 = at entry (static)
    real(dp)            , intent(in)    :: rate   ! ||r_k||/||r_{k-1}|| when pass > 1

  end subroutine tune

  !===================================================================!
  ! Largest eigenvalue size of the operator by power iteration - the
  ! solver-level objective measurement, built from products and norms only.
  !===================================================================!

  impure subroutine estimate_spectral_radius(this, system, mu, max_iter)

    class(linear_solver), intent(in)  :: this
    class(assembler)    , intent(in)  :: system
    real(dp)            , intent(out) :: mu
    integer             , intent(in)  :: max_iter

    real(dp), allocatable :: v(:), w(:)
    integer  :: iter
    real(dp) :: wnorm

    call system % create_vector(v)
    call random_number(v)
    v = v/sqrt(system % inner_product(v, v))

    call system % create_vector(w)

    power_iteration: do iter = 1, max_iter
       call system % get_jacobian_residual_product(w, v)
       wnorm = sqrt(system % inner_product(w, w))
       v  = w/wnorm
       mu = system % inner_product(v, w)
       if (this % print_level .gt. 1) write(*,*) iter, mu
    end do power_iteration

    deallocate(v, w)

  end subroutine estimate_spectral_radius

  !===================================================================!
  ! Name of the iteration-history trace file (print_level == -1)
  !===================================================================!

  pure function history_file(this) result(fname)

    class(linear_solver), intent(in) :: this
    character(len=:), allocatable    :: fname

    if (allocated(this % res_file)) then
       fname = this % res_file
    else
       fname = 'solver.res'
    end if

  end function history_file

end module interface_linear_solver
