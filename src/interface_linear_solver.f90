!=====================================================================!
! The generic linear solver, plus the preconditioner contract it
! consumes. A linear solver takes a guess and improves it until the
! system's imbalance (the residual) is gone. It speaks matrix, vector
! and inner product only - it asks the system for the residual and the
! product, and never names a source or a correction scheme.
!
! The provided solve (bound to converge) is the residual-minimization
! march - one lap of the cycle state -> residual -> correction -> state
! per pass, terminating on the true imbalance. The deferred iterate is
! the one thing a concrete solver must supply: one application of the
! inverse arrow, driving A dx = r from dx = 0. This mirrors
! interface_integrator, where the provided integrate marches the
! deferred step.
!
! Two optional operator slots accelerate the kernels: pre_conditioner
! (acts on imbalances - the left preconditioner / smoother) and
! post_conditioner (acts on answers - the right preconditioner). The
! slots are storage, not outer-loop actions: each kernel applies its
! slot where its own algorithm requires (cg inside every iteration,
! gmres inside its arnoldi loop with the back-map at the end). Under
! REVERSE the two must swap roles and apply their transposes; until
! that contract is implemented, converge REFUSES a preconditioned
! REVERSE solve rather than silently mis-preconditioning.
!
! One solve override exists in the family: cg's newton/bdf linearized
! path, where the linearization (coefficients and an external
! right-hand side) defines a different frozen operator. That state
! belongs on the system and moves there in the linearization commit;
! until then the override is the documented exception, and solve
! carries one meaning - drive the system's imbalance to zero.
!
! linear_solver extends the common algebraic_solver base.
!
! Author : Komahan Boopathy
!=====================================================================!

module interface_linear_solver

  use iso_fortran_env           , only : dp => REAL64
  use interface_algebraic_solver, only : algebraic_solver
  use interface_assembler       , only : assembler
  use module_solve_mode         , only : FORWARD, REVERSE

  implicit none

  private
  public :: linear_solver
  public :: preconditioner

  !===================================================================!
  ! Abstract preconditioner: applies an approximate inverse z = M^-1 r -
  ! a cheap version of the residual -> correction arrow. A solver's
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

  type, abstract, extends(algebraic_solver) :: linear_solver

     real(dp) :: max_tol
     integer  :: max_it
     integer  :: print_level = 0

     ! iteration-history trace written when print_level == -1
     ! (machine-readable columns, not prose)
     character(len=:), allocatable :: res_file

     ! the two operator slots (storage; kernels apply them)
     class(preconditioner), allocatable :: pre_conditioner
     class(preconditioner), allocatable :: post_conditioner

   contains

     ! provided solve: the residual-minimization march. converge is also
     ! bound by name so the one documented override can delegate to it.
     procedure :: solve => converge
     procedure :: converge

     ! deferred: the one-step map (the twin of the integrator's step)
     procedure(iterate_interface), deferred :: iterate

     ! passthrough apply helpers for the two slots
     procedure :: apply_pre_conditioner
     procedure :: apply_post_conditioner

     ! shared reporting, residual-based
     procedure :: residual_norm
     procedure :: monitor_step

  end type linear_solver

  !===================================================================!
  ! Deferred interfaces
  !===================================================================!

  abstract interface

     ! One application of the inverse arrow: drive A dx = r starting
     ! from dx = 0, where A is the system's frozen matrix and r the
     ! imbalance handed in by the march. iter reports the kernel's
     ! iterations for the monitor.
     impure subroutine iterate_interface(this, system, r, dx, iter)
       import linear_solver
       import assembler
       import dp
       class(linear_solver), intent(in)  :: this
       class(assembler)    , intent(in)  :: system
       real(dp)            , intent(in)  :: r(:)
       real(dp)            , intent(out) :: dx(:)
       integer             , intent(out) :: iter
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
  ! The provided solve: measure the imbalance, stop when it has dropped
  ! below tolerance relative to where it started, otherwise apply the
  ! sweep to the correction equation and add the correction. When the
  ! system's residual does not depend on the answer this collapses to
  ! solve-once-and-confirm; when it does (a non-orthogonal mesh), the
  ! march converges the true imbalance.
  !===================================================================!

  impure subroutine converge(this, system, x, mode)

    class(linear_solver)  , intent(in)           :: this
    class(assembler)      , intent(in)           :: system
    real(dp), allocatable , intent(out)          :: x(:)
    integer               , intent(in), optional :: mode  ! FORWARD (default) / REVERSE

    real(dp), allocatable :: r(:), dx(:)
    real(dp) :: rnorm, rnorm0, rnorm_prev, dxnorm
    integer  :: pass, inner

    ! Refuse a preconditioned REVERSE solve: under the transpose the two
    ! slots swap roles and must apply their transposes, and that
    ! contract is not implemented yet. Silence is the one option not
    ! permitted.
    if (present(mode)) then
       if (mode .eq. REVERSE .and. &
            & (allocated(this % pre_conditioner) .or. &
            &  allocated(this % post_conditioner))) then
          error stop "linear_solver: preconditioned REVERSE solve refused - " // &
               & "the swap-and-transpose contract is not implemented"
       end if
    end if

    allocate(x(system % num_state_vars))
    x = 0.0_dp
    allocate(r, dx, mold = x)

    if (this % print_level .eq. -1) then
       open(13, file = history_file(this), action = 'write')
       write(13,'(a)') "# pass  inner  rnorm  rnorm_rel  dxnorm"
    end if

    rnorm0 = 1.0_dp
    dxnorm = 0.0_dp
    inner  = 0

    marching: do pass = 1, this % max_it + 1

       ! the imbalance at the current answer
       call system % get_residual(r, x)
       rnorm = sqrt(system % inner_product(r, r))

       if (pass .eq. 1) then
          rnorm0 = rnorm
          if (rnorm0 .le. tiny(1.0_dp)) then
             write(*,*) "warning: zero imbalance at the zero answer; x = 0"
             exit marching
          end if
       end if

       if (this % print_level .eq. -1) then
          write(13,'(2(1x,i6),3(1x,es22.14))') pass, inner, rnorm, rnorm/rnorm0, dxnorm
       end if
       if (this % print_level .gt. 0) then
          call this % monitor_step(pass, inner, rnorm, rnorm0, dxnorm, &
               & sqrt(system % inner_product(x, x)))
       end if

       ! termination on the true imbalance
       if (rnorm .le. this % max_tol * rnorm0) exit marching

       ! stagnation: the sweep can no longer reduce the imbalance (its
       ! floor - typically machine precision). this is the criterion the
       ! old update-norm loop applied implicitly.
       if (pass .gt. 1) then
          if (rnorm .ge. 0.999_dp*rnorm_prev) then
             if (this % print_level .gt. 0) then
                write(*,'(1x,a,es12.5)') &
                     & "converge: stagnated at the sweep's floor, ||r||/||r0|| = ", rnorm/rnorm0
             end if
             exit marching
          end if
       end if
       rnorm_prev = rnorm

       if (pass .gt. this % max_it) then
          write(*,*) "warning: converge hit max_it without meeting tolerance; ||r||/||r0|| =", &
               & rnorm/rnorm0
          exit marching
       end if

       ! one application of the inverse arrow: A dx = r from dx = 0
       call this % iterate(system, r, dx, inner)

       x      = x + dx
       dxnorm = sqrt(system % inner_product(dx, dx))

    end do marching

    if (this % print_level .eq. -1) close(13)

  end subroutine converge

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
  ! Norm of the system's imbalance at x
  !===================================================================!

  pure function residual_norm(this, sys, x) result(rnorm)

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
  ! One row of the convergence table: the imbalance (absolute and as a
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
