!=====================================================================!
! Conjugate Gradient linear solver class that uses the functionalities of assembler class
! in the iterative solution process.
!=====================================================================!

module class_conjugate_gradient

  use iso_fortran_env         , only : dp => REAL64
  use interface_linear_solver , only : linear_solver, preconditioner
  use interface_assembler     , only : assembler
  use module_solve_mode       , only : FORWARD, REVERSE

  implicit none

  ! Expose only the linear solver datatype
  private
  public :: conjugate_gradient
  public :: cg_last_iters          ! inner CG iterations of the last solve (diagnostics)

  ! Accumulated inner CG iterations of the most recent solve. Written by
  ! solve/iterate (which take `this` as intent(in), so this cannot live on
  ! the object); read by tests comparing plain-CG vs PCG-AMG iteration counts.
  integer :: cg_last_iters = 0

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!

  type, extends(linear_solver) :: conjugate_gradient

     ! stateless w.r.t. the system: solve takes the assembler as an argument
     class(preconditioner), allocatable :: precond     ! optional; absent => plain CG
     integer                            :: print_level

     ! Newton/BDF linearized mode: when lin_coeff is set, solve drives the
     ! matrix-free system  J dq = rhs,  J v = add_jacobian_vector_product(.,.,
     ! lin_coeff) (transpose in REVERSE), with rhs the external rhs if set
     ! else the system residual -R. Unset => the steady A x = b path below.
     real(dp), allocatable :: lin_coeff(:)
     real(dp), allocatable :: rhs(:)

   contains

     ! type bound procedures
     procedure :: solve
     procedure :: iterate
     procedure, private :: apply_precond
     procedure, private :: solve_linearized

     ! destructor
     final :: destroy

  end type conjugate_gradient

  !===================================================================!
  ! Interface for multiple constructors
  !===================================================================!

  interface conjugate_gradient
     module procedure construct
  end interface conjugate_gradient

contains

  !===================================================================!
  ! Constructor for linear solver
  !===================================================================!

  pure type(conjugate_gradient) function construct(max_it, &
       & max_tol, print_level, precond) result (this)

    type(integer)        , intent(in)           :: max_it
    type(real(dp))       , intent(in)           :: max_tol
    type(integer)        , intent(in)           :: print_level
    class(preconditioner), intent(in), optional :: precond

    this % max_it      = max_it
    this % max_tol     = max_tol
    this % print_level = print_level

    ! optional preconditioner (absent => plain CG, bit-identical to before)
    if (present(precond)) allocate(this % precond, source = precond)

  end function construct

  !===================================================================!
  ! Destructor for linear solver
  !===================================================================!

  pure subroutine destroy(this)

    type(conjugate_gradient), intent(inout) :: this

    if (allocated(this % precond)) deallocate(this % precond)

  end subroutine destroy

  !===================================================================!
  ! Iterative linear solution using conjugate gradient method
  !===================================================================!

  impure subroutine solve(this, system, x, mode)

    class(conjugate_gradient), intent(in)       :: this
    class(assembler)         , intent(in)       :: system
    real(dp), allocatable    , intent(out)      :: x(:)
    integer              , intent(in), optional :: mode  ! FORWARD (default) / REVERSE

    ! Locals
    real(dp), allocatable :: xold(:), ss(:)
    real(dp) :: update_norm, rnorm0
    integer  :: iter, num_inner_iters

    ! reset the inner-iteration counter for this solve
    cg_last_iters = 0

    ! Newton/BDF linearized inner solve (J dq = rhs) takes a separate path
    if (allocated(this % lin_coeff)) then
       call this % solve_linearized(system, x, mode)
       return
    end if

    ! Initial guess vector for the subspace is "b"
    allocate(x(system % num_state_vars))
    call system % get_source(x)
    if (norm2(x) .lt. epsilon(1.0_dp)) then
       print *, 'zero rhs?'
       ! error stop
    end if

    allocate(xold(system % num_state_vars))
    allocate(ss(system % num_state_vars))

    if (this % print_level .eq. -1) then
       open(13, file='cg.res', action='write')
       write(13,*) "iteration ", " residual"
    end if

    rnorm0 = 0.0_dp
    if (this % print_level .gt. 0) rnorm0 = this % residual_norm(system, x)

    update_norm = huge(1.0d0);
    iter = 1;
    outer_iterations: do while (update_norm .gt. this % max_tol .and. iter .le. this % max_it)

       xold = x

       ! Inner iterations with CG
       if ( iter .eq. 1) then
          ss = 0.0d0
          call this % iterate(system, x, ss, num_inner_iters)
       else
          call system % get_skew_source(ss, x)
          call this % iterate(system, x, ss, num_inner_iters)
       end if

       update_norm = norm2(x - xold)
       if (this % print_level .eq. -1) then
          write(13, *) iter, update_norm
       end if
       if (this % print_level .gt. 0) then
          call this % monitor_step(system, iter, num_inner_iters, x, xold, rnorm0)
       end if
       iter = iter + 1

    end do outer_iterations

    if (this % print_level .eq. -1) then
       close(13)
    end if

    ! Safety net: report if the outer (deferred-correction) loop was
    ! capped at max_it without meeting the tolerance.
    if (update_norm .gt. this % max_tol) then
       write(*,*) "warning: outer correction loop hit max_it without convergence; update_norm =", update_norm
    end if

    deallocate(xold)
    deallocate(ss)

  end subroutine solve

  !===================================================================!
  ! Iterative linear solution using conjugate gradient method
  !===================================================================!

  impure subroutine iterate(this, system, x, ss, iter)

    class(conjugate_gradient) , intent(in)    :: this
    class(assembler)          , intent(in)    :: system
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(in)    :: ss(:)
    integer                   , intent(out)   :: iter

    ! Create local data
    real(dp), allocatable :: p(:), r(:), w(:), Ax(:), tmp(:), z(:)
    real(dp), allocatable :: b(:)
    real(dp)              :: alpha, beta
    real(dp)              :: bnorm, rnorm
    real(dp)              :: tol
    real(dp)              :: rho(2)


    ! Start the iteration counter
    iter = 1

    ! Memory allocations. z holds the preconditioned residual M^-1 r
    ! (z = r when no preconditioner is attached -> reduces to plain CG).
    allocate(b,p,r,w,Ax,tmp,z,mold=x)

    ! Norm of the right hand side
    call system % get_source(tmp)
    ! Add the additional right hand side supplied
    tmp = tmp + ss
    b = tmp
    bnorm = norm2(b)

    ! Homogeneous case
    if (bnorm .le. this % max_tol) then
       x = 0.0d0
       return
    end if

    ! Norm of the initial residual
    call system % get_jacobian_vector_product(tmp, x)
    Ax = tmp
    r         = b - Ax ! could directly form this residual using get_residual_call
    rnorm     = norm2(r)
    tol       = rnorm/bnorm
    ! preconditioned residual z = M^-1 r; rho = <r, z> (PCG). z = r when
    ! unpreconditioned, so rho = <r, r> and this is exactly plain CG.
    call this % apply_precond(r, z)
    rho(2)    = dot_product(r, z)

    !open(13, file='cg.log', action='write', position='append')

    ! Apply Iterative scheme until tolerance is achieved
    do while ((tol .gt. this % max_tol) .and. (iter .lt. this % max_it))

       ! step (a) compute the descent direction (on the preconditioned residual)
       if ( iter .eq. 1) then
          ! steepest descent direction p
          p = z
       else
          ! take a conjugate direction
          beta = rho(2)/rho(1)
          p = z + beta*p
       end if

       ! step (b) compute the solution update
       call system % get_jacobian_vector_product(tmp, p)
       w = tmp

       ! step (c) compute the step size for update
       alpha = rho(2)/dot_product(p, w)

       ! step (d) Add dx to the old solution
       x = x + alpha*p

       ! step (e) compute the new residual (true residual drives the tolerance)
       r = r - alpha*w

       ! step(f) update values before next iteration
       rnorm = norm2(r)
       tol = rnorm/bnorm

       ! write(13,*) iter, tol
       if (this % print_level .gt. 1) then
          write(*,*) iter, tol, rnorm, rho ! causes valgrind errors
       end if

       iter = iter + 1

       call this % apply_precond(r, z)
       rho(1) = rho(2)
       rho(2) = dot_product(r, z)

    end do

    !close(13)

    ! record inner iterations done (iter starts at 1; iter-1 = CG steps)
    cg_last_iters = cg_last_iters + (iter - 1)

    deallocate(r, p, w, b, Ax, tmp, z)

  end subroutine iterate

  !===================================================================!
  ! Apply the preconditioner z = M^-1 r, or z = r if none is attached
  !===================================================================!

  impure subroutine apply_precond(this, r, z)

    class(conjugate_gradient), intent(in)  :: this
    real(dp)                 , intent(in)  :: r(:)
    real(dp)                 , intent(out) :: z(:)

    if (allocated(this % precond)) then
       call this % precond % apply(r, z)
    else
       z = r
    end if

  end subroutine apply_precond

  !===================================================================!
  ! Matrix-free CG for the linearized system  J dq = rhs, where
  ! J v = add_jacobian_vector_product(., v, lin_coeff) (its transpose in
  ! REVERSE for the adjoint). rhs is the external rhs member if set (the
  ! adjoint's -df/du), else the system residual -R at the current state.
  ! Drives the increment dq (x0 = 0); the Newton/BDF outer loop owns the
  ! state update. This is the inner solve the nonlinear solver delegates.
  !===================================================================!

  impure subroutine solve_linearized(this, system, x, mode)

    class(conjugate_gradient), intent(in)       :: this
    class(assembler)         , intent(in)       :: system
    real(dp), allocatable    , intent(out)      :: x(:)
    integer              , intent(in), optional :: mode

    real(dp), allocatable :: r(:), p(:), Jp(:), b(:), res(:)
    real(dp)              :: rs_old, rs_new, alpha, pJp
    integer               :: nvars, it, max_it, dir

    dir = FORWARD
    if (present(mode)) dir = mode

    nvars  = system % get_num_state_vars()
    max_it = nvars + 100

    allocate(x(nvars), r(nvars), p(nvars), Jp(nvars), b(nvars))

    ! right-hand side: external (adjoint -df/du) if set, else -residual
    if (allocated(this % rhs)) then
       b = this % rhs
    else
       allocate(res(nvars))
       res = 0.0_dp
       call system % add_residual(res)
       b = -res
    end if

    x = 0.0_dp
    r = b
    p = r
    rs_old = dot_product(r, r)

    cg: do it = 1, max_it

       if (norm2(r) .le. 1.0d-14) exit cg

       Jp = 0.0_dp
       if (dir .eq. REVERSE) then
          call system % add_jacobian_vector_product_transpose(Jp, p, this % lin_coeff)
       else
          call system % add_jacobian_vector_product(Jp, p, this % lin_coeff)
       end if

       pJp = dot_product(p, Jp)
       if (abs(pJp) .le. tiny(1.0_dp)) exit cg

       alpha = rs_old/pJp

       x = x + alpha*p
       r = r - alpha*Jp

       rs_new = dot_product(r, r)
       p      = r + (rs_new/rs_old)*p
       rs_old = rs_new

    end do cg

    cg_last_iters = it

  end subroutine solve_linearized

end module class_conjugate_gradient
