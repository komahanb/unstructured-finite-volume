!=====================================================================!
! Conjugate-gradient linear solver: supplies the sweep (`iterate`) -
! one preconditioned cg pass driving the correction equation A dx = r -
! and inherits the residual-minimization march from linear_solver. Its
! optional preconditioner fills the inherited pre_conditioner slot
! (applied inside every cg iteration).
!
! One documented override: the newton/bdf linearized path (lin_coeff
! set) solves the different frozen system J dq = rhs; that state moves
! to the system in the linearization commit.
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

     ! Newton/BDF linearized mode: when lin_coeff is set, solve drives the
     ! matrix-free system  J dq = rhs,  J v = add_jacobian_vector_product(.,.,
     ! lin_coeff) (transpose in REVERSE), with rhs the external rhs if set
     ! else the system residual -R. Unset => the inherited march below.
     real(dp), allocatable :: lin_coeff(:)
     real(dp), allocatable :: rhs(:)

   contains

     ! type bound procedures
     procedure :: solve
     procedure :: iterate
     procedure, private :: solve_linearized

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
    this % res_file    = 'cg.res'

    ! optional preconditioner fills the pre-operator slot (plain CG without)
    if (present(precond)) allocate(this % pre_conditioner, source = precond)

  end function construct

  !===================================================================!
  ! Solve. Resets the iteration counter, takes the newton/bdf linearized
  ! path when lin_coeff is set (the documented override), and otherwise
  ! runs the inherited residual-minimization march around the cg sweep.
  !===================================================================!

  impure subroutine solve(this, system, x, mode)

    class(conjugate_gradient), intent(in)       :: this
    class(assembler)         , intent(in)       :: system
    real(dp), allocatable    , intent(out)      :: x(:)
    integer              , intent(in), optional :: mode  ! FORWARD (default) / REVERSE

    ! reset the inner-iteration counter for this solve
    cg_last_iters = 0

    ! Newton/BDF linearized inner solve (J dq = rhs) takes a separate path
    if (allocated(this % lin_coeff)) then
       call this % solve_linearized(system, x, mode)
       return
    end if

    call this % converge(system, x, mode)

  end subroutine solve

  !===================================================================!
  ! The sweep: one (preconditioned) conjugate-gradient pass driving the
  ! correction equation A dx = r from dx = 0.
  !===================================================================!

  impure subroutine iterate(this, system, r, dx, iter)

    class(conjugate_gradient) , intent(in)  :: this
    class(assembler)          , intent(in)  :: system
    real(dp)                  , intent(in)  :: r(:)
    real(dp)                  , intent(out) :: dx(:)
    integer                   , intent(out) :: iter

    ! Local data. z holds the preconditioned residual M^-1 res
    ! (z = res when no pre-operator is attached -> plain CG).
    real(dp), allocatable :: p(:), res(:), w(:), z(:)
    real(dp)              :: alpha, beta
    real(dp)              :: bnorm, rnorm
    real(dp)              :: tol
    real(dp)              :: rho(2)

    iter = 1
    dx   = 0.0_dp

    allocate(p, res, w, z, mold = dx)

    ! the imbalance handed in is the right-hand side of the correction
    ! equation; at dx = 0 it is also the whole residual
    res   = r
    bnorm = sqrt(system % inner_product(res, res))

    ! Homogeneous case
    if (bnorm .le. this % max_tol) then
       iter = 0
       return
    end if

    rnorm = bnorm
    tol   = 1.0_dp

    ! preconditioned residual z = M^-1 res; rho = <res, z> (PCG). z = res
    ! when unpreconditioned, so rho = <res, res> and this is plain CG.
    call this % apply_pre_conditioner(res, z)
    rho(2) = system % inner_product(res, z)

    ! Apply the iterative scheme until tolerance is achieved
    do while ((tol .gt. this % max_tol) .and. (iter .lt. this % max_it))

       ! step (a) descent direction (on the preconditioned residual)
       if ( iter .eq. 1) then
          p = z
       else
          beta = rho(2)/rho(1)
          p = z + beta*p
       end if

       ! step (b) the operator action on the direction
       call system % get_jacobian_residual_product(w, p)

       ! step (c) the step size
       alpha = rho(2)/system % inner_product(p, w)

       ! step (d) update the correction
       dx = dx + alpha*p

       ! step (e) the new residual (the true residual drives the tolerance)
       res = res - alpha*w

       ! step (f) update values before the next iteration
       rnorm = sqrt(system % inner_product(res, res))
       tol   = rnorm/bnorm

       if (this % print_level .gt. 1) then
          write(*,*) iter, tol, rnorm, rho
       end if

       iter = iter + 1

       call this % apply_pre_conditioner(res, z)
       rho(1) = rho(2)
       rho(2) = system % inner_product(res, z)

    end do

    ! record inner iterations done (iter starts at 1; iter-1 = CG steps)
    cg_last_iters = cg_last_iters + (iter - 1)

    deallocate(res, p, w, z)

  end subroutine iterate

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
    rs_old = system % inner_product(r, r)

    cg: do it = 1, max_it

       if (sqrt(system % inner_product(r, r)) .le. 1.0d-14) exit cg

       Jp = 0.0_dp
       if (dir .eq. REVERSE) then
          call system % add_jacobian_vector_product_transpose(Jp, p, this % lin_coeff)
       else
          call system % add_jacobian_vector_product(Jp, p, this % lin_coeff)
       end if

       pJp = system % inner_product(p, Jp)
       if (abs(pJp) .le. tiny(1.0_dp)) exit cg

       alpha = rs_old/pJp

       x = x + alpha*p
       r = r - alpha*Jp

       rs_new = system % inner_product(r, r)
       p      = r + (rs_new/rs_old)*p
       rs_old = rs_new

    end do cg

    cg_last_iters = it

  end subroutine solve_linearized

end module class_conjugate_gradient
