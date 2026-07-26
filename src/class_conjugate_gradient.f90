!=====================================================================!
! The conjugate-gradient linear solver supplies the sweep (`iterate`)
! - one preconditioned cg pass driving the correction equation
! A dx = r - and inherits the residual-minimization iteration from
! linear_solver. Its optional preconditioner fills the inherited
! pre_conditioner slot (applied inside every cg iteration).
!
!    r --> z --> p --> A p --> alpha --> dx, res
!          ^                              |
!          '<------------ beta ----------'
!
! Matvecs and dots, nothing else. The newton/bdf linearized path that
! once lived here as a march override now lives on the system
! (linearize / clear_linearization): a frozen linear system answers
! the same two questions every system answers, so this solver cannot
! tell it apart from any other - and needs no override.
!=====================================================================!

module class_conjugate_gradient

  use iso_fortran_env         , only : dp => REAL64
  use interface_linear_solver , only : linear_solver, preconditioner
  use interface_assembler     , only : assembler

  implicit none

  ! Expose only the linear solver datatype.
  private
  public :: conjugate_gradient

  !===================================================================!
  ! The linear solver datatype.
  !===================================================================!

  type, extends(linear_solver) :: conjugate_gradient

   contains

     ! The sweep consumed by the inherited outer iteration.
     procedure :: iterate

  end type conjugate_gradient

  !===================================================================!
  ! The interface for multiple constructors.
  !===================================================================!

  interface conjugate_gradient
     module procedure construct
  end interface conjugate_gradient

contains

  !===================================================================!
  ! The constructor for the linear solver.
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

    ! The optional preconditioner fills the pre-operator slot (plain CG without one).
    if (present(precond)) allocate(this % pre_conditioner, source = precond)

  end function construct

  !===================================================================!
  ! The sweep: one (preconditioned) conjugate-gradient pass driving the
  ! correction equation A dx = r from dx = 0.
  !===================================================================!

  impure subroutine iterate(this, system, r, dx, iter)

    class(conjugate_gradient) , intent(inout) :: this
    class(assembler)          , intent(in)    :: system
    real(dp)                  , intent(in)    :: r(:)
    real(dp)                  , intent(out)   :: dx(:)
    integer                   , intent(out)   :: iter

    !-----------------------------------------------------------------!
    ! Local data. The vector z holds the preconditioned residual
    ! M^-1 res (z = res when no pre-operator is attached, which gives
    ! plain CG).
    !-----------------------------------------------------------------!

    real(dp), allocatable :: p(:), res(:), w(:), z(:)
    real(dp)              :: alpha, beta
    real(dp)              :: bnorm, rnorm
    real(dp)              :: tol
    real(dp)              :: rho(2)

    iter = 1
    dx   = 0.0_dp

    allocate(p, res, w, z, mold = dx)

    !-----------------------------------------------------------------!
    ! The residual handed in is the right-hand side of the correction
    ! equation; at dx = 0 it is also the whole residual.
    !-----------------------------------------------------------------!

    res   = r
    bnorm = sqrt(system % inner_product(res, res))

    ! The homogeneous case leaves nothing to do.
    if (bnorm .le. this % max_tol) then
       iter = 0
       return
    end if

    rnorm = bnorm
    tol   = 1.0_dp

    !-----------------------------------------------------------------!
    ! The preconditioned residual is z = M^-1 res and rho = <res, z>
    ! (PCG). When unpreconditioned z = res, so rho = <res, res> and
    ! this is plain CG.
    !-----------------------------------------------------------------!

    call this % apply_pre_conditioner(res, z)
    rho(2) = system % inner_product(res, z)

    ! Apply the iterative scheme until the tolerance is achieved.
    do while ((tol .gt. this % max_tol) .and. (iter .lt. this % max_it))

       ! Step (a): form the descent direction (on the preconditioned residual).
       if ( iter .eq. 1) then
          p = z
       else
          beta = rho(2)/rho(1)
          p = z + beta*p
       end if

       ! Step (b): apply the operator action on the direction.
       call system % get_jacobian_residual_product(w, p)

       ! Step (c): compute the step size.
       alpha = rho(2)/system % inner_product(p, w)

       ! Step (d): update the correction.
       dx = dx + alpha*p

       ! Step (e): form the new residual (the true residual drives the tolerance).
       res = res - alpha*w

       ! Step (f): update values before the next iteration.
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

    !-----------------------------------------------------------------!
    ! Report the CG steps taken (iter starts at 1, so iter - 1 equals
    ! the steps); the outer iteration accumulates this on the object.
    !-----------------------------------------------------------------!

    iter = iter - 1

    deallocate(res, p, w, z)

  end subroutine iterate

end module class_conjugate_gradient
