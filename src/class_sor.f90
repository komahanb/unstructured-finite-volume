!=====================================================================!
! Relaxed-triangle-sweep linear solver (traditionally: successive
! over-relaxation, SOR): supplies only the sweep (`iterate`) - the
! relaxed lower-triangle solve (D+wL)y = R is itself done iteratively -
! and inherits the residual-minimization march from linear_solver.
! omega is the pseudo-time step size of the march.
!=====================================================================!

module class_sor

  use iso_fortran_env         , only : dp => REAL64
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

   contains

     ! the sweep consumed by the inherited march
     procedure :: iterate
     procedure :: estimate_spectral_radius

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
       & max_tol, print_level) result (this)

    type(real(dp)) , intent(in) :: omega
    type(integer)  , intent(in) :: max_it
    type(real(dp)) , intent(in) :: max_tol
    type(integer)  , intent(in) :: print_level

    this % max_it      = max_it
    this % max_tol     = max_tol
    this % print_level = print_level
    this % omega       = omega
    this % res_file    = 'sor.res'

  end function construct

  !===================================================================!
  ! Estimate spectral radius using power iteration (maximum absolute
  ! eigen value)
  !===================================================================!

  impure subroutine estimate_spectral_radius(this, system, mu, max_iter)

    ! Arguments
    class(sor)      , intent(in)  :: this
    class(assembler), intent(in)  :: system
    real(dp)        , intent(out) :: mu
    integer         , intent(in)  :: max_iter

    ! Random vector
    real(dp) , allocatable :: v(:), w(:)
    integer  :: iter
    real(dp) :: wnorm

    ! Create a random unit vector
    call system % create_vector(v)
    call random_number(v)
    v = v/sqrt(system % inner_product(v, v))

    ! A temp vector for processing
    call system % create_vector(w)

    power_iteration: do iter = 1, max_iter
       call system % get_jacobian_residual_product(w, v)
       wnorm = sqrt(system % inner_product(w, w))
       v = w/wnorm
       mu = system % inner_product(v, w)
       print *, iter, mu
    end do power_iteration

    deallocate(v, w)

  end subroutine estimate_spectral_radius

  !===================================================================!
  ! The sweep: relaxed (omega) triangle sweep on the correction
  ! equation A dx = r from dx = 0. The relaxed lower-triangle solve
  ! (D+wL)y = R is itself iterative.
  !===================================================================!

  impure subroutine iterate(this, system, r, dx, iter)

    class(sor)          , intent(in)  :: this
    class(assembler)    , intent(in)  :: system
    real(dp)            , intent(in)  :: r(:)
    real(dp)            , intent(out) :: dx(:)
    integer             , intent(out) :: iter

    ! Locals
    real(dp) :: tol, bnorm
    real(dp) , allocatable :: Ux(:), D(:), Ddx(:), R2(:), xnew(:), identity(:)
    real(dp) , allocatable :: Ly(:), y(:), ynew(:)

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
