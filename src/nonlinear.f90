#include "scalar.fpp"

!=====================================================================!
! Module that implements solving a nonlinear problem using different
! methods such as Newton's method.
!=====================================================================!

module nonlinear_algebra

  ! import dependencies
  use iso_fortran_env           , only : dp => REAL64
  use dynamic_physics_interface , only : dynamics
  use utils                     , only : norm

  ! disable implicit datatypes
  implicit none

  ! Define constants used  
  real(dp) :: abs_res_tol          = 1.0d-14
  real(dp) :: rel_res_tol          = 1.0d-10

  integer  :: max_newton_iters     = 15
  integer  :: print_level          = 0
  
  public :: solve

  !-------------------------------------------------------------------!
  ! Interface for nonlinear solution problems
  !-------------------------------------------------------------------!

  interface solve
     module procedure newton_solve_condensed
  end interface solve

contains

  !==================================================================!
  ! Newton solve for condensed form of equations
  !==================================================================!
  
  subroutine newton_solve_condensed(system, coeff, t, U, X, rhs)

    class(dynamics) , intent(inout) :: system
    type(scalar)    , intent(in)    :: coeff(:)
    type(scalar)    , intent(in)    :: t
    type(scalar)    , intent(inout) :: U(:,:)
    type(scalar)   , intent(in)     :: X(:,:)
    type(scalar)   , intent(in), optional :: rhs(:)
   
    ! Norms for tracking progress
    real(dp)                                  :: abs_res_norm = 0
    real(dp)                                  :: rel_res_norm = 0
    real(dp)                                  :: init_norm    = 0

    ! Other Local variables
    type(scalar), allocatable, dimension(:)   :: res, dq
    type(scalar), allocatable, dimension(:,:) :: jac

    integer                                   :: n, nvars, jj
    logical                                   :: conv = .false.

    ! find the size of the linear system based on the calling object
    nvars = system % get_num_state_vars()

    if (.not. allocated(res)) allocate(res(nvars))
    if (.not. allocated(dq)) allocate(dq(nvars))
    if (system % sparse .eqv. .true.) then
       if (.not.allocated(jac)) allocate(jac(nvars,3)) ! system % get_bandwidth() #TODO
    else
       if (.not.allocated(jac)) allocate(jac(nvars,nvars))
    end if

    if ( print_level .ge. 1 ) then
       write(*,'(A11, 2A12)') "Newton-Iter", "|R|", "|R|/|R1|"
    end if

    newton: do n = 1, max_newton_iters

       res = 0.0d0
       call system % add_residual(res, U, X)
       
       ! Handle arbitrary rhs external to system
       if (present(rhs)) then
          res = res + rhs
       end if

       jac = 0.0d0
       if ( (system % approximate_jacobian .eqv. .true.) .and. (system % sparse .eqv. .false.) ) then
          call approximate_jacobian(system, jac, coeff, U, X)
       else
          call system % add_jacobian(jac, coeff, U, X)
       end if

       ! Find norm of the residual
       abs_res_norm = norm(res)
       if ( n .eq. 1) init_norm = abs_res_norm
       rel_res_norm = abs_res_norm/init_norm

       if ( print_level .eq. 2) then
          write(*, "(I10,2ES12.2)") n, abs_res_norm, rel_res_norm
       end if

       ! Check stopping
       if ((abs_res_norm .le. abs_res_tol) .or. (rel_res_norm .le. rel_res_tol)) then
          conv = .true.
          exit newton
       else if ((abs_res_norm .ne. abs_res_norm) .or. (rel_res_norm .ne. rel_res_norm) ) then
          conv = .false.
          exit newton
       end if
              
       ! Solve the linear system if residual is not below tolerance

!!$       
!!$       if (system % sparse .eqv. .true.) then
!!$          ! works currently only for tridiagonal systems
!!$          res = -res
!!$          call tdma(jac, res, dq)
!!$       else
!!$          dq = linsolve(jac, -res)
!!$       end if       

       ! Apply BCs on dq
       ! call system % apply_bc(dq)?
       ! print *, 'update', dq(1), dq(nvars)
       ! dq(1) = 0.0d0
       ! dq(nvars) = 0.0d0

       forall(jj=1:system % get_differential_order() + 1)
          U(jj,:) = U(jj,:) + coeff(jj)*dq
       end forall
              
    end do newton

    if (print_level .eq. 1) then 
       write(*, "(I10,2ES12.2)") n, abs_res_norm, rel_res_norm
    end if

    ! Print warning message if not converged
    if (.not. conv) then
       write(*,'(A5, 2A12)') "ITER", "|R|", "|R|/|R1|"
       write(*, "(I5,2ES12.2)") n, abs_res_norm, rel_res_norm
       stop "Newton Solve Failed"
    end if

    if (allocated(res))    deallocate(res)
    if (allocated(dq))     deallocate(dq)
    if (allocated(jac))    deallocate(jac)
    
  end subroutine newton_solve_condensed
  
  !===================================================================! 
  ! Routine that approximates the Jacobian based on finite differences
  ! [d{R}/d{q}] = alpha*[dR/dq] + beta*[dR/dqdot] + gamma*[dR/dqddot]
  !===================================================================!

  subroutine approximate_jacobian( system, jac, coeff, U, X )

    class(dynamics)                              :: system
    type(scalar) , intent(inout) :: jac(:,:)
    type(scalar) , intent(inout)      :: U(:,:)                  ! states
    type(scalar)   , intent(in)    :: X(:,:)
!@    type(scalar) , allocatable, dimension(:)     :: pstate           ! perturbed ates
    type(scalar) , allocatable, dimension(:)     :: R, Rtmp            ! original residual and perturbed residual

    ! Scalars
    type(scalar)                                 :: dh = 1.0d-6       ! finite-diff step size
    type(scalar) , intent(in)                    :: coeff(:)          ! linearization coefficients
    integer                                      :: m                 ! loop variables
    integer :: nvars
    
    if (system % get_differential_order() .gt. 2) then
       print *, 'error: using approximate jacobian for system order > 2'
       stop
    end if

    !  Zero the supplied jacobian matrix for safety (as we are
    !  computing everything newly here)
    jac = 0.0d0

    nvars = system % get_num_state_vars()
       
    ! Allocate required arrays
    !   allocate(pstate(nvars)); pstate = 0.0d0;
    allocate(R(nvars));      R = 0.0d0;
    allocate(Rtmp(nvars));   Rtmp = 0.0d0;

    ! Make a residual call with original variables
    R = 0.0d0
    call system % add_residual(R,  U, X)

    !-----------------------------------------------------------!
    ! Derivative of R WRT Q: dR/dQ
    !-----------------------------------------------------------!
    
    associate (pstate => U(1,:), alpha=> coeff(1))

      loop_vars: do m = 1, nvars

         ! Perturb the k-th variable
         pstate(m) = pstate(m) + dh

         ! Make a residual call with the perturbed variable
         rtmp = 0
         call system % add_residual(Rtmp, U, X)

         ! Unperturb (restore) the k-th variable
         pstate(m) =  pstate(m) - dh

         ! Approximate the jacobian with respect to the k-th variable
         jac(:,m) = jac(:,m) + alpha*(Rtmp-R)/dh

      end do loop_vars
      
    end associate

    !-----------------------------------------------------------!
    ! Derivative of R WRT QDOT: dR/dQDOT
    !-----------------------------------------------------------!

    associate (pstate => U(2,:), beta => coeff(2))

      do m = 1, nvars

         ! Perturb the k-th variable
         pstate(m) = pstate(m) + dh

         ! Make a residual call with the perturbed variable
         rtmp = 0
         call system % add_residual(Rtmp, U, X)

         ! Unperturb (restore) the k-th variable
         pstate(m) = pstate(m) - dh

         ! Approximate the jacobian with respect to the k-th variable
         Jac(:,m) = Jac(:,m) + beta*(Rtmp-R)/dh

      end do

    end associate

    ! Second order equations have an extra block to add
    if (system % get_differential_order() == 2) then

       !-----------------------------------------------------------!
       ! Derivative of R WRT QDDOT: dR/dQDDOT
       !-----------------------------------------------------------!     
       associate (pstate => U(3,:), gamma => coeff(3))

       do m = 1, nvars

          ! Perturb the k-th variable
          pstate(m) = pstate(m) + dh

          ! Make a residual call with the perturbed variable
          rtmp = 0
          call system % add_residual(Rtmp, U, X)

          ! Unperturb (restore) the k-th variable
          pstate(m) = pstate(m) - dh

          ! Approximate the jacobian with respect to the k-th variable
          Jac(:,m) = Jac(:,m) + gamma*(Rtmp-R)/dh

       end do
       
     end associate

    end if ! first or second order

!    deallocate(pstate)
    deallocate(R,Rtmp)

  end subroutine approximate_jacobian

end module nonlinear_algebra
