!=====================================================================!
! Conjugate Gradient linear solver class that uses the functionalities of assembler class 
! in the iterative solution process.
!=====================================================================!

module class_gauss_seidel

  use iso_fortran_env         , only : dp => REAL64
  use interface_linear_solver , only : linear_solver
  use interface_assembler     , only : assembler, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE

  implicit none
  
  ! Expose only the linear solver datatype
  private
  public :: gauss_seidel

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!

  type, extends(linear_solver) :: gauss_seidel

     ! stateless w.r.t. the system: solve takes the assembler as an argument
     integer                       :: print_level

   contains

     ! type bound procedures
     procedure :: solve
     procedure :: iterate

     ! destructor
     final :: destroy

  end type gauss_seidel

  !===================================================================!
  ! Interface for multiple constructors
  !===================================================================!

  interface gauss_seidel
     module procedure construct
  end interface gauss_seidel

contains
  
  !===================================================================!
  ! Constructor for linear solver
  !===================================================================!
  
  pure type(gauss_seidel) function construct(max_it, &
       & max_tol, print_level) result (this)

    type(integer)  , intent(in) :: max_it
    type(real(dp)) , intent(in) :: max_tol
    type(integer)  , intent(in) :: print_level

    this % max_it      = max_it
    this % max_tol     = max_tol
    this % print_level = print_level

  end function construct

  !===================================================================!
  ! Destructor for linear solver
  !===================================================================!

  pure subroutine destroy(this)

    type(gauss_seidel), intent(inout) :: this

!!$    end if
!!$    

  end subroutine destroy

  !===================================================================!
  ! Iterative linear solution using conjugate gradient method
  !===================================================================!
  
  impure subroutine solve(this, system, x, mode)

    class(gauss_seidel), intent(in)           :: this
    class(assembler)   , intent(in)           :: system
    real(dp), allocatable    , intent(out) :: x(:)
    integer              , intent(in), optional :: mode  ! FORWARD (default) / REVERSE

    ! Locals
    real(dp), allocatable :: xold(:), ss(:)
    real(dp) :: update_norm, rnorm0
    integer  :: iter, num_inner_iters

    ! Initial guess vector for the subspace is "b"
    allocate(x(system % num_state_vars))
    call system % get_source(x)
    if (norm2(x) .lt. epsilon(1.0_dp)) then
       print *, 'zero rhs? stopping'
       error stop
    end if

    allocate(xold(system % num_state_vars))
    allocate(ss(system % num_state_vars))

    if (this % print_level .eq. -1) then
       open(13, file='gs.res', action='write')
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

    class(gauss_seidel) , intent(in)    :: this
    class(assembler)    , intent(in)    :: system
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(in)    :: ss(:)
    integer                   , intent(out)   :: iter

    ! Locals
    real(dp) :: tol
    real(dp) , allocatable :: b(:)
    real(dp) , allocatable :: Ux(:)
    real(dp) , allocatable :: D(:)
    real(dp) , allocatable :: R(:)
    real(dp) , allocatable :: xnew(:)
    real(dp) , allocatable :: identity(:)

    real(dp) :: bnorm

    real(dp) , allocatable :: Ly(:)
    real(dp) , allocatable :: y(:), ynew(:)

    allocate(b, Ux, D, R, xnew, identity, mold=x)
    allocate(y, ynew, Ly, mold=x)

    ! Identity vector
    identity = 1.0d0
    
    ! Extract the diagonal entries
    call system % get_jacobian_vector_product(&
         & D, identity, filter = DIAGONAL)

    ! Assemble RHS of the linear system (source + boundary terms)
    call system % get_source(b)

    ! Add the skew source terms if supplied 
    b = b + ss

    bnorm = norm2(b)

    ! Homogeneous case (nothing to do)
    if (bnorm .le. this % max_tol) then
       x = 0.0d0
       return
    end if

    !-----------------------------------------------------------------!
    ! Apply Gauss Jacobi Iterative scheme until tolerance is achieved
    !-----------------------------------------------------------------!

    iter = 1; tol  = huge(1.0d0);
    do while ((tol .gt. this % max_tol) .and. (iter .lt. this % max_it))

       ! Form the residual (partial) after split
       call system % get_jacobian_vector_product(&
            & Ux, x, filter = UPPER_TRIANGLE)

       R = b - Ux

       ! call system % apply_preconditioner(x, D)
       ! Invert diagonal
       !xnew = R/D ! (D+L)^{-1}(b-Ux)       

       !--------------------------------------------------------------!
       ! Solve the linear system: By=R ; (D+L)y=R ;  Dy=R-Ly
       !--------------------------------------------------------------!
       
       solve_lower_triangle: block

         real(dp) :: tol2
         integer :: iter2

         ! Initial guess is the current solution
         y = x

         iter2 = 1; tol2 = huge(1.0d0)
         do while ((tol2 .gt. this % max_tol) .and. (iter2 .lt. this % max_it))

            call system % get_jacobian_vector_product(&
                 & Ly, y, filter = LOWER_TRIANGLE)

            ynew  = (R - Ly)/D
            tol2  = norm2(y-ynew)
            
            if (this % print_level .gt. 2) then
               write(*,*) "inner (2)", iter2, tol2
            end if
            
            y = ynew
            iter2 = iter2 + 1

         end do
         
       end block solve_lower_triangle
       
       xnew = y 
       tol = norm2(x-xnew)
       
       if (this % print_level .gt. 1) then
          write(*,*) "inner (1)", iter, tol
       end if

       x = xnew
       iter = iter + 1

    end do
    
    deallocate(b, Ux, D, R, xnew, identity)
    deallocate(y, ynew, Ly)

  end subroutine iterate

end module class_gauss_seidel
