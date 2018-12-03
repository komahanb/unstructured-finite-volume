!=====================================================================!
! Conjugate Gradient linear solver class that uses the functionalities of assembler class 
! in the iterative solution process.
!=====================================================================!

module class_sor

  use iso_fortran_env         , only : dp => REAL64
  use interface_linear_solver , only : linear_solver
  use class_assembler         , only : assembler

  implicit none
  
  ! Expose only the linear solver datatype
  private
  public :: sor

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!

  type, extends(linear_solver) :: sor

     !type(assembler), pointer :: FVAssembler
     class(assembler), allocatable :: FVAssembler
     integer                       :: print_level
     real(dp)                      :: omega

   contains

     ! type bound procedures
     procedure :: solve
     procedure :: iterate
     procedure :: estimate_spectral_radius

     ! destructor
     final :: destroy

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
  
  type(sor) function construct(FVAssembler, omega, max_it, &
       & max_tol, print_level) result (this)

    type(assembler), intent(in) :: FVAssembler
    type(real(dp)) , intent(in) :: omega
    type(integer)  , intent(in) :: max_it
    type(real(dp)) , intent(in) :: max_tol
    type(integer)  , intent(in) :: print_level

    allocate(this % FVassembler, source = FVAssembler)
    this % max_it      = max_it
    this % max_tol     = max_tol
    this % print_level = print_level
    this % omega       = omega

  end function construct

  !===================================================================!
  ! Destructor for linear solver
  !===================================================================!

  pure subroutine destroy(this)

    type(sor), intent(inout) :: this

!!$    if(associated(this % FVAssembler)) then
!!$       deallocate(this % FVAssembler)
!!$       nullify(this % FVAssembler)
!!$    end if
!!$    
    if (allocated(this % FVAssembler)) deallocate(this % FVAssembler)

  end subroutine destroy
  
  !===================================================================!
  ! Estimate spectral radius using power iteration (maximum absolute
  ! eigen value)
  !===================================================================!
  
  subroutine estimate_spectral_radius(this, mu, max_iter)
    
    ! Arguments
    class(sor), intent(in)  :: this
    real(dp)  , intent(out) :: mu
    integer   , intent(in)  :: max_iter

    ! Random vector
    real(dp) , allocatable :: v(:), w(:)
    integer  :: iter
    real(dp) :: wnorm

    ! Create a random unit vector
    call this % FVAssembler % create_vector(v)
    call random_number(v)
    v = v/norm2(v)

    ! A temp vector for processing
    call this % FVAssembler % create_vector(w)
    
    power_iteration: do iter = 1, max_iter       
       call this % FVAssembler % get_jacobian_vector_product(w, v)
       wnorm = norm2(w)
       v = w/wnorm
       mu = dot_product(v,w)
       print *, iter, mu
    end do power_iteration
    
    deallocate(v, w)

  end subroutine estimate_spectral_radius

  !===================================================================!
  ! Iterative linear solution using conjugate gradient method
  !===================================================================!
  
  subroutine solve(this, x)

    class(sor), intent(in)  :: this
    real(dp), allocatable    , intent(out) :: x(:)

    ! Locals
    real(dp), allocatable :: xold(:), ss(:)
    real(dp) :: update_norm
    integer  :: iter, num_inner_iters

    ! Initial guess vector for the subspace is "b"
    allocate(x(this % FVAssembler % num_state_vars))
    call this % FVAssembler % get_source(x)
    if (norm2(x) .lt. epsilon(1.0_dp)) error stop

    allocate(xold(this % FVAssembler % num_state_vars))
    allocate(ss(this % FVAssembler % num_state_vars))

    if (this % print_level .eq. -1) then
       open(13, file='sor.res', action='write')
       write(13,*) "iteration ", " residual"
    end if

    update_norm = huge(1.0d0);
    iter = 1;
    outer_iterations: do while (update_norm .gt. this % max_tol)

       xold = x

       ! Inner iterations with CG
       if ( iter .eq. 1) then
          ss = 0.0d0
          call this % iterate(x, ss, num_inner_iters)
       else
          call this % FVAssembler % get_skew_source(ss, x)
          call this % iterate(x, ss, num_inner_iters)
       end if
       
       update_norm = norm2(x - xold)
       if (this % print_level .eq. -1) then
          write(13, *) iter, update_norm
       end if
       if (this % print_level .gt. 0) then
          print *, "outer iter", iter, "num_inner_iters", num_inner_iters,  "update norm", update_norm
       end if
       iter = iter + 1

    end do outer_iterations

    if (this % print_level .eq. -1) then
       close(13)
    end if

    deallocate(xold)
    deallocate(ss)

  end subroutine solve
  
  !===================================================================!
  ! Iterative linear solution using conjugate gradient method
  !===================================================================!

  subroutine iterate(this, x, ss, iter)

    class(sor) , intent(in)    :: this
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(in)    :: ss(:)
    integer                   , intent(out)   :: iter

    ! Locals
    real(dp) :: tol
    real(dp) , allocatable :: b(:)
    real(dp) , allocatable :: Ux(:)
    real(dp) , allocatable :: D(:)
    real(dp) , allocatable :: Dx(:)
    real(dp) , allocatable :: R(:)
    real(dp) , allocatable :: xnew(:)
    real(dp) , allocatable :: identity(:)

    real(dp) :: bnorm

    real(dp) , allocatable :: Ly(:)
    real(dp) , allocatable :: y(:), ynew(:)

    allocate(b, Ux, D, Dx, R, xnew, identity, mold=x)
    allocate(y, ynew, Ly, mold=x)

    ! Identity vector
    identity = 1.0d0
    
    ! Extract the diagonal entries
    call this % FVAssembler % get_jacobian_vector_product(&
         & D, identity, filter = this % FVAssembler % DIAGONAL)

    ! Assemble RHS of the linear system (source + boundary terms)
    call this % FVAssembler % get_source(b)

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

       ! Form Ux
       call this % FVAssembler % get_jacobian_vector_product(&
            & Ux, x, filter = this % FVAssembler % UPPER_TRIANGLE)

       ! Form Dx
       call this % FVAssembler % get_jacobian_vector_product(&
            & Dx, x, filter = this % FVAssembler % DIAGONAL)

       ! R = w(b-Ux_k)+(1-w)Dx_k
       R = this % omega * (b - Ux) + (1.0_dp-this % omega)*Dx
      
       !--------------------------------------------------------------!
       ! Solve the linear system: By=R ; (D+wL)y=R ;  Dy=R-Ly
       !--------------------------------------------------------------!
       
       solve_lower_triangle: block

         real(dp) :: tol2
         integer :: iter2

         ! Initial guess is the current solution
         y = x

         iter2 = 1; tol2 = huge(1.0d0)
         do while ((tol2 .gt. this % max_tol) .and. (iter2 .lt. this % max_it))

            call this % FVAssembler % get_jacobian_vector_product(&
                 & Ly, y, filter = this % FVAssembler % LOWER_TRIANGLE)

            ynew  = (R - this % omega*Ly)/D
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
    
    deallocate(b, Ux, D, Dx, R, xnew, identity)
    deallocate(y, ynew, Ly)

  end subroutine iterate

end module class_sor

!!$    real(dp) :: rho
!!$    call this % estimate_spectral_radius(rho, 100)
!!$    print *, rho
!!$    stop
