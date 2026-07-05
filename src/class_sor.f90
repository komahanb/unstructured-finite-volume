!=====================================================================!
! Successive over-relaxation (SOR) linear solver: supplies only the
! correction sweep (`iterate`); the deferred-correction outer solve is
! inherited from linear_solver.
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

     ! stateless w.r.t. the system: solve takes the assembler as an argument
     real(dp)                      :: omega

   contains

     ! correction sweep consumed by the inherited deferred-correction solve
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
  ! Destructor for linear solver
  !===================================================================!

  pure subroutine destroy(this)

    type(sor), intent(inout) :: this

!!$    end if
!!$    

  end subroutine destroy
  
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
    v = v/norm2(v)

    ! A temp vector for processing
    call system % create_vector(w)
    
    power_iteration: do iter = 1, max_iter       
       call system % get_jacobian_vector_product(w, v)
       wnorm = norm2(w)
       v = w/wnorm
       mu = dot_product(v,w)
       print *, iter, mu
    end do power_iteration
    
    deallocate(v, w)

  end subroutine estimate_spectral_radius

  !===================================================================!
  ! Correction sweep: SOR on the skew-corrected system. The relaxed
  ! lower-triangle solve (D+wL)y = R is itself done iteratively.
  !===================================================================!

  impure subroutine iterate(this, system, x, ss, iter)

    class(sor)                , intent(in)    :: this
    class(assembler)          , intent(in)    :: system
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
    ! Apply the SOR sweep until tolerance is achieved
    !-----------------------------------------------------------------!

    iter = 1; tol  = huge(1.0d0);
    do while ((tol .gt. this % max_tol) .and. (iter .lt. this % max_it))

       ! Form Ux
       call system % get_jacobian_vector_product(&
            & Ux, x, filter = UPPER_TRIANGLE)

       ! Form Dx
       call system % get_jacobian_vector_product(&
            & Dx, x, filter = DIAGONAL)

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

            call system % get_jacobian_vector_product(&
                 & Ly, y, filter = LOWER_TRIANGLE)

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
