!=====================================================================!
! Gauss-Seidel linear solver: supplies only the correction sweep
! (`iterate`); the deferred-correction outer solve is inherited from
! linear_solver.
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

   contains

     ! correction sweep consumed by the inherited deferred-correction solve
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
    this % res_file    = 'gs.res'

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
  ! Correction sweep: Gauss-Seidel on the skew-corrected system. The
  ! lower-triangle solve (D+L)y = R is itself done iteratively.
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
    ! Apply the Gauss-Seidel sweep until tolerance is achieved
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
