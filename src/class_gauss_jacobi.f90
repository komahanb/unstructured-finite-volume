!=====================================================================!
! Gauss-Jacobi linear solver: supplies only the correction sweep
! (`iterate`); the deferred-correction outer solve is inherited from
! linear_solver.
!=====================================================================!

module class_gauss_jacobi

  use iso_fortran_env         , only : dp => REAL64
  use interface_linear_solver , only : linear_solver
  use interface_assembler     , only : assembler, DIAGONAL, LOWER_TRIANGLE, UPPER_TRIANGLE

  implicit none
  
  ! Expose only the linear solver datatype
  private
  public :: gauss_jacobi

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!

  type, extends(linear_solver) :: gauss_jacobi

     ! stateless w.r.t. the system: solve takes the assembler as an argument

   contains

     ! correction sweep consumed by the inherited deferred-correction solve
     procedure :: iterate

     ! destructor
     final :: destroy

  end type gauss_jacobi

  !===================================================================!
  ! Interface for multiple constructors
  !===================================================================!

  interface gauss_jacobi
     module procedure construct
  end interface gauss_jacobi

contains
  
  !===================================================================!
  ! Constructor for linear solver
  !===================================================================!
  
  pure type(gauss_jacobi) function construct(max_it, &
       & max_tol, print_level) result (this)

    type(integer)  , intent(in) :: max_it
    type(real(dp)) , intent(in) :: max_tol
    type(integer)  , intent(in) :: print_level

    this % max_it      = max_it
    this % max_tol     = max_tol
    this % print_level = print_level
    this % res_file    = 'gj.res'

  end function construct

  !===================================================================!
  ! Destructor for linear solver
  !===================================================================!

  pure subroutine destroy(this)

    type(gauss_jacobi), intent(inout) :: this

!!$    end if
!!$    

  end subroutine destroy

  !===================================================================!
  ! Correction sweep: Jacobi iteration on the skew-corrected system
  !===================================================================!

  impure subroutine iterate(this, system, x, ss, iter)

    class(gauss_jacobi) , intent(in)    :: this
    class(assembler)    , intent(in)    :: system
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(in)    :: ss(:)
    integer                   , intent(out)   :: iter

    ! Locals
    real(dp) :: tol
    real(dp) , allocatable :: b(:)
    real(dp) , allocatable :: Ux(:)
    real(dp) , allocatable :: Lx(:)
    real(dp) , allocatable :: D(:)
    real(dp) , allocatable :: R(:)
    real(dp) , allocatable :: xnew(:)
    real(dp) , allocatable :: identity(:)

    real(dp) :: bnorm

    allocate(b, Ux, Lx, D, R, xnew, identity, mold=x)

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
       
       call system % get_jacobian_vector_product(&
            & Lx, x, filter = LOWER_TRIANGLE)
       
       R = b - Lx - Ux

       ! call system % apply_preconditioner(x, D)
       ! Invert diagonal
       xnew = R/D ! D^{-1}(b-Lx-Ux)
       tol = norm2(x-xnew)

       if (this % print_level .gt. 1) then
          write(*,*) "inner (1)", iter, tol
       end if

       x = xnew
       iter = iter + 1

    end do
    
    deallocate(b, Ux, Lx, D, R, xnew, identity)

  end subroutine iterate

end module class_gauss_jacobi
