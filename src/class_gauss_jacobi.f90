!=====================================================================!
! Conjugate Gradient linear solver class that uses the functionalities of assembler class 
! in the iterative solution process.
!=====================================================================!

module class_gauss_jacobi

  use iso_fortran_env         , only : dp => REAL64
  use interface_linear_solver , only : linear_solver
  use class_assembler         , only : assembler

  implicit none
  
  ! Expose only the linear solver datatype
  private
  public :: gauss_jacobi

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!

  type, extends(linear_solver) :: gauss_jacobi

     !type(assembler), pointer :: FVAssembler
     class(assembler), allocatable :: FVAssembler
     integer                       :: print_level

   contains

     ! type bound procedures
     procedure :: solve
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
  
  type(gauss_jacobi) function construct(FVAssembler, max_it, &
       & max_tol, print_level) result (this)

    type(assembler), intent(in) :: FVAssembler
    type(integer)  , intent(in) :: max_it
    type(real(dp)) , intent(in) :: max_tol
    type(integer)  , intent(in) :: print_level

    allocate(this % FVassembler, source = FVAssembler)
    this % max_it      = max_it
    this % max_tol     = max_tol
    this % print_level = print_level

  end function construct

  !===================================================================!
  ! Destructor for linear solver
  !===================================================================!

  pure subroutine destroy(this)

    type(gauss_jacobi), intent(inout) :: this

!!$    if(associated(this % FVAssembler)) then
!!$       deallocate(this % FVAssembler)
!!$       nullify(this % FVAssembler)
!!$    end if
!!$    
    if (allocated(this % FVAssembler)) deallocate(this % FVAssembler)

  end subroutine destroy

  !===================================================================!
  ! Iterative linear solution using conjugate gradient method
  !===================================================================!
  
  subroutine solve(this, x)

    class(gauss_jacobi), intent(in)  :: this
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
       if (this % print_level .gt. 0) then
          print *, "outer iter", iter, "num_inner_iters", num_inner_iters,  "update norm", update_norm
       end if
       iter = iter + 1

    end do outer_iterations

    deallocate(xold)
    deallocate(ss)

  end subroutine solve
  
  !===================================================================!
  ! Iterative linear solution using conjugate gradient method
  !===================================================================!

  subroutine iterate(this, x, ss, iter)

    class(gauss_jacobi) , intent(in)    :: this
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(in)    :: ss(:)
    integer                   , intent(out)   :: iter

    ! Locals
    real(dp) :: tol

    iter = 1
    tol  = huge(1.0d0)

    ! Apply Gauss Jacobi Iterative scheme until tolerance is achieved
    do while ((tol .gt. this % max_tol) .and. (iter .lt. this % max_it))

       print *, 'doing gauss jacobi iteration'

       return

    end do

  end subroutine iterate

end module class_gauss_jacobi
