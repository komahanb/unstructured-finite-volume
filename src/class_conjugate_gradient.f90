!=====================================================================!
! Conjugate Gradient linear solver class that uses the functionalities of assembler class 
! in the iterative solution process.
!=====================================================================!

module class_conjugate_gradient

  use iso_fortran_env         , only : dp => REAL64
  use interface_linear_solver , only : linear_solver
  use class_assembler         , only : assembler

  implicit none
  
  ! Expose only the linear solver datatype
  private
  public :: conjugate_gradient

  !===================================================================!
  ! Linear solver datatype
  !===================================================================!

  type, extends(linear_solver) :: conjugate_gradient

     class(assembler), pointer :: FVAssembler

   contains

     ! type bound procedures
     procedure :: solve

     ! destructor
     final :: destroy

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
  
  type(conjugate_gradient) function construct(FVAssembler) result (this)

    class(assembler), intent(in) :: FVAssembler

    allocate(this % FVassembler, source = FVAssembler)  

  end function construct

  !===================================================================!
  ! Destructor for linear solver
  !===================================================================!

  pure subroutine destroy(this)

    type(conjugate_gradient), intent(inout) :: this

    if(associated(this % FVAssembler)) then
       deallocate(this % FVAssembler)
       nullify(this % FVAssembler)
    end if

  end subroutine destroy

  !===================================================================!
  ! Iterative linear solution using conjugate gradient method
  !===================================================================!

  subroutine solve(this, x, ss)

    class(conjugate_gradient) , intent(in)    :: this
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(in)    :: ss(:)

    ! Create local data
    real(dp), allocatable :: p(:), r(:), w(:), Ax(:), tmp(:)
    real(dp), allocatable :: b(:)
    real(dp)              :: alpha, beta
    real(dp)              :: bnorm, rnorm
    real(dp)              :: tol
    integer               :: iter
    real(dp)              :: rho(2)

    ! Start the iteration counter
    iter = 1

    ! Memory allocations
    allocate(b,p,r,w,Ax,tmp,mold=x)

    ! Norm of the right hand side
    call this % FVAssembler % get_source(tmp)
    if (this % FVAssembler % symmetry .eqv. .false.) then
       call this % FVAssembler % get_transpose_jacobian_vector_product(b, tmp)
    else
       b = tmp
    end if
    ! Add the additional right hand side supplied
    b = b + ss
    bnorm = norm2(b)

    ! Homogeneous case
    if (bnorm .le. this % max_tol) then
       x = 0.0d0
       return
    end if

    ! Norm of the initial residual
    call this % FVAssembler % get_jacobian_vector_product(tmp, x)
    if (this % FVAssembler % symmetry .eqv. .false.) then
       call this % FVAssembler % get_transpose_jacobian_vector_product(Ax, tmp)
    else
       Ax = tmp
    end if
    r         = b - Ax ! could directly form this residual using get_residual_call
    rnorm     = norm2(r)
    tol       = rnorm/bnorm
    rho(2)    = rnorm*rnorm

    open(13, file='cg.log', action='write', position='append')

    ! Apply Iterative scheme until tolerance is achieved
    do while ((tol .gt. this % max_tol) .and. (iter .lt. this % max_it))

       ! step (a) compute the descent direction
       if ( iter .eq. 1) then
          ! steepest descent direction p
          p = r
       else
          ! take a conjugate direction
          beta = rho(2)/rho(1)
          p = r + beta*p
       end if

       ! step (b) compute the solution update
       call this % FVAssembler % get_jacobian_vector_product(tmp, p)
       if (this % FVAssembler % symmetry .eqv. .false.) then
          call this % FVAssembler % get_transpose_jacobian_vector_product(w, tmp)
       else
          w = tmp
       end if

       ! step (c) compute the step size for update
       alpha = rho(2)/dot_product(p, w)

       ! step (d) Add dx to the old solution
       x = x + alpha*p

       ! step (e) compute the new residual
       r = r - alpha*w

       ! step(f) update values before next iteration
       rnorm = norm2(r)
       tol = rnorm/bnorm

       write(13,*) iter, tol
       write(*,*) iter, tol, rnorm, rho ! causes valgrind errors

       iter = iter + 1

       rho(1) = rho(2)
       rho(2) = rnorm*rnorm

    end do

    close(13)

    deallocate(r, p, w, b, Ax, tmp)

  end subroutine solve

end module class_conjugate_gradient