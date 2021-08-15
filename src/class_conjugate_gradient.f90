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

     !type(assembler), pointer :: FVAssembler
     class(assembler), allocatable :: FVAssembler
     integer                       :: print_level

   contains

     ! type bound procedures
     procedure :: solve
     procedure :: iterate

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

  type(conjugate_gradient) function construct(FVAssembler, max_it, &
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

    type(conjugate_gradient), intent(inout) :: this

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

    class(conjugate_gradient), intent(in)  :: this
    real(dp), allocatable    , intent(out) :: x(:)

    ! Locals
    real(dp), allocatable :: xold(:), ss(:)
    real(dp) :: update_norm
    integer  :: iter, num_inner_iters

    ! Initial guess vector for the subspace is "b"
    allocate(x(this % FVAssembler % num_state_vars))
    call this % FVAssembler % get_source(x)
    if (norm2(x) .lt. epsilon(1.0_dp)) then
       print *, 'zero rhs?'
       ! error stop
    end if

    allocate(xold(this % FVAssembler % num_state_vars))
    allocate(ss(this % FVAssembler % num_state_vars))

    if (this % print_level .eq. -1) then
       open(13, file='cg.res', action='write')
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

    class(conjugate_gradient) , intent(in)    :: this
    real(dp)                  , intent(inout) :: x(:)
    real(dp)                  , intent(in)    :: ss(:)
    integer                   , intent(out)   :: iter

    ! Create local data
    real(dp), allocatable :: p(:), r(:), w(:), Ax(:), tmp(:)
    real(dp), allocatable :: b(:)
    real(dp)              :: alpha, beta
    real(dp)              :: bnorm, rnorm
    real(dp)              :: tol
    real(dp)              :: rho(2)


    ! Start the iteration counter
    iter = 1

    ! Memory allocations
    allocate(b,p,r,w,Ax,tmp,mold=x)

    ! Norm of the right hand side
    call this % FVAssembler % get_source(tmp)
    ! Add the additional right hand side supplied
    tmp = tmp + ss
    if (this % FVAssembler % symmetry .eqv. .false.) then
       call this % FVAssembler % get_transpose_jacobian_vector_product(b, tmp)
    else
       b = tmp
    end if
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

    !open(13, file='cg.log', action='write', position='append')

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

       ! write(13,*) iter, tol
       if (this % print_level .gt. 1) then
          write(*,*) iter, tol, rnorm, rho ! causes valgrind errors
       end if

       iter = iter + 1

       rho(1) = rho(2)
       rho(2) = rnorm*rnorm

    end do

    !close(13)

    deallocate(r, p, w, b, Ax, tmp)

  end subroutine iterate

end module class_conjugate_gradient
