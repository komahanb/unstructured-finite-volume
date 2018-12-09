!=====================================================================!
! Abstract class for linear solvers to extend and provide iterate
! method.
! 
! Author: Komahan Boopathy
!=====================================================================!

module interface_linear_solver

  use iso_fortran_env           , only : dp => REAL64
  use interface_algebraic_solver, only : algebraic_solver
  use interface_assembler       , only : assembler

  implicit none
  
  private
  public :: linear_solver

  !===================================================================!
  ! Abstract linear solver datatype
  !===================================================================!
  
  type, abstract, extends(algebraic_solver) :: linear_solver
     
     ! Tolerances and print control for iterative solvers
     real(dp) :: max_tol
     integer  :: max_it
     integer  :: print_level
     
   contains
     
     ! Generic solve procedure
     procedure :: solve

     ! Iteration behavior is deferred (specific to the type of solver)
     procedure(iterate_interface), deferred :: iterate
     
  end type linear_solver
  
  !===================================================================!
  ! Interface for deferred iterate function 
  !===================================================================!
  
  interface
     subroutine iterate_interface(this, system, x)
       import linear_solver
       import assembler
       import dp
       class(linear_solver)  , intent(in)  :: this
       class(assembler)      , intent(in)  :: system
       real(dp), allocatable , intent(out) :: x(:)
     end subroutine iterate_interface
  end interface

contains

  !===================================================================!
  ! Iterative nonlinear solution A x = b(x)
  !===================================================================!
  
  subroutine solve(this, system)

    class(linear_solver)   , intent(in)  :: this
    class(assembler)       , intent(in)  :: system

    ! Locals
    real(dp), allocatable :: xold(:), ss(:)
    real(dp) :: update_norm
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
       open(13, file='sor.res', action='write')
       write(13,*) "iteration ", " residual"
    end if

    update_norm = huge(1.0d0); iter = 1;
    outer_iterations: do while (update_norm .gt. this % max_tol)

       xold = x

       ! Inner iterations with CG (linear solver)
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
          print *, &
               & "outer iter      ", iter, &
               & "num_inner_iters ", num_inner_iters,  &
               & "update norm     ", update_norm
       end if
       iter = iter + 1
       
    end do outer_iterations

    if (this % print_level .eq. -1) then
       close(13)
    end if

    deallocate(xold)
    deallocate(ss)

  end subroutine solve
  
end module interface_linear_solver
