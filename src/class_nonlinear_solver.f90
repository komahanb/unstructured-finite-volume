!=====================================================================!
! Abstract class for nonlinear solvers to extend and provide iterate
! method.
! 
! Author : Komahan Boopathy
!=====================================================================!

module class_nonlinear_solver

  use iso_fortran_env           , only : dp => REAL64
  use interface_algebraic_solver, only : algebraic_solver

  implicit none
  
  ! Expose only the linear solver interface
  private
  public :: nonlinear_solver

  !===================================================================!
  ! Abstract nonlinear solver datatype
  !===================================================================!
  
  type, abstract, extends(algebraic_solver) :: nonlinear_solver

     class(algebraic_solver), allocatable :: linear_solver

     ! Nonlinear solvers need a linear solver to solve a linearized
     ! nonlinear system, and they implement the solve method required
     ! by algebraic solver interface
     integer  :: print_level
     real(dp) :: max_tol
     integer  :: max_it
     
   contains
     
     ! Generic solve procedure
     ! procedure :: solve
     
  end type nonlinear_solver
  
contains
  
  !===================================================================!
  ! Constructor for nonlinear solver
  !===================================================================!
  
  type(nonlinear_solver) function construct(FVAssembler, &
       & max_it, max_tol, print_level) result (this)
    
    type(assembler), intent(in) :: FVAssembler
    type(real(dp)) , intent(in) :: omega
    type(integer)  , intent(in) :: max_it
    type(real(dp)) , intent(in) :: max_tol
    type(integer)  , intent(in) :: print_level
    
    allocate(this % FVassembler, source = FVAssembler)
    this % max_it      = max_it
    this % max_tol     = max_tol
    this % print_level = print_level
    
  end function construct

  !===================================================================!
  ! Destructor for nonlinear solver
  !===================================================================!
  
  pure subroutine destroy(this)

    type(nonlinear_solver), intent(inout) :: this
    
!!$    if(associated(this % FVAssembler)) then
!!$       deallocate(this % FVAssembler)
!!$       nullify(this % FVAssembler)
!!$    end if
!!$    
    if (allocated(this % FVAssembler)) deallocate(this % FVAssembler)
    
  end subroutine destroy
  
end module class_nonlinear_solver
