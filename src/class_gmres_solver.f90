!=====================================================================!
! linear_solver wrapper around restarted GMRES, so the config-driven
! driver can pick "gmres" the same way it picks "cg". GMRES needs the
! assembled operator (it is not matrix-free here), so solve() assembles the
! csr once via get_operator_csr, takes the rhs from get_source, and runs
! class_gmres % gmres. Use it for nonsymmetric operators (advection) where
! CG does not apply.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_gmres_solver

  use iso_fortran_env         , only : dp => REAL64
  use interface_linear_solver , only : linear_solver
  use class_assembler         , only : assembler
  use class_csr               , only : csr_matrix
  use class_gmres             , only : gmres

  implicit none

  private
  public :: gmres_solver

  type, extends(linear_solver) :: gmres_solver
     class(assembler), allocatable :: FVAssembler
     integer                       :: restart     = 200
     integer                       :: print_level = 0
   contains
     procedure :: solve => gmres_solve
  end type gmres_solver

  interface gmres_solver
     module procedure construct
  end interface gmres_solver

contains

  type(gmres_solver) function construct(FVAssembler, max_it, max_tol, restart, print_level) result(this)
    type(assembler), intent(in)           :: FVAssembler
    integer        , intent(in)           :: max_it
    real(dp)       , intent(in)           :: max_tol
    integer        , intent(in), optional :: restart, print_level
    allocate(this % FVAssembler, source = FVAssembler)
    this % max_it  = max_it
    this % max_tol = max_tol
    if (present(restart))     this % restart     = restart
    if (present(print_level)) this % print_level = print_level
  end function construct

  subroutine gmres_solve(this, x)
    class(gmres_solver)  , intent(in)  :: this
    real(dp), allocatable, intent(out) :: x(:)
    type(csr_matrix)      :: A
    real(dp), allocatable :: b(:)
    integer :: n
    call this % FVAssembler % get_operator_csr(A)
    n = this % FVAssembler % num_state_vars
    allocate(b(n), x(n))
    call this % FVAssembler % get_source(b)
    call gmres(A, b, x, this % max_it, this % restart, this % max_tol, this % print_level)
  end subroutine gmres_solve

end module class_gmres_solver
