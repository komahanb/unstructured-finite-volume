!=====================================================================!
! Diagonal-sweep linear solver (traditionally: Gauss-Jacobi): supplies
! only the sweep (`iterate`) - inverting the self-loop subgraph of the
! operator - and inherits the residual-minimization iteration from
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

   contains

     ! the sweep consumed by the inherited outer iteration
     procedure :: iterate

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
  ! The sweep: jacobi iteration on the correction equation A dx = r
  ! from dx = 0, inverting the diagonal (the self-loops) each pass.
  !===================================================================!

  impure subroutine iterate(this, system, r, dx, iter)

    class(gauss_jacobi) , intent(inout) :: this
    class(assembler)    , intent(in)  :: system
    real(dp)            , intent(in)  :: r(:)
    real(dp)            , intent(out) :: dx(:)
    integer             , intent(out) :: iter

    ! Locals
    real(dp) :: tol, bnorm
    real(dp) , allocatable :: Ux(:), Lx(:), D(:), R2(:), xnew(:), identity(:)

    dx = 0.0_dp
    allocate(Ux, Lx, D, R2, xnew, identity, mold = dx)

    ! Extract the diagonal entries (the self-loop subgraph on ones)
    identity = 1.0d0
    call system % get_jacobian_residual_product(D, identity, part = DIAGONAL)

    bnorm = sqrt(system % inner_product(r, r))

    ! Homogeneous case (nothing to do)
    if (bnorm .le. this % max_tol) then
       iter = 0
       return
    end if

    !-----------------------------------------------------------------!
    ! Apply the jacobi sweep until tolerance is achieved
    !-----------------------------------------------------------------!

    iter = 1; tol  = huge(1.0d0)
    do while ((tol .gt. this % max_tol) .and. (iter .lt. this % max_it))

       ! Form the residual (partial) after the split
       call system % get_jacobian_residual_product(Ux, dx, part = UPPER_TRIANGLE)
       call system % get_jacobian_residual_product(Lx, dx, part = LOWER_TRIANGLE)

       R2 = r - Lx - Ux

       ! Invert the diagonal
       xnew = R2/D ! D^{-1}(r - L dx - U dx)
       tol  = sqrt(system % inner_product(dx - xnew, dx - xnew))

       if (this % print_level .gt. 1) then
          write(*,*) "inner (1)", iter, tol
       end if

       dx   = xnew
       iter = iter + 1

    end do

    deallocate(Ux, Lx, D, R2, xnew, identity)

  end subroutine iterate

end module class_gauss_jacobi
