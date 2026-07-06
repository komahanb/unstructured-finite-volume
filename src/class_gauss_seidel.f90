!=====================================================================!
! Forward-triangle-sweep linear solver (traditionally: Gauss-Seidel):
! supplies only the sweep (`iterate`) - the lower-triangle solve
! (D+L)y = R is itself done iteratively - and inherits the
! residual-minimization iteration from linear_solver.
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

   contains

     ! the sweep consumed by the inherited outer iteration
     procedure :: iterate

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
  ! The sweep: gauss-seidel on the correction equation A dx = r from
  ! dx = 0. The lower-triangle solve (D+L)y = R is itself iterative.
  !===================================================================!

  impure subroutine iterate(this, system, r, dx, iter)

    class(gauss_seidel) , intent(in)  :: this
    class(assembler)    , intent(in)  :: system
    real(dp)            , intent(in)  :: r(:)
    real(dp)            , intent(out) :: dx(:)
    integer             , intent(out) :: iter

    ! Locals
    real(dp) :: tol, bnorm
    real(dp) , allocatable :: Ux(:), D(:), R2(:), xnew(:), identity(:)
    real(dp) , allocatable :: Ly(:), y(:), ynew(:)

    dx = 0.0_dp
    allocate(Ux, D, R2, xnew, identity, mold = dx)
    allocate(y, ynew, Ly, mold = dx)

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
    ! Apply the gauss-seidel sweep until tolerance is achieved
    !-----------------------------------------------------------------!

    iter = 1; tol  = huge(1.0d0)
    do while ((tol .gt. this % max_tol) .and. (iter .lt. this % max_it))

       ! Form the residual (partial) after the split
       call system % get_jacobian_residual_product(Ux, dx, part = UPPER_TRIANGLE)

       R2 = r - Ux

       !--------------------------------------------------------------!
       ! Solve the linear system: By=R ; (D+L)y=R ;  Dy=R-Ly
       !--------------------------------------------------------------!

       solve_lower_triangle: block

         real(dp) :: tol2
         integer :: iter2

         ! Initial guess is the current correction
         y = dx

         iter2 = 1; tol2 = huge(1.0d0)
         do while ((tol2 .gt. this % max_tol) .and. (iter2 .lt. this % max_it))

            call system % get_jacobian_residual_product(Ly, y, part = LOWER_TRIANGLE)

            ynew  = (R2 - Ly)/D
            tol2  = sqrt(system % inner_product(y - ynew, y - ynew))

            if (this % print_level .gt. 2) then
               write(*,*) "inner (2)", iter2, tol2
            end if

            y = ynew
            iter2 = iter2 + 1

         end do

       end block solve_lower_triangle

       xnew = y
       tol  = sqrt(system % inner_product(dx - xnew, dx - xnew))

       if (this % print_level .gt. 1) then
          write(*,*) "inner (1)", iter, tol
       end if

       dx   = xnew
       iter = iter + 1

    end do

    deallocate(Ux, D, R2, xnew, identity)
    deallocate(y, ynew, Ly)

  end subroutine iterate

end module class_gauss_seidel
