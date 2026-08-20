!=====================================================================!
! Newton's iteration: a minimizer like every other.
!
! One family, one story: attach a statement, drive its residual to
! zero. The linear members reach the answer in the statement's own
! space; newton reaches it by linearizing where it stands,
!
!      J(q) dq = rhs - action(q)          q <- q + dq
!
! and it is not a different kind of thing for that - it extends the
! same base, wears the same operation face, and answers the same
! solve(rhs, x, achieved). A network's training loop would join the
! family the same way: another concretion of the one creed.
!
! Newton owns no derivative mathematics. The tangent is a level-1
! citizen - the linearization operator - and newton merely governs:
! freeze the linearization at the standing state, hand the linear
! question to the inner minimizer, step. The seat is filled by what
! the statement IS: a differentiable statement linearizes itself
! exactly, anything else is differenced - the promise the family
! made, kept by one dispatch, and the governance below it never
! changes.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_newton

  use iso_fortran_env           , only : dp => REAL64
  use operation_minimization        , only : minimizer
  use operation_discretization            , only : linearization_operator
  use operation_linearization, only : tangent_of

  implicit none

  private
  public :: newton

  !===================================================================!
  ! Newton: one component beyond the family - the minimizer it
  ! governs, one rank down, handed one linear question per step.
  !===================================================================!

  type, extends(minimizer) :: newton

     class(minimizer), allocatable :: inner

   contains

     procedure :: name => newton_name
     procedure :: solve

  end type newton

contains

  pure function newton_name(this) result(name)

    class(newton), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'newton'

  end function newton_name

  !===================================================================!
  ! Drive action(q) toward rhs from the given q.
  !===================================================================!

  subroutine solve(this, rhs, x, achieved)

    class(newton), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    class(linearization_operator), allocatable :: jacobian
    real(dp), allocatable :: residual(:), g(:), y(:), dq(:)
    real(dp) :: linear_achieved
    integer :: it

    allocate(dq(size(x)))

    call this % constant(g)

    ! the seat is filled by what the statement is - the shelf's
    ! own chooser answers, and no dispatch lives here
    jacobian = tangent_of(this % action)

    do it = 1, this % max_iterations

       ! Where we stand: the full statement, whatever its linearity.
       call this % matvec(x, y)
       residual = y + g - rhs

       achieved = this % norm(residual)
       if (achieved < this % tolerance) return

       ! The linear question at this point, answered by the governed
       ! minimizer.
       call jacobian % freeze(x, base=y + g)

       call this % inner % attach(jacobian, this % on, this % unknown_domain, &
            & this % n_unknown_domain, ncomp = this % ncomp)
       dq = 0.0_dp
       call this % inner % solve(-residual, dq, linear_achieved)

       x = x + dq

    end do

    call this % matvec(x, y)
    achieved = this % norm(y + g - rhs)

  end subroutine solve

end module operation_newton
