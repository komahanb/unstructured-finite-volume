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
  use field_stored  , only : stored_field
  use operation_linearization, only : linearization, tangent_of

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

    type(linearization) :: jacobian
    type(stored_field), allocatable :: inputs(:)
    real(dp), allocatable :: residual(:), g(:), y(:), dq(:)
    real(dp) :: linear_achieved
    integer :: it

    allocate(dq(size(x)))

    call this % constant(g)

    ! the tangent in the unknown's argument; which road it takes is
    ! the statement's own answer, and no dispatch lives here
    jacobian = tangent_of(this % action, this % action % argument(1))

    do it = 1, this % max_iterations

       ! Where we stand: the full statement, whatever its linearity.
       call this % matvec(x, y)
       residual = y + g - rhs

       achieved = this % norm(residual)
       if (achieved < this % tolerance) return

       ! The linear question at this point, answered by the governed
       ! minimizer: the Jacobian is frozen at the same input tuple the
       ! residual was evaluated on, held inputs included.
       call this % evaluation_inputs(x, inputs)
       call jacobian % freeze(inputs, base=y + g)

       call this % inner % attach(jacobian, this % on, this % unknown_domain, &
            & this % num_unknowns, num_components = this % num_components)
       dq = 0.0_dp
       call this % inner % solve(-residual, dq, linear_achieved)

       x = x + dq

    end do

    call this % matvec(x, y)
    achieved = this % norm(y + g - rhs)

  end subroutine solve

end module operation_newton
