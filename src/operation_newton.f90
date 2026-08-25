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
!                        WHEN IT IS NOT WORKING
!
! A statement need not have a solution near where it was started, and
! one that does not sends the iterate away rather than toward it. Left
! to run, that reaches a jacobian the inner minimizer cannot factor,
! and a direct one stops the program there - so a caller loses the
! whole run to one statement that was never going to be solved.
!
! So a residual that grows far past the one the first guess gave, or
! that is no longer a number at all, is taken as a statement
! diverging, and the iteration returns with that residual rather than
! pursuing it. The caller then sees what it
! would have seen from any other failure to converge, which is a
! number too large, and may say so and carry on.
!
! The factor is wide on purpose. Newton's first steps on a hard
! statement often climb before they fall, and this is not a line
! search; it is only the difference between reporting a failure and
! becoming one.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_newton

  use util_precision  , only : dp
  use operation_minimization        , only : minimizer
  use operation_stencil     , only : stencil
  use util_tally, only : tally_record, newton_solves, primal_loops
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

     ! Whether the tangent is taken compiled where the statement
     ! offers it. Off, the linearization is attached - a matvec, no
     ! matrix anywhere - and the inner minimizer must iterate.
     logical :: compiled = .true.

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
    type(stencil) :: compiled
    type(stored_field), allocatable :: inputs(:)
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:)
    logical :: available
    real(dp), allocatable :: residual(:), g(:), y(:), dq(:)
    real(dp) :: linear_achieved
    integer :: it

    call tally_record(newton_solves)

    allocate(dq(size(x)))

    call this % begin_imbalance()

    call this % constant(g)

    ! the tangent in the unknown's argument; which road it takes is
    ! the statement's own answer, and no dispatch lives here
    jacobian = tangent_of(this % action, this % action % argument(1))

    do it = 1, this % max_iterations

       call tally_record(primal_loops)

       ! Where we stand: the full statement, whatever its linearity.
       call this % matvec(x, y)
       residual = y + g - rhs

       achieved = this % norm(residual)

       ! Met; or no longer a number; or past what the arithmetic
       ! holds; or diverging; or, where the budget is taken from the
       ! rate, flattened - the slope of the residual's logarithm
       ! against its own scatter, no floor named. One question.
       if (this % halted(achieved, it)) return

       ! The linear question at this point, answered by the governed
       ! minimizer: the Jacobian is frozen at the same input tuple the
       ! residual was evaluated on, held inputs included.
       call this % evaluation_inputs(x, inputs)
       call jacobian % freeze(inputs, base=y + g)

       ! A statement that compiles its tangent hands the inner
       ! minimizer a stencil, whose pattern is then the coupling a
       ! structured minimizer sweeps by; any other is handed the
       ! linearization, a matvec.
       available = .false.
       if (this % compiled) then
          call this % action % compiled_tangent(this % on, inputs, 1, rows, columns, &
               & weights, available)
       end if
       if (available) then
          compiled = stencil(rows, columns, weights, &
               & spread(0.0_dp, 1, this % num_unknowns), 'compiled tangent')
          call compiled % stamped(this % action % stamp())
          call this % inner % attach(compiled, compiled % pattern, this % unknown_domain, &
               & this % num_unknowns, num_components = this % num_components, &
               & coupling = compiled % pattern)
       else
          call this % inner % attach(jacobian, this % on, this % unknown_domain, &
               & this % num_unknowns, num_components = this % num_components)
       end if
       dq = 0.0_dp
       call this % inner % solve(-residual, dq, linear_achieved)

       ! An inner minimizer that met a singular tangent reports a
       ! residual no completed solve produces, and one that overflowed
       ! reports no number at all. Neither leaves a step to take.
       if (linear_achieved /= linear_achieved) return
       if (linear_achieved > huge(1.0_dp) / 2.0_dp) return

       x = x + dq

    end do

    call this % matvec(x, y)
    achieved = this % norm(y + g - rhs)

  end subroutine solve

end module operation_newton
