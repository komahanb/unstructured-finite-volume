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
!                  THE HALLEY-CHEBYSHEV FAMILY
!
! A statement that reports max_degree above one can be asked for more
! than its tangent, and higher_order_jacobian_product is how much of
! it a run asks for. Above the plain Newton step delta_1, solving
! J delta_1 = -R, each further order s = 2, ..., p adds
!
!      J delta_s  =  - B^(s) / s! ,
!
! B^(s) the s-th total derivative of R composed with the path whose
! k-th derivative is k! delta_k for k < s - the chain rule's own
! composition, assembled by operation_chain_rule - with delta_s left
! unoccupied, so the one term that would need it is not assembled;
! that missing term is exactly the J delta_s being solved for. The
! step taken is delta_1 + delta_2 + ... + delta_p. p = 1 is Newton
! unchanged; p = 2 is Halley's method, cubically convergent; each
! further p adds one more derivative of R and one more solve against
! the SAME jacobian, already frozen and, where the inner minimizer
! factors, already factored - the higher orders are additional right
! hand sides against one linear system, not a second one.
!
! The expansion is asymptotic: it is trusted only where delta_1 is
! already a fair local model of the root, and far from there a
! higher-order term can be larger than the one before it rather than
! smaller, which is the expansion leaving the regime it describes
! rather than refining within it. Each delta_s is kept only while it
! is no bigger than delta_(s-1); the first that is not stops the
! correction there, so the step taken is never worse than the Newton
! step this extends, and the check is on the correction's own
! decline, not a magnitude chosen from outside it.
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
  use field_calculus, only : field
  use operation_linearization, only : linearization, tangent_of
  use operation_chain_rule   , only : total_derivative, derivative_of, argument_path

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

     ! The order of the Halley-Chebyshev family taken: one is Newton
     ! unchanged, and the statement's own max_degree is the ceiling
     ! on how far above one this may be asked to go.
     integer :: higher_order_jacobian_product = 1

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
    real(dp), allocatable :: residual(:), y(:), dq(:)
    real(dp) :: linear_achieved
    integer :: it

    call tally_record(newton_solves)

    allocate(dq(size(x)))

    call this % begin_imbalance()

    ! the tangent in the unknown's argument; which road it takes is
    ! the statement's own answer, and no dispatch lives here
    jacobian = tangent_of(this % action, this % action % argument(1))

    do it = 1, this % max_iterations

       call tally_record(primal_loops)

       ! Where we stand: the full statement, whatever its linearity,
       ! on the input tuple every tangent below is frozen on.
       call this % raw_apply(x, y, inputs)
       residual = y - rhs

       achieved = this % norm(residual)

       ! Met; or no longer a number; or past what the arithmetic
       ! holds; or diverging; or, where the budget is taken from the
       ! rate, flattened - the slope of the residual's logarithm
       ! against its own scatter, no floor named. One question.
       if (this % halted(achieved, it)) return

       ! The linear question at this point, answered by the governed
       ! minimizer: the Jacobian is frozen at the same input tuple the
       ! residual was evaluated on, held inputs included.
       call jacobian % freeze(inputs, base=y)

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
          call compiled % stamped(this % action % stamp(), this % action % stamp_transposed())
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

       if (this % higher_order_jacobian_product > 1) then
          call halley_correction(this, inputs, dq, linear_achieved)
          if (linear_achieved /= linear_achieved) return
          if (linear_achieved > huge(1.0_dp) / 2.0_dp) return
       end if

       x = x + dq

    end do

    call this % raw_apply(x, y)
    achieved = this % norm(y - rhs)

  end subroutine solve

  !===================================================================!
  ! Add delta_2, ..., delta_p to the Newton step delta already
  ! solved, each against the same frozen jacobian this % inner is
  ! already attached to. achieved is the worst of the extra solves,
  ! read by the same guard the caller applies to the Newton one.
  !===================================================================!

  subroutine halley_correction(this, inputs, delta, achieved)

    class(newton)      , intent(inout) :: this
    type(stored_field) , intent(in)    :: inputs(:)
    real(dp)           , intent(inout) :: delta(:)
    real(dp)           , intent(out)   :: achieved

    type(total_derivative) :: total
    type(argument_path) :: path
    class(field), allocatable :: out
    type(stored_field) :: seeded
    real(dp), allocatable :: b(:), correction(:), individual(:,:)
    real(dp) :: fact, one_achieved
    integer :: s, p

    p = this % higher_order_jacobian_product
    achieved = 0.0_dp


    allocate(individual(size(delta), p))
    individual(:, 1) = delta

    allocate(correction(size(delta)))
    path = argument_path(this % action % argument(1), p - 1)
    fact = 1.0_dp

    do s = 2, p

       ! derivative(m) carries m! delta_m; fact holds (s-1)! at the
       ! point derivative(s-1) is set, so this line and no other
       ! needs to know a factorial's value.
       fact = fact * real(s - 1, dp)
       seeded = stored_field('correction', this % unknown_domain, size(delta))
       call seeded % set_real_vector(fact * individual(:, s - 1))
       path % derivative(s - 1) % occupied  = .true.
       path % derivative(s - 1) % direction = seeded

       ! the derivative of the statement, of this order, along this
       ! path: an operation, applied like any other
       total = derivative_of(this % action, s, [path])
       call total % apply(this % on, inputs, out)
       call out % real_vector(b)

       correction = 0.0_dp
       call this % inner % solve(-b / (fact * real(s, dp)), correction, one_achieved)
       achieved = max(achieved, one_achieved)
       if (one_achieved /= one_achieved) return
       if (one_achieved > huge(1.0_dp) / 2.0_dp) return

       ! The series is trusted only while it is shrinking: a
       ! correction no smaller than the one before it says the local
       ! model has left the regime a truncated expansion describes,
       ! and adding it would perturb rather than refine. What was
       ! already accumulated is kept - at s = 2 that is delta_1,
       ! Newton's own step, so this can never do worse than Newton.
       if (norm2(correction) >= norm2(individual(:, s - 1))) return

       individual(:, s) = correction
       delta = delta + correction

    end do

  end subroutine halley_correction

end module operation_newton
