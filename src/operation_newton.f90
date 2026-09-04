!=====================================================================!
! Newton's iteration: a minimizer like every other.
!
! One family, one protocol: state a statement, drive its residual to
! zero. The linear members reach the solution in the statement's own
! space; newton reaches it by linearising at the current iterate,
!
!      J(q) dq = rhs - action(q)          q <- q + dq
!
! and it is not a different kind of object for that - it extends the
! same base, implements the same operation interface, and implements
! the same solve(rhs, x, achieved). A network's training loop would
! join the family the same way: another concrete type of the one
! interface.
!
! Newton contains no derivative mathematics. The tangent is a level-1
! member - the linearisation operator - and newton only controls the
! iteration: freeze the linearisation at the current iterate, pass
! the linear system to the inner minimizer, step. The linearisation
! is determined by the statement's type: a differentiable statement
! linearises itself exactly, any other is differenced - the contract
! the family declares, implemented by one dispatch, and the iteration
! control below it never changes.
!
!                  THE HALLEY-CHEBYSHEV FAMILY
!
! A statement that reports max_degree above one can return more than
! its tangent, and higher_order_jacobian_product is the order a run
! requests. Above the plain Newton step delta_1, solving
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
! The expansion is asymptotic: it is valid only where delta_1 is
! already an accurate local model of the root, and far from there a
! higher-order term can be larger than the one before it rather than
! smaller, which is the expansion leaving the regime it describes
! rather than refining within it. Each delta_s is retained only while
! its norm is no larger than that of delta_(s-1); the first that is
! not ends the correction there, so the step taken is never longer
! than the Newton step this extends, and the check is on the
! correction's own decrease, not a magnitude chosen from outside it.
!
!                     WHEN THE ITERATION DIVERGES
!
! A statement need not have a solution near the initial iterate, and
! one that does not moves the iterate away from any root rather than
! toward one. Left to run, that reaches a jacobian the inner
! minimizer cannot factor, and a direct one stops the program there -
! so a caller loses the whole run to one statement that has no
! solution.
!
! So a residual that grows far past the residual of the initial
! iterate, or that is not a number, is classified as a diverging
! statement, and the iteration returns with that residual rather than
! continuing. The caller then receives what it would have received
! from any other failure to converge, which is a residual above
! tolerance, and may report so and continue.
!
! The factor is large by design. Newton's first steps on a difficult
! statement often increase the residual before it decreases, and this
! is not a line search; it is only the difference between reporting
! a failure and stopping the program.
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
  ! controls, one level down, passed one linear system per step.
  !===================================================================!

  type, extends(minimizer) :: newton

     ! Whether the tangent is taken explicit where the statement
     ! provides it. When false, the linearization is stated - a
     ! matrix-vector product, no stored matrix - and the inner
     ! minimizer must iterate.
     logical :: explicit = .true.

     ! The order of the Halley-Chebyshev family taken: one is Newton
     ! unchanged, and the statement's own max_degree is the upper
     ! bound on this order.
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
    type(stencil) :: tangent
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

    ! the tangent in the unknown's argument; which mode it uses is
    ! determined by the statement, and no dispatch is defined here
    jacobian = tangent_of(this % action, this % action % argument(1))

    do it = 1, this % max_iterations

       call tally_record(primal_loops)

       ! The current iterate: the full statement, of any linearity,
       ! evaluated on the input tuple every tangent below is frozen on.
       call this % evaluate(x, y, inputs)
       residual = y - rhs

       achieved = this % norm(residual)

       ! Converged; or not a number; or past the range of the
       ! arithmetic; or diverging; or, where the iteration limit is
       ! derived from the rate, stagnant - the slope of the residual's
       ! logarithm compared with its own variance, no absolute lower
       ! bound named. One predicate.
       if (this % halted(achieved, it)) return

       ! The linear system at this iterate, solved by the inner
       ! minimizer: the Jacobian is frozen at the same input tuple the
       ! residual was evaluated on, stored inputs included.
       call jacobian % freeze(inputs, base=y)

       ! A statement whose tangent is explicit passes the inner
       ! minimizer a stencil, whose pattern is then the coupling a
       ! structured minimizer sweeps by; any other is passed the
       ! linearization, a matrix-vector product.
       available = .false.
       if (this % explicit) then
          call this % action % explicit_tangent(this % graph, this % action % bind(inputs), &
               & 1, rows, columns, &
               & weights, available)
       end if
       if (available) then
          tangent = stencil(rows, columns, weights, &
               & spread(0.0_dp, 1, this % num_unknowns), 'explicit tangent')
          call tangent % versioned(this % action % version(), this % action % transpose_version())
          call this % inner % state(tangent, tangent % pattern, this % unknown_domain, &
               & this % num_unknowns, num_components = this % num_components, &
               & coupling = tangent % pattern)
       else
          call this % inner % state(jacobian, this % graph, this % unknown_domain, &
               & this % num_unknowns, num_components = this % num_components)
       end if
       dq = 0.0_dp
       call this % inner % solve(-residual, dq, linear_achieved)

       ! An inner minimizer that encountered a singular tangent reports
       ! a residual no completed solve produces, and one that
       ! overflowed reports a value that is not a number. Neither
       ! yields a step.
       if (linear_achieved /= linear_achieved) return
       if (linear_achieved > huge(1.0_dp) / 2.0_dp) return

       if (this % higher_order_jacobian_product > 1) then
          call halley_correction(this, inputs, dq, linear_achieved)
          if (linear_achieved /= linear_achieved) return
          if (linear_achieved > huge(1.0_dp) / 2.0_dp) return
       end if

       x = x + dq

    end do

    call this % evaluate(x, y)
    achieved = this % norm(y - rhs)

  end subroutine solve

  !===================================================================!
  ! Add delta_2, ..., delta_p to the Newton step delta already
  ! solved, each against the same frozen jacobian this % inner is
  ! already stated with. achieved is the maximum over the additional
  ! solves, checked by the same guard condition the caller applies to
  ! the Newton solve.
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

       ! derivative(m) stores m! delta_m; fact stores (s-1)! when
       ! derivative(s-1) is set, so this line and no other requires a
       ! factorial's value.
       fact = fact * real(s - 1, dp)
       seeded = stored_field('correction', this % unknown_domain, size(delta))
       call seeded % set_real_vector(fact * individual(:, s - 1))
       path % derivative(s - 1) % occupied  = .true.
       path % derivative(s - 1) % direction = seeded

       ! the derivative of the statement, of this order, along this
       ! path: an operation, applied like any other
       total = derivative_of(this % action, s, [path])
       call total % apply(this % graph, total % bind(inputs), out)
       call out % real_vector(b)

       correction = 0.0_dp
       call this % inner % solve(-b / (fact * real(s, dp)), correction, one_achieved)
       achieved = max(achieved, one_achieved)
       if (one_achieved /= one_achieved) return
       if (one_achieved > huge(1.0_dp) / 2.0_dp) return

       ! The series is valid only while it is decreasing: a
       ! correction no smaller than the one before it indicates the
       ! local model has left the regime a truncated expansion
       ! describes, and adding it would perturb rather than refine.
       ! What was already accumulated is retained - at s = 2 that is
       ! delta_1, Newton's own step, so this is never larger than
       ! Newton.
       if (norm2(correction) >= norm2(individual(:, s - 1))) return

       individual(:, s) = correction
       delta = delta + correction

    end do

  end subroutine halley_correction

end module operation_newton
