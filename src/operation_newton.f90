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
! Newton contains no derivative mathematics. The jacobian is a level-1
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
! its jacobian, and higher_order_jacobian_product is the order a run
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

  use iso_fortran_env , only : int64
  use util_precision  , only : dp
  use operation_minimization        , only : minimizer, restrict, solve_result, SOLVE_INNER_FAILED
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use operation_stencil     , only : stencil
  use util_tally, only : newton_solves, primal_loops, tally
  use util_verbosity, only : verbosity
  use operation_dense_direct, only : dense_direct
  use field_stored  , only : stored_field
  use field_calculus, only : field
  use operation_linearization, only : linearization, jacobian_of
  use operation_chain_rule   , only : total_derivative, derivative_of, argument_path

  implicit none

  private
  public :: newton

  !===================================================================!
  ! Newton: one component beyond the family - the minimizer it
  ! controls, one level down, passed one linear system per step.
  !===================================================================!

  ! how many halvings of a Newton step the line search tries before
  ! it takes the last trial
  integer, parameter :: max_halvings = 8

  type, extends(minimizer) :: newton

     ! Whether the jacobian is taken explicit where the statement
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
     procedure :: restrict => newton_restrict
     procedure :: bind_account => newton_bind_account
     procedure :: storage_entries => newton_storage_entries
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
  ! The linear system of every step is over the same unknowns as the
  ! statement, so the inner minimizer is restricted by the selection
  ! itself.
  !===================================================================!

  subroutine newton_restrict(this, selected)

    class(newton), intent(inout) :: this
    integer      , intent(in)    :: selected(:)

    call restrict(this, selected)
    if (allocated(this % inner)) call this % inner % restrict(selected)

  end subroutine newton_restrict

  ! the inner minimizer's entries over the same unknowns
  pure integer(int64) function newton_storage_entries(this, num_unknowns) result(entries)

    class(newton), intent(in) :: this
    integer      , intent(in) :: num_unknowns

    entries = 0_int64
    if (allocated(this % inner)) entries = this % inner % storage_entries(num_unknowns)

  end function newton_storage_entries

  !===================================================================!
  ! Drive action(q) toward rhs from the given q.
  !===================================================================!

  subroutine solve(this, rhs, x, achieved)

    class(newton), intent(inout) :: this
    real(dp), intent(in)    :: rhs(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out)   :: achieved

    type(linearization) :: linearized
    type(stencil) :: jacobian
    type(stored_field), allocatable :: inputs(:)
    integer , allocatable :: rows(:), columns(:)
    real(dp), allocatable :: weights(:)
    logical :: jacobian_defined
    real(dp), allocatable :: residual(:), y(:), dq(:)
    real(dp) :: linear_achieved, initial, factor, previous, step_norm
    integer :: nonzeros
    real(dp), allocatable :: trial(:)
    integer :: it, halving
    type(solve_result) :: outcome

    call this % record_event(newton_solves)

    allocate(dq(size(x)))

    call this % initialize_residual_history()

    call this % evaluate(x, y, inputs)
    residual = y - rhs
    achieved = this % norm(residual)
    initial  = achieved
    if (verbosity >= 1) then
       print '(a)', ' step          |F|     |F|/|F0|          |dq|   factor    rows  nonzeros  jacobian  halley' // &
            & '   cycles          |r|     |r|/|r0|  linear solve'
       print '(i5,es13.4)', 0, achieved
    end if
    if (this % terminated(achieved, 0)) return

    ! the jacobian in the unknown's argument; which mode it uses is
    ! determined by the statement, and no dispatch is defined here
    linearized = jacobian_of(this % action, this % action % argument(1))

    do it = 1, this % max_iterations

       call this % record_event(primal_loops)

       ! The linear system at this iterate, solved by the inner
       ! minimizer: the Jacobian is frozen at the same input tuple the
       ! residual was evaluated on, stored inputs included.
       call linearized % freeze(inputs, base=y)

       ! A statement whose jacobian is explicit passes the inner
       ! minimizer a stencil, whose pattern is then the coupling a
       ! structured minimizer sweeps by; any other is passed the
       ! linearization, a matrix-vector product.
       jacobian_defined = .false.
       if (this % explicit) then
          call this % action % explicit_jacobian(this % graph, this % action % bind(inputs), &
               & 1, rows, columns, &
               & weights, jacobian_defined)
       end if
       nonzeros = 0
       if (jacobian_defined) nonzeros = size(rows)
       if (jacobian_defined) then
          jacobian = stencil(rows, columns, weights, &
               & spread(0.0_dp, 1, this % num_unknowns), 'explicit jacobian')
          call jacobian % versioned(this % action % version(), this % action % transpose_version())
          if (allocated(this % distribution)) this % inner % distribution = this % distribution
          call this % inner % state(jacobian, jacobian % pattern, this % unknown_domain, &
               & this % num_unknowns, num_components = this % num_components, &
               & coupling = jacobian % pattern)
       else
          call this % inner % state(linearized, this % graph, this % unknown_domain, &
               & this % num_unknowns, num_components = this % num_components)
       end if
       dq = 0.0_dp
       call this % inner % solve(-residual, dq, linear_achieved)
       ! the step of a distributed solve is owned in parts: every
       ! image gathers the whole before the residual is evaluated
       if (allocated(this % distribution)) call this % distribution % gather(dq)

       ! An inner minimizer that encountered a singular jacobian reports
       ! a residual no completed solve produces, and one that
       ! overflowed reports a value that is not a number. Neither
       ! yields a step.
       outcome = this % inner % result()
       if (outcome % failed() .or. .not. ieee_is_finite(linear_achieved) .or. &
            & linear_achieved > huge(1.0_dp) / 2.0_dp) then
          call this % record_result(achieved, it - 1, SOLVE_INNER_FAILED)
          return
       end if

       if (this % higher_order_jacobian_product > 1) then
          call halley_correction(this, inputs, dq, linear_achieved)
          outcome = this % inner % result()
          if (outcome % failed() .or. .not. ieee_is_finite(linear_achieved) .or. &
               & linear_achieved > huge(1.0_dp) / 2.0_dp) then
             call this % record_result(achieved, it - 1, SOLVE_INNER_FAILED)
             return
          end if
       end if

       ! THE STEP IS DAMPED where the full step does not decrease the
       ! residual: the factor is halved until |F(x + factor dq)| is
       ! below (1 - 1e-4 factor) |F(x)|, the Armijo condition, or the
       ! halvings are exhausted, after which the last trial is taken.
       ! The reported residual and the next jacobian use the same
       ! updated state. The count is the number of completed steps.
       previous = achieved
       factor   = 1.0_dp
       do halving = 0, max_halvings
          trial = x + factor * dq
          call this % evaluate(trial, y, inputs)
          residual = y - rhs
          achieved = this % norm(residual)
          if (achieved <= (1.0_dp - 1.0e-4_dp * factor) * previous) exit
          if (halving == max_halvings) exit
          factor = 0.5_dp * factor
       end do
       x  = trial
       dq = factor * dq
       ! the norm of the step is a reduction over the images, taken by
       ! every image whether or not it prints
       step_norm = this % norm(dq)
       if (verbosity >= 1) call reported(this, it, achieved, initial, step_norm, factor, jacobian_defined, nonzeros, outcome)
       if (this % terminated(achieved, it)) return

    end do

  end subroutine solve

  !===================================================================!
  ! One row per Newton step at verbosity one, under a header printed
  ! with step zero: the residual after the step, absolute and relative
  ! to the first, the step's norm and the factor the line search took
  ! it by, the jacobian's order (rows, one per unknown), its nonzeros
  ! and its representation - sparse, the explicit stencil formed at
  ! this step; dense, that stencil compiled to a full matrix by a
  ! direct inner solver; free, no matrix, the directional derivatives
  ! of the linearization - the Halley correction's order or its
  ! absence, and the inner solve's restart cycles, residual absolute
  ! and relative to its own first, and outcome.
  !===================================================================!

  subroutine reported(this, it, achieved, initial, step_norm, factor, jacobian_defined, nonzeros, outcome)

    class(newton)     , intent(in) :: this
    integer           , intent(in) :: it
    real(dp)          , intent(in) :: achieved, initial, step_norm, factor
    logical           , intent(in) :: jacobian_defined
    integer           , intent(in) :: nonzeros
    type(solve_result), intent(in) :: outcome

    character(len=12) :: jacobian, halley
    real(dp) :: relative, linear_relative

    jacobian = 'free'
    if (jacobian_defined) then
       jacobian = 'sparse'
       select type (inner => this % inner)
       type is (dense_direct)
          jacobian = 'dense'
       end select
    end if
    halley = 'none'
    if (this % higher_order_jacobian_product > 1) then
       write(halley,'(a,i0)') 'order ', this % higher_order_jacobian_product
    end if
    relative = 0.0_dp
    if (initial > 0.0_dp) relative = achieved / initial
    linear_relative = 0.0_dp
    if (outcome % initial_residual > 0.0_dp) linear_relative = outcome % residual / outcome % initial_residual
    print '(i5,3es13.4,f9.4,i8,i10,2x,a8,2x,a6,i7,2es13.4,2x,a)', it, achieved, relative, step_norm, factor, &
         & this % num_unknowns * this % num_components, nonzeros, jacobian, halley, outcome % iterations, &
         & outcome % residual, linear_relative, trim(outcome % description())

  end subroutine reported

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
    type(solve_result) :: outcome

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
       outcome = this % inner % result()
       achieved = max(achieved, one_achieved)
       if (outcome % failed()) return
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

  !===================================================================!
  ! Bind the account of this minimizer and of its children.
  !===================================================================!

  subroutine newton_bind_account(this, account)

    class(newton), intent(inout) :: this
    type(tally), pointer, intent(in) :: account

    this % account => account
    if (allocated(this % inner)) call this % inner % bind_account(account)

  end subroutine newton_bind_account

end module operation_newton
