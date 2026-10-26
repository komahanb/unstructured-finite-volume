!=====================================================================!
! Module that states Newton's iteration: a minimizer like every other.
!
! One family, one protocol: state a statement, drive its residual to
! zero. The linear members reach the solution in the statement's own
! space; newton reaches it by linearising at the current iterate,
!
!      J(q) dq = rhs - action(q)          q <- q + dq
!
! and it is not a different kind of object for that: it extends the
! same base, implements the same operation interface, and implements
! the same solve(rhs, x, achieved). A network's training loop would
! join the family the same way, as another concrete type of the one
! interface.
!
! Newton contains no derivative mathematics. The jacobian is a level-1
! member, the linearisation operator, and newton only controls the
! iteration: freeze the linearisation at the current iterate, pass the
! linear system to the inner minimizer, step.
!
!      q --> F(q) = action(q) - rhs --> freeze J(q) --> inner solve
!      ^                                                 J dq = -F
!      |                                                     |
!      +---- q + factor dq <---- halve factor until the -----+
!                                residual decreases (Armijo)
!
! The linearisation is determined by the statement's type: a
! differentiable statement linearises itself exactly, any other is
! differenced. That is the contract the family declares, implemented
! by one dispatch, and the iteration control below it never changes.
!
! The Halley-Chebyshev family. A statement that reports max_degree
! above one can return more than its jacobian, and
! higher_order_jacobian_product is the order a run requests. Above the
! plain Newton step delta_1, solving J delta_1 = -R, each further order
! s = 2, ..., p adds
!
!      J delta_s  =  - B^(s) / s! ,
!
! B^(s) the s-th total derivative of R composed with the path whose
! k-th derivative is k! delta_k for k < s: the chain rule's own
! composition, assembled by operation_chain_rule, with delta_s left
! unoccupied, so the one term that would need it is not assembled. That
! missing term is exactly the J delta_s being solved for. The step
! taken is delta_1 + delta_2 + ... + delta_p.
!
!      J delta_1  =  -R              +-----------------------------+
!      J delta_2  =  -B^(2) / 2!     |  one jacobian, frozen, and  |
!         ...                        |  where the inner minimizer  |
!      J delta_p  =  -B^(p) / p!     |  factors, factored once     |
!                                    +-----------------------------+
!
! p = 1 is Newton unchanged; p = 2 is Halley's method, cubically
! convergent; each further p adds one more derivative of R and one
! more solve against the same jacobian. The higher orders are
! additional right hand sides against one linear system, not a second
! one.
!
! The expansion is asymptotic: it is valid only where delta_1 is
! already an accurate local model of the root, and far from there a
! higher-order term can be larger than the one before it rather than
! smaller, which is the expansion leaving the regime it describes
! rather than refining within it. Each delta_s is retained only while
! its norm is no larger than that of delta_(s-1); the first that is not
! ends the correction there, so the step taken is never longer than
! the Newton step this extends, and the check is on the correction's
! own decrease, not a magnitude chosen from outside it.
!
! When the iteration diverges. A statement need not have a solution
! near the initial iterate, and one that does not moves the iterate
! away from any root rather than toward one. Left to run, that reaches
! a jacobian the inner minimizer cannot factor, and a direct one stops
! the program there, so a caller loses the whole run to one statement
! that has no solution.
!
! So a residual that grows far past the residual of the initial
! iterate, or that is not a number, is classified by the family's
! termination test as a diverging statement, and the iteration returns
! with that residual rather than continuing. The caller then receives
! what it would have received from any other failure to converge,
! which is a residual above tolerance, and may report so and continue.
!
! The factor is large by design. Newton's first steps on a difficult
! statement often increase the residual before it decreases, and this
! is not a line search; it is only the difference between reporting a
! failure and stopping the program.
!
! The chord method. With chord, the jacobian formed at the first step,
! and under a direct inner minimizer its factorisation, serves the
! steps after it, J_0 dq = -F(q), each a residual and a substitution
! where Newton forms and factorises the jacobian of every iterate. The
! rate is linear, ||I - J_0^-1 J(q)||, where Newton's is quadratic. A
! step the line search shortens is one along which J_0 no longer gives
! a decrease of the residual, and the jacobian of the step after it is
! formed again:
!
!      step          1      2      3      4      5
!      jacobian      J_0    J_0    J_0    J_3    J_3
!      factor        1      1      1/2    1      1
!                                   |      ^
!                                   +------+  shortened, so formed again
!
! The Halley correction is of the jacobian at the iterate, and is not
! taken with it.
!
! At the rounding of the residual. The tolerance is relative to the
! residual the iteration began at. An iterate already near the
! solution, the solution at the instant before in a march, begins at a
! residual so small that a reduction by the tolerance lies below the
! rounding of the residual's own evaluation, and no step reaches it.
!
!      |F|
!       | *
!       |    *  *  *  *  *      the rounding of F's evaluation
!       |
!       |- - - - - - - - - -    tolerance x |F_0|: no step reaches it
!       +--------------------- step
!
! There no halving of the step decreases the residual, and the full
! step is below half the digits of the iterate: the iterate is the
! solution to the rounding, and the iteration ends converged at it.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_newton

  use iso_fortran_env        , only : int64
  use util_precision         , only : dp, half_digits
  use operation_minimizer    , only : minimizer, restrict, solve_result, SOLVE_INNER_FAILED, SOLVE_CONVERGED
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use operation_stencil      , only : stencil
  use util_tally             , only : newton_solves, primal_loops, tally
  use util_verbosity         , only : verbosity
  use operation_dense_direct , only : dense_direct
  use field_stored           , only : stored_field
  use field_calculus         , only : field
  use operation_linearization, only : linearization, jacobian_of
  use operation_chain_rule   , only : total_derivative, derivative_of, argument_path

  implicit none

  private

  public :: newton

  ! How many halvings of a Newton step the line search tries before it
  ! takes the last trial
  integer, parameter :: max_halvings = 8

  ! The jacobians the chord method forms, numbered from one: a direct
  ! inner minimizer retains a factorisation for the jacobian of its own
  ! number alone, and factorises one of number zero at every solve
  integer, save :: jacobians_formed = 0

  !-------------------------------------------------------------------!
  ! Type that models Newton's iteration: one component beyond the
  ! family, the minimizer it controls one level down, passed one linear
  ! system per step
  !
  !      newton            the iteration: linearise, step, damp
  !        `-- inner       the minimizer of each step's linear system
  !-------------------------------------------------------------------!

  type, extends(minimizer) :: newton

     ! Whether the jacobian is taken explicit where the statement
     ! provides it. When false, the linearization is stated, a
     ! matrix-vector product with no stored matrix, and the inner
     ! minimizer must iterate.
     logical :: explicit = .true.

     ! The order of the Halley-Chebyshev family taken: one is Newton
     ! unchanged, and the statement's own max_degree is the upper bound
     ! on this order
     integer :: higher_order_jacobian_product = 1

     ! The chord method: the jacobian of the first step for the steps
     ! after it, formed again after a step the line search shortens
     logical :: chord = .false.

     ! The minimizer of each step's linear system
     class(minimizer), allocatable :: inner

   contains

     procedure :: name            => newton_name
     procedure :: restrict        => newton_restrict
     procedure :: bind_account    => newton_bind_account
     procedure :: storage_entries => newton_storage_entries
     procedure :: solve

  end type newton

  ! A newton stated by its inner minimizer and its options
  interface newton
     module procedure newton_of
  end interface newton

contains

  !===================================================================!
  ! Returns the name of this minimizer
  !===================================================================!

  pure function newton_name(this) result(name)

    class(newton)   , intent(in)  :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'newton'

  end function newton_name

  !===================================================================!
  ! Restricts the minimizer to the selected unknowns. The linear system
  ! of every step is over the same unknowns as the statement, so the
  ! inner minimizer is restricted by the selection itself.
  !===================================================================!

  subroutine newton_restrict(this, selected)

    class(newton), intent(inout) :: this
    integer      , intent(in)    :: selected(:)

    call restrict(this, selected)
    if (allocated(this % inner)) call this % inner % restrict(selected)

  end subroutine newton_restrict

  !===================================================================!
  ! Returns the inner minimizer's stored entries over the same unknowns
  !===================================================================!

  pure integer(int64) function newton_storage_entries(this, num_unknowns) result(entries)

    class(newton), intent(in) :: this
    integer      , intent(in) :: num_unknowns

    entries = 0_int64
    if (allocated(this % inner)) entries = this % inner % storage_entries(num_unknowns)

  end function newton_storage_entries

  !===================================================================!
  ! Drives action(q) toward rhs from the given q: Newton's iteration,
  ! each step's linear system passed to the inner minimizer, the step
  ! damped by halving, until the family's termination test is met
  !===================================================================!

  subroutine solve(this, rhs, x, achieved)

    class(newton), intent(inout) :: this
    real(dp)     , intent(in)    :: rhs(:)
    real(dp)     , intent(inout) :: x(:)
    real(dp)     , intent(out)   :: achieved

    type(linearization)             :: linearized
    type(stencil)                   :: jacobian
    type(stored_field), allocatable :: inputs(:)
    integer           , allocatable :: rows(:), columns(:)
    real(dp)          , allocatable :: weights(:)
    logical                         :: jacobian_defined, failed, decreased, forming
    real(dp)          , allocatable :: residual(:), y(:), dq(:)
    real(dp)                        :: linear_achieved, initial, factor, previous, step_norm, full_step
    integer                         :: nonzeros
    real(dp)          , allocatable :: trial(:)
    integer                         :: it, halving
    type(solve_result)              :: outcome

    call this % record_event(newton_solves)

    if (this % chord .and. this % higher_order_jacobian_product > 1) then
       error stop 'operation_newton: the chord method solves against the jacobian of an earlier step; the &
            &Halley correction is formed with the jacobian at the step'
    end if

    allocate(dq(size(x)))

    call this % initialize_residual_history()

    call this % evaluate(x, y, inputs)
    residual = y - rhs
    achieved = this % norm(residual)
    initial  = achieved
    if (verbosity >= 1) then
       print '(a)', ' step          |F|     |F|/|F0|          |dq|   factor    rows  nonzeros  jacobian' // &
            & '  factorised  halley   cycles          |r|     |r|/|r0|  linear solve'
       print '(i5,es13.4)', 0, achieved
    end if
    if (this % terminated(achieved, 0)) return

    ! The jacobian in the unknown's argument; which mode it uses is
    ! determined by the statement, and no dispatch is defined here
    linearized = jacobian_of(this % action, this % action % argument(1))
    jacobian_defined = .false.
    nonzeros = 0
    factor   = 1.0_dp

    do it = 1, this % max_iterations

       call this % record_event(primal_loops)

       ! The jacobian is formed at this step: at every step, and under
       ! the chord method at the first and after a step the line search
       ! shortened; otherwise the inner minimizer solves the system it
       ! was stated with
       forming = .not. this % chord .or. it == 1 .or. factor < 1.0_dp

       if (forming) then
          ! The linear system at this iterate, solved by the inner
          ! minimizer: the jacobian is frozen at the same input tuple the
          ! residual was evaluated on, stored inputs included
          call linearized % freeze(inputs, base=y)

          ! A statement whose jacobian is explicit passes the inner
          ! minimizer a stencil. The coupling a structured minimizer
          ! sweeps by is the one this solver was stated with, when one
          ! was, and the stencil's pattern otherwise: a row that reads
          ! every unknown, the integral of a field over the region, lies
          ! in the pattern but couples no two points, and a colouring
          ! taken from the pattern gives its block a colour of its own
          ! and the others one colour fewer, under which the sweep
          ! converges at half the rate (255 against 117 preconditioner
          ! applications on the 16 x 16 vortex). Any other statement is
          ! passed the linearization, a matrix-vector product.
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
             if (this % chord) then
                jacobians_formed = jacobians_formed + 1
                call jacobian % versioned(jacobians_formed)
             else
                call jacobian % versioned(this % action % version(), this % action % transpose_version())
             end if
             if (allocated(this % distribution)) this % inner % distribution = this % distribution
             if (allocated(this % coupling)) then
                call this % inner % state(jacobian, jacobian % pattern, this % unknown_domain, &
                     & this % num_unknowns, num_components = this % num_components, &
                     & coupling = this % coupling)
             else
                call this % inner % state(jacobian, jacobian % pattern, this % unknown_domain, &
                     & this % num_unknowns, num_components = this % num_components, &
                     & coupling = jacobian % pattern)
             end if
          else
             call this % inner % state(linearized, this % graph, this % unknown_domain, &
                  & this % num_unknowns, num_components = this % num_components)
          end if
       end if
       dq = 0.0_dp
       call this % inner % solve(-residual, dq, linear_achieved)
       ! The step of a distributed solve is owned in parts: every image
       ! gathers the whole before the residual is evaluated
       if (allocated(this % distribution)) call this % distribution % gather(dq)

       ! An inner minimizer that encountered a singular jacobian reports
       ! a residual no completed solve produces, and one that overflowed
       ! reports a value that is not a number. Neither yields a step.
       outcome = this % inner % result()
       failed = outcome % failed() .or. .not. ieee_is_finite(linear_achieved) .or. &
            & linear_achieved > huge(1.0_dp) / 2.0_dp
       if (.not. failed .and. this % higher_order_jacobian_product > 1) then
          call halley_correction(this, inputs, dq, linear_achieved)
          outcome = this % inner % result()
          failed = outcome % failed() .or. .not. ieee_is_finite(linear_achieved) .or. &
               & linear_achieved > huge(1.0_dp) / 2.0_dp
       end if
       if (failed) then
          call this % record_result(achieved, it - 1, SOLVE_INNER_FAILED)
          return
       end if

       ! The step is damped where the full step does not decrease the
       ! residual: the factor is halved until |F(x + factor dq)| is
       ! below (1 - 1e-4 factor) |F(x)|, the Armijo condition, or the
       ! halvings are exhausted, after which the last trial is taken.
       ! The reported residual and the next jacobian use the same
       ! updated state. The count is the number of completed steps.
       !
       !      factor    1     1/2     1/4    ...    1/2^8
       !      trial     x + factor dq, until the residual decreases
       previous  = achieved
       factor    = 1.0_dp
       full_step = this % norm(dq)
       do halving = 0, max_halvings
          trial = x + factor * dq
          call this % evaluate(trial, y, inputs)
          residual  = y - rhs
          achieved  = this % norm(residual)
          decreased = achieved <= (1.0_dp - 1.0e-4_dp * factor) * previous
          if (decreased) exit
          if (halving == max_halvings) exit
          factor = 0.5_dp * factor
       end do
       x  = trial
       dq = factor * dq
       ! The norm of the step is a reduction over the images, taken by
       ! every image whether or not it prints
       step_norm = this % norm(dq)
       if (verbosity >= 1) call reported(this, it, achieved, initial, step_norm, factor, jacobian_defined, nonzeros, &
            & forming, outcome)
       ! At the rounding of the residual: no halving of the step
       ! decreases the residual, and the full step is below half the
       ! digits of the iterate. The residual is the rounding of its own
       ! evaluation, which no step reduces, and the iterate before the
       ! step is the solution: converged, whatever the residual the
       ! iteration began at
       if (.not. decreased .and. full_step <= half_digits * this % norm(x)) then
          x        = x - dq
          achieved = previous
          call this % record_result(achieved, it, SOLVE_CONVERGED)
          return
       end if
       if (this % terminated(achieved, it)) return

    end do

  end subroutine solve

  !===================================================================!
  ! Prints one row per Newton step at verbosity one, under a header
  ! printed with step zero:
  !
  !      step             the step, zero for the initial residual
  !      |F|, |F|/|F0|    the residual after the step, absolute and
  !                       relative to the first
  !      |dq|, factor     the step's norm, and the factor the line
  !                       search took it by
  !      rows             the jacobian's order, one row per unknown
  !                       component
  !      nonzeros         its entries
  !      jacobian         sparse   the explicit stencil formed at this
  !                                step
  !                       dense    that stencil compiled to a full
  !                                matrix by a direct inner minimizer
  !                       free     no matrix: the directional
  !                                derivatives of the linearization
  !      factorised       yes      formed at this step, and by a direct
  !                                inner minimizer factorised
  !                       no       solved against the jacobian of an
  !                                earlier step (the chord method)
  !      halley           the Halley correction's order, or none
  !      cycles, |r|,     the inner solve's restart cycles, residual
  !      |r|/|r0|,        absolute and relative to its own first, and
  !      linear solve     outcome
  !===================================================================!

  subroutine reported(this, it, achieved, initial, step_norm, factor, jacobian_defined, nonzeros, formed, outcome)

    class(newton)     , intent(in) :: this
    integer           , intent(in) :: it
    real(dp)          , intent(in) :: achieved, initial, step_norm, factor
    logical           , intent(in) :: jacobian_defined, formed
    integer           , intent(in) :: nonzeros
    type(solve_result), intent(in) :: outcome

    character(len=12) :: jacobian, halley, factorised
    real(dp)          :: relative, linear_relative

    jacobian = 'free'
    if (jacobian_defined) then
       jacobian = 'sparse'
       select type (inner => this % inner)
       type is (dense_direct)
          jacobian = 'dense'
       end select
    end if
    factorised = merge('yes', 'no ', formed)
    halley = 'none'
    if (this % higher_order_jacobian_product > 1) then
       write(halley,'(a,i0)') 'order ', this % higher_order_jacobian_product
    end if
    relative = 0.0_dp
    if (initial > 0.0_dp) relative = achieved / initial
    linear_relative = 0.0_dp
    if (outcome % initial_residual > 0.0_dp) linear_relative = outcome % residual / outcome % initial_residual
    print '(i5,3es13.4,f9.4,i8,i10,2x,a8,2x,a10,2x,a6,i7,2es13.4,2x,a)', it, achieved, relative, step_norm, factor, &
         & this % num_unknowns * this % num_components, nonzeros, jacobian, factorised, halley, outcome % iterations, &
         & outcome % residual, linear_relative, trim(outcome % description())

  end subroutine reported

  !===================================================================!
  ! Adds delta_2, ..., delta_p to the Newton step delta already solved,
  ! each against the same frozen jacobian this % inner is already
  ! stated with. achieved is the maximum over the additional solves,
  ! checked by the same guard condition the caller applies to the
  ! Newton solve.
  !
  !      s       path seeded                  right hand side
  !      2       1! delta_1                   -B^(2) / 2!
  !      3       1! delta_1, 2! delta_2       -B^(3) / 3!
  !      ...
  !===================================================================!

  subroutine halley_correction(this, inputs, delta, achieved)

    class(newton)     , intent(inout) :: this
    type(stored_field), intent(in)    :: inputs(:)
    real(dp)          , intent(inout) :: delta(:)
    real(dp)          , intent(out)   :: achieved

    type(total_derivative)    :: total
    type(argument_path)       :: path
    class(field), allocatable :: out
    type(stored_field)        :: seeded
    real(dp)    , allocatable :: b(:), correction(:), individual(:,:)
    real(dp)                  :: fact, one_achieved
    integer                   :: s, p
    type(solve_result)        :: outcome

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
       ! factorial's value
       fact = fact * real(s - 1, dp)
       seeded = stored_field('correction', this % unknown_domain, size(delta))
       call seeded % set_real_vector(fact * individual(:, s - 1))
       path % derivative(s - 1) % occupied  = .true.
       path % derivative(s - 1) % direction = seeded

       ! The derivative of the statement, of this order, along this
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

       ! The series is valid only while it is decreasing: a correction
       ! no smaller than the one before it indicates the local model has
       ! left the regime a truncated expansion describes, and adding it
       ! would perturb rather than refine. What was already accumulated
       ! is retained; at s = 2 that is delta_1, Newton's own step, so
       ! this is never larger than Newton. The norm is the family's,
       ! reduced over the images, so every image takes the same
       ! decision.
       if (this % norm(correction) >= this % norm(individual(:, s - 1))) return

       individual(:, s) = correction
       delta = delta + correction

    end do

  end subroutine halley_correction

  !===================================================================!
  ! Binds the account of this minimizer and of its children
  !===================================================================!

  subroutine newton_bind_account(this, account)

    class(newton)       , intent(inout) :: this
    type(tally), pointer, intent(in)    :: account

    this % account => account
    if (allocated(this % inner)) call this % inner % bind_account(account)

  end subroutine newton_bind_account

  !===================================================================!
  ! Creates a newton over its inner minimizer, to the tolerance, by the
  ! chord method where chord
  !===================================================================!

  function newton_of(inner, tolerance, chord) result(this)

    class(minimizer), intent(in)           :: inner
    real(dp)        , intent(in), optional :: tolerance
    logical         , intent(in), optional :: chord
    type(newton)                           :: this

    allocate(this % inner, source=inner)
    if (present(tolerance)) this % tolerance = tolerance
    if (present(chord))     this % chord     = chord

  end function newton_of

end module operation_newton
