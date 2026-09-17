!=====================================================================!
! LEVEL 7 OF THE NEW TOWER . THE MINIMIZATION
!
! The first level with an objective, stated once: given a
! residual map R : U -> Y, vary the values on the UNKNOWN domain U
! to drive the values on the RESIDUAL domain Y toward zero. U and Y
! are set graph identities - never assumed to be the vertices of any
! graph; the graph argument remains only as the legacy
! operation host the compatibility apply() signature still requires. This
! module stores the minimizer base - ONE family for one purpose:
! state a statement, drive its residual to zero. Linear solvers,
! newton, and every other operation that minimizes a residual are its
! concretions; their differences are governance inside the family,
! never a second taxonomy. The solver vocabulary is defined here as
! thin delegations to engine entries, so the engine defines no
! solver term and the solver never calls apply or measure directly:
!
!      matvec ········· the operation applied, minus its constant
!      inner_product ·· a sum reduction with the second field as
!                       the measure
!      norm ··········· the norm reduction
!      sweep_order ···· the colouring traversal
!      block_diagonal · matvec evaluated by colour: applied to the
!                       indicator of one colour, the result at a
!                       member is that member's block column, because
!                       no two coupled blocks share a colour
!
! A concrete solver works in plain arrays - read once, updated,
! written back once, as the field banner specifies - and defines only
! its iteration. Everything else is inherited from here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_minimization

  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use iso_fortran_env , only : int64
  use util_precision  , only : dp, half_digits
  use util_norm, only : euclidean_norm
  use operation_action  , only : operation, contract
  use operation_action, only : binding, bound_value
  use operation_action, only : emit
  use view_directed   , only : directed_graph
  use field_calculus  , only : field, FIELD_REAL
  use graph_fractal      , only : graph
  use field_calculus        , only : functional
  use field_stored     , only : stored_field, typed_field_domain
  use operation_reduction , only : reduction, REDUCE_SUM, REDUCE_NORM
  use operation_traversal      , only : traversal, TRAVERSAL_COLOURING
  use operation_stencil        , only : stencil
  use util_tally               , only : tally

  implicit none

  private

  !-------------------------------------------------------------------!
  ! What a tolerance is measured against, and where an iteration that
  ! has not converged stops.
  !-------------------------------------------------------------------!

  integer, parameter, public :: relative = 1
  integer, parameter, public :: absolute = 2

  integer, parameter, public :: by_count = 1
  integer, parameter, public :: by_rate  = 2

  integer, parameter, public :: SOLVE_NOT_STARTED = -1, SOLVE_CONTINUE = 0
  integer, parameter, public :: SOLVE_CONVERGED = 1, SOLVE_EXHAUSTED = 2
  integer, parameter, public :: SOLVE_STAGNATED = 3, SOLVE_NONFINITE = 4
  integer, parameter, public :: SOLVE_DIVERGED = 5, SOLVE_SINGULAR = 6
  integer, parameter, public :: SOLVE_BREAKDOWN = 7, SOLVE_INNER_FAILED = 8
  integer, parameter, public :: SOLVE_EVALUATED = 9
  integer, parameter, public :: SOLVE_STORAGE_EXCEEDED = 10

  type, public :: solve_result
     real(dp) :: residual = 0.0_dp, initial_residual = 0.0_dp
     integer :: iterations = 0, reason = SOLVE_NOT_STARTED
   contains
     procedure :: converged => result_converged
     procedure :: failed => result_failed
     procedure :: description => result_description
  end type solve_result

  ! How many of the last imbalances a rate is fitted over. Three would
  ! leave one degree of freedom for the scatter and no more, so five is
  ! the smallest window that gives the slope a usable standard error.
  integer, parameter :: window = 5
  public :: minimizer
  public :: state
  public :: restrict
  public :: compact_labels
  public :: state_tuple
  public :: saturated_sum

  !===================================================================!
  ! The base: a stated operation, the graph it reads, and the
  ! tolerances every iteration honours.
  !===================================================================!

  type, abstract, extends(operation) :: minimizer

     class(operation), allocatable :: action

     ! THE EXECUTION CONTEXT. The graph passed to the action when the
     ! action is applied, and nothing more. The graph is whatever the
     ! action needs to compute with - a mesh, a compatibility host, an
     ! interface - and the graph determines nothing of the solver's
     ! own structure.
     class(directed_graph)          , allocatable :: graph

     ! THE DEPENDENT-VARIABLE COUPLING. Which unknowns depend on which:
     ! the stencil the structural algorithms need, and the ONLY graph
     ! sweep_order is permitted to colour.
     !
     ! The coupling is OPTIONAL WHEN STATED AND HAS NO FALLBACK. A
     ! structure-free minimizer - gmres, conjugate gradient, newton -
     ! never reads the coupling and need never supply one. A structured
     ! minimizer - jacobi, gauss-seidel - receives the coupling from a
     ! caller that records which object owns the dependent axis, and
     ! stops the program if the coupling was not supplied.
     !
     ! `coupling := graph` is the error this component exists to prevent:
     ! the graph an action executes over does not determine which
     ! unknowns are coupled. Where the two are the same graph, the
     ! CALLER states so, at its own call site.
     class(directed_graph)          , allocatable :: coupling

     ! The unknown domain U: the domain of the solution, explicit
     ! when stated, identity preserved - never inferred from the host.
     type(graph) :: unknown_domain
     integer         :: num_unknowns = 0

     ! The residual domain Y: the codomain of the action, read from
     ! the action itself when stated.
     type(graph) :: residual_domain
     integer         :: num_residuals = 0

     ! How wide an entry is. One number per cell is the common case
     ! and the default; a state with several numbers per cell - a
     ! complex plane point, a species vector, a whole spatial field
     ! at one instant - declares the width when stated, and every
     ! procedure below then measures the entire vector instead of its
     ! first component.
     integer :: num_components = 1

     ! THE BLOCK WIDTH. Unknowns that come in consecutive blocks of
     ! this width - a point's components, for example - and are coupled
     ! within the block far more strongly than across it, are smoothed
     ! one block at a time: the coupling stated is then over blocks,
     ! one vertex each, and the indicator reads a block's square
     ! submatrix. One is the default case, each unknown its own block.
     integer :: block_width = 1

     real(dp), allocatable :: affine(:)

     ! THE FIXED INPUTS. Inputs fixed during this solve - a
     ! scheme's history, an action's parameters - applied after the
     ! unknown on every evaluation. Solve-context data: the solver
     ! varies the unknown and nothing else.
     type(stored_field), allocatable :: stored(:)

     integer  :: max_iterations = 1000
     real(dp) :: tolerance      = half_digits

     ! WHAT THE TOLERANCE IS MEASURED AGAINST. relative divides the
     ! imbalance by the imbalance the iteration began with, which is
     ! the criterion whenever the target is a reduction. absolute
     ! compares the imbalance itself, which is the criterion only
     ! when the number has a meaning of its own - a functional driven
     ! under a stated value, and nothing else.
     integer :: criterion = relative

     ! WHERE THE ITERATION STOPS WHEN IT HAS NOT CONVERGED. by_count
     ! stops at max_iterations. by_rate also stops where the imbalance
     ! has flattened, which is where the slope of its logarithm over
     ! the last few iterations is no longer distinguishable from zero
     ! against its own scatter.
     integer :: limit_kind = by_count

     ! The imbalance the iteration began at, and the last few
     ! recorded. Updated only by the residual-history procedures.
     real(dp), private :: initial_residual = 0.0_dp
     real(dp), private :: residual_history(window) = 0.0_dp
     integer , private :: num_residual_samples = 0

     ! Whether any window has yet had a slope significantly below
     ! zero. An iteration that has never descended has not flattened
     ! either, whatever the slope over a window of its early iterations.
     logical , private :: residual_decreased = .false.

     ! false whenever the operator is stated: a cached block diagonal
     ! is valid only for the operator it was evaluated from
     logical :: diagonal_valid = .false.

     type(solve_result), private :: final_result

     ! THE ACCOUNT. The tally of the execution whose solve this is,
     ! bound by the caller for the duration of one solve and null
     ! between solves, so a stored minimizer copied with an execution
     ! references no other execution's tally.
     type(tally), pointer :: account => null()

   contains

     procedure :: bind_account
     procedure :: record_event

     procedure :: initialize_residual_history
     procedure :: record_residual_norm
     procedure :: converged
     procedure :: stagnated
     procedure :: diverging
     procedure :: initial_residual_norm
     procedure, private :: fit_log_residual
     procedure :: exhausted
     procedure :: terminated
     procedure :: record_result
     procedure :: result => minimizer_result

     procedure :: state
     procedure :: restrict
     procedure :: storage_entries
     procedure :: evaluate
     procedure :: matvec
     procedure :: imbalance
     procedure :: inner_product
     procedure :: norm
     procedure :: sweep_order
     procedure :: block_diagonal

     ! The operation interface: a solver IS an operation - the one
     ! that solves its own stated action. apply solves from zero, so a
     ! solver composes wherever operations compose; a preconditioner is
     ! exactly this interface of an inner solver.
     procedure :: domain => solver_domain
     procedure :: apply  => solver_apply

     procedure(solve_interface), deferred :: solve

  end type minimizer

  abstract interface

     !----------------------------------------------------------------!
     ! Drive || rhs - matvec(x) || under the tolerance, within the
     ! iteration limit. The achieved norm is reported whether or not
     ! the tolerance was met.
     !----------------------------------------------------------------!

     subroutine solve_interface(this, rhs, x, achieved)
       import :: minimizer, dp
       class(minimizer), intent(inout) :: this
       real(dp), intent(in)    :: rhs(:)
       real(dp), intent(inout) :: x(:)
       real(dp), intent(out)   :: achieved
     end subroutine solve_interface

  end interface

contains

  !===================================================================!
  ! Bind the tally the solves that follow record into, or null. A
  ! minimizer with children binds theirs as well.
  !===================================================================!

  subroutine bind_account(this, account)

    class(minimizer), intent(inout) :: this
    type(tally), pointer, intent(in) :: account

    this % account => account

  end subroutine bind_account

  !===================================================================!
  ! Record one event into the bound tally; nothing when unbound.
  !===================================================================!

  subroutine record_event(this, event)

    class(minimizer), intent(in) :: this
    integer         , intent(in) :: event

    if (associated(this % account)) call this % account % record(event)

  end subroutine record_event

  !===================================================================!
  ! The imbalance an iteration begins at, against which a relative
  ! tolerance is measured, and the last few recorded, from which a
  ! rate is fitted. Initialization clears the history,
  ! and the first recorded residual defines the relative tolerance.
  !===================================================================!

  subroutine initialize_residual_history(this)

    class(minimizer), intent(inout) :: this

    this % initial_residual     = 0.0_dp
    this % residual_history     = 0.0_dp
    this % num_residual_samples     = 0
    this % residual_decreased   = .false.
    this % final_result = solve_result()

  end subroutine initialize_residual_history

  subroutine record_residual_norm(this, imbalance)

    class(minimizer), intent(inout) :: this
    real(dp)        , intent(in)    :: imbalance

    real(dp) :: slope, error
    logical  :: regression_defined
    integer  :: i

    if (this % num_residual_samples == 0) then
       this % initial_residual = imbalance
    end if

    do i = 1, window - 1
       this % residual_history(i) = this % residual_history(i + 1)
    end do
    this % residual_history(window) = imbalance

    this % num_residual_samples = this % num_residual_samples + 1

    if (.not. all(ieee_is_finite(this % residual_history))) return
    call this % fit_log_residual(slope, error, regression_defined)
    if (regression_defined .and. slope < -error) this % residual_decreased   = .true.

  end subroutine record_residual_norm

  !===================================================================!
  ! Whether the imbalance meets the tolerance. Relative divides
  ! by the imbalance the iteration began at; absolute does not divide
  ! at all. A criterion that is neither stops the program.
  !===================================================================!

  logical function converged(this, imbalance) result(criterion_met)

    class(minimizer), intent(in) :: this
    real(dp)        , intent(in) :: imbalance
    integer :: relative_exponent
    character(len=250) :: message

    criterion_met = .false.
    if (.not. ieee_is_finite(imbalance) .or. .not. ieee_is_finite(this % initial_residual)) return
    if (imbalance < 0.0_dp) return
    if (.not. ieee_is_finite(this % tolerance)) return
    if (this % tolerance < 0.0_dp) return
    select case (this % criterion)
    case (relative)
       if (imbalance == 0.0_dp) then
          criterion_met = .true.
          return
       end if
       if (this % initial_residual <= 0.0_dp .or. this % tolerance == 0.0_dp) return
       ! Compare r <= tolerance*r0 without an underflowing product
       ! or an overflowing ratio. The fractions lie in [1/2, 1).
       relative_exponent = exponent(imbalance) - exponent(this % tolerance) - exponent(this % initial_residual)
       if (relative_exponent <= -2) then
          criterion_met = .true.
       else if (relative_exponent <= 0) then
          criterion_met = scale(fraction(imbalance), relative_exponent) <= &
               & fraction(this % tolerance) * fraction(this % initial_residual)
       end if
    case (absolute)
       criterion_met = imbalance <= this % tolerance
    case default
       write(message,'(a,i0)') 'minimizer: a tolerance must be measured relative or absolute; &
            &this % criterion = ', this % criterion
       error stop trim(message)
    end select

  end function converged

  !===================================================================!
  ! Whether the imbalance has stopped falling: the slope of its
  ! logarithm over the last window against the standard error of that
  ! slope. A descent still under way has a slope far outside its own
  ! error; scatter about a lower bound has a slope inside it. No
  ! constant is chosen here except the width of the window, and no scale enters,
  ! the slope of a logarithm being dimensionless.
  !===================================================================!

  logical function stagnated(this) result(is_stagnating)

    class(minimizer), intent(in) :: this

    real(dp) :: slope, error
    logical  :: regression_defined

    call this % fit_log_residual(slope, error, regression_defined)
    is_stagnating = regression_defined .and. this % residual_decreased .and. abs(slope) < error

  end function stagnated

  !===================================================================!
  ! The slope of the logarithm of the last window of imbalances, and
  ! the standard error of that slope. Not usable until the window is
  ! full, and not usable where an imbalance has reached zero, there
  ! being no logarithm of it.
  !===================================================================!

  subroutine fit_log_residual(this, slope, error, regression_defined)

    class(minimizer), intent(in)  :: this
    real(dp)        , intent(out) :: slope, error
    logical         , intent(out) :: regression_defined

    real(dp) :: x(window), y(window), mx, my, sxx, sxy, scatter
    integer  :: i

    slope  = 0.0_dp
    error  = 0.0_dp
    regression_defined = .false.

    if (this % num_residual_samples < window) return
    if (.not. all(ieee_is_finite(this % residual_history))) return

    do i = 1, window
       if (this % residual_history(i) <= 0.0_dp) return
       x(i) = real(i, dp)
       y(i) = log(this % residual_history(i))
    end do

    mx  = sum(x) / real(window, dp)
    my  = sum(y) / real(window, dp)
    sxx = sum((x - mx) ** 2)
    sxy = sum((x - mx) * (y - my))

    slope   = sxy / sxx
    scatter = sum((y - (my + slope * (x - mx))) ** 2)
    error   = sqrt(scatter / real(window - 2, dp) / sxx)

    regression_defined = .true.

  end subroutine fit_log_residual

  !===================================================================!
  ! Whether the imbalance is diverging: above the imbalance the
  ! iteration began at, growing, and at the least rate consistent with
  ! the window - the slope less its error - projected to exceed the
  ! largest representable number before the iteration limit is
  ! reached. An early rise that later descends is not divergence: its
  ! projection over the remaining iterations stays finite. The one
  ! limit named is the arithmetic's own.
  !===================================================================!

  logical function diverging(this, imbalance) result(satisfied)

    class(minimizer), intent(in) :: this
    real(dp)        , intent(in) :: imbalance

    real(dp) :: slope, error, remaining
    logical  :: regression_defined

    satisfied = .false.
    if (.not. ieee_is_finite(imbalance) .or. .not. ieee_is_finite(this % initial_residual)) return
    if (imbalance <= this % initial_residual .or. imbalance <= 0.0_dp) return

    call this % fit_log_residual(slope, error, regression_defined)
    if (.not. regression_defined .or. slope <= error) return

    remaining = real(max(this % max_iterations - this % num_residual_samples, 0), dp)
    satisfied = log(imbalance) + (slope - error) * remaining > log(huge(1.0_dp))

  end function diverging

  pure real(dp) function initial_residual_norm(this) result(imbalance)

    class(minimizer), intent(in) :: this

    imbalance = this % initial_residual

  end function initial_residual_norm

  !===================================================================!
  ! Whether the iteration has reached its limit. A limit setting
  ! that is neither by_count nor by_rate stops the program.
  !===================================================================!

  logical function exhausted(this, iteration) result(criterion_met)

    class(minimizer), intent(in) :: this
    integer         , intent(in) :: iteration
    character(len=250) :: message

    criterion_met = iteration >= this % max_iterations

    select case (this % limit_kind)
    case (by_count)
       continue
    case (by_rate)
       criterion_met = criterion_met .or. this % stagnated()
    case default
       write(message,'(a,i0)') 'minimizer: an iteration limit must be counted or taken from the &
            &rate; this % limit_kind = ', this % limit_kind
       error stop trim(message)
    end select

  end function exhausted

  !===================================================================!
  ! Whether an iteration stops at this imbalance: recorded, then
  ! tested - converged; or not a number; or above the largest
  ! representable number; or diverging; or the iteration limit
  ! reached. Every iteration in the hierarchy evaluates this one
  ! predicate on the imbalance just measured, and the loop around it
  ! belongs to the member.
  !===================================================================!

  logical function terminated(this, imbalance, iteration) result(is_terminated)

    class(minimizer), intent(inout) :: this
    real(dp)        , intent(in)    :: imbalance
    integer         , intent(in)    :: iteration

    call this % record_residual_norm(imbalance)
    call this % record_result(imbalance, iteration)
    is_terminated = this % final_result % reason /= SOLVE_CONTINUE

  end function terminated

  ! A result records the true residual and the number of completed
  ! iterations. Recording it does not add a second convergence sample.
  subroutine record_result(this, imbalance, iteration, reason)
    class(minimizer), intent(inout) :: this
    real(dp), intent(in) :: imbalance
    integer, intent(in) :: iteration
    integer, intent(in), optional :: reason
    character(len=250) :: message
    this % final_result % residual = imbalance
    this % final_result % initial_residual = this % initial_residual
    this % final_result % iterations = iteration
    if (present(reason)) then
       this % final_result % reason = reason
    else if (.not. ieee_is_finite(imbalance)) then
       this % final_result % reason = SOLVE_NONFINITE
    else if (this % converged(imbalance)) then
       this % final_result % reason = SOLVE_CONVERGED
    else if (imbalance > huge(1.0_dp) / 2.0_dp) then
       this % final_result % reason = SOLVE_DIVERGED
    else if (this % diverging(imbalance)) then
       this % final_result % reason = SOLVE_DIVERGED
    else if (iteration >= this % max_iterations) then
       this % final_result % reason = SOLVE_EXHAUSTED
    else if (this % limit_kind == by_rate .and. this % stagnated()) then
       this % final_result % reason = SOLVE_STAGNATED
    else
       this % final_result % reason = SOLVE_CONTINUE
    end if
    if (this % limit_kind /= by_count .and. this % limit_kind /= by_rate) then
       write(message,'(a,i0)') 'minimizer: an iteration limit must be counted or taken from the &
            &rate; this % limit_kind = ', this % limit_kind
       error stop trim(message)
    end if
  end subroutine record_result

  pure function minimizer_result(this) result(outcome)
    class(minimizer), intent(in) :: this
    type(solve_result) :: outcome
    outcome = this % final_result
  end function minimizer_result

  pure logical function result_converged(this) result(satisfied)
    class(solve_result), intent(in) :: this
    satisfied = this % reason == SOLVE_CONVERGED
  end function result_converged

  ! Exhaustion and stagnation can supply an approximate correction.
  ! These reasons instead report a numerical failure of that
  ! correction, or a construction refused within its storage limit,
  ! which supplies no correction at all.
  pure logical function result_failed(this) result(satisfied)
    class(solve_result), intent(in) :: this
    satisfied = (this % reason >= SOLVE_NONFINITE .and. this % reason <= SOLVE_INNER_FAILED) &
         & .or. this % reason == SOLVE_STORAGE_EXCEEDED
  end function result_failed

  pure function result_description(this) result(description)
    class(solve_result), intent(in) :: this
    character(len=:), allocatable :: description
    select case (this % reason)
    case (SOLVE_NOT_STARTED); description = 'not started'
    case (SOLVE_CONTINUE); description = 'iteration in progress'
    case (SOLVE_CONVERGED); description = 'converged'
    case (SOLVE_EXHAUSTED); description = 'iteration limit reached'
    case (SOLVE_STAGNATED); description = 'residual stagnated'
    case (SOLVE_NONFINITE); description = 'nonfinite residual'
    case (SOLVE_DIVERGED); description = 'residual diverged'
    case (SOLVE_SINGULAR); description = 'singular matrix'
    case (SOLVE_BREAKDOWN); description = 'linear iteration breakdown'
    case (SOLVE_INNER_FAILED); description = 'inner solve failed'
    case (SOLVE_EVALUATED); description = 'schedule evaluated'
    case (SOLVE_STORAGE_EXCEEDED); description = 'storage limit exceeded'
    case default; description = 'invalid solve result'
    end select
  end function result_description

  !===================================================================!
  ! Store the operation and the graph it reads. The affine part is
  ! evaluated here, once: the operation applied to the zero state is
  ! the contribution of its boundary values and sources alone.
  !===================================================================!

  subroutine state(this, action, context, unknown_domain, num_unknowns, &
       & num_components, coupling, stored_inputs)

    class(minimizer)  , intent(inout) :: this
    class(operation), intent(in)    :: action
    class(directed_graph)          , intent(in)    :: context
    type(graph)       , intent(in)    :: unknown_domain
    integer               , intent(in)    :: num_unknowns
    integer, intent(in), optional         :: num_components
    class(directed_graph)  , intent(in), optional  :: coupling
    type(stored_field), intent(in), optional :: stored_inputs(:)

    real(dp), allocatable :: zero(:)
    integer :: n
    character(len=250) :: message

    call this % initialize_residual_history()
    if (allocated(this % action)) deallocate(this % action)
    allocate(this % action, source=action)
    if (allocated(this % graph)) deallocate(this % graph)
    allocate(this % graph, source=context)

    ! the inputs fixed during this solve follow the unknown in
    ! every evaluation, the affine part's included
    if (allocated(this % stored)) deallocate(this % stored)
    if (present(stored_inputs)) allocate(this % stored, source=stored_inputs)

    ! The dependent-variable coupling arrives EXPLICIT or not at all.
    ! No fallback to the execution context: a solver that needs
    ! structure and was given none stops the program when it reads
    ! the coupling, rather than colouring the execution context.
    if (allocated(this % coupling)) deallocate(this % coupling)
    if (present(coupling)) allocate(this % coupling, source=coupling)

    this % num_components = 1
    if (present(num_components)) this % num_components = max(num_components, 1)

    ! The solver reads one right-hand side with the configured width.
    call this % declare_arguments(1, [contract(FIELD_REAL, this % num_components)])

    ! The unknown domain arrives EXPLICIT and identity-preserving.
    ! No hidden fallback to the host's vertices: a caller that
    ! means vertices passes them at its own call site.
    this % unknown_domain   = unknown_domain
    this % num_unknowns = num_unknowns

    ! The residual domain is the action's own codomain.
    call action % domain(context, this % residual_domain, this % num_residuals)

    n = this % num_unknowns

    ! THE CONSTANT PART IS THE STATEMENT AT THE ZERO STATE. A(0) is
    ! the contribution of boundary values and sources alone when
    ! A(x) = A(0) + L x, and matvec and the operations built on it -
    ! the operations of a linear solve - are its only readers. A
    ! statement whose domain excludes the zero state, as the radial
    ! oscillator's nu / q**3 does, is not evaluated there and declares
    ! no constant part.
    if (allocated(this % affine)) deallocate(this % affine)
    if (action % defined_at_zero()) then
       allocate(zero(n * this % num_components))
       zero = 0.0_dp
       call evaluate(this, zero, this % affine)

       ! Classification is not admissibility: U and Y stay distinct
       ! identities, but THIS solver family is square - the scalar
       ! dimensions must agree. A rectangular least-squares family may
       ! support R^n -> R^m later; none exists yet.
       if (size(this % affine) /= n * this % num_components) then
          write(message,'(a,i0,a,i0)') 'minimization: the current solver family requires equal &
               &unknown and residual value dimensions; size(affine) = ', size(this % affine), &
               & ', n * num_components = ', n * this % num_components
          error stop trim(message)
       end if
    else
       allocate(this % affine(0))
    end if

    this % diagonal_valid = .false.

  end subroutine state

  !===================================================================!
  ! RESTRICTION TO A SELECTED DOMAIN. selected(i) is the whole index
  ! of the i-th unknown of the selected domain: an injective map of
  ! the selected domain into the whole. The family is square, so the
  ! selected residual rows are the rows of the same indices.
  !
  ! The base restricts what every minimizer owns. The block layout is
  ! preserved: a selection is whole blocks of the stated width, aligned
  ! as the whole numbers them, in block order. Every quantity evaluated
  ! on the whole domain - the statement, its affine part, the coupling,
  ! the stored inputs, the block diagonal, the residual history - is
  ! invalid on the selected domain and is discarded. The operator on
  ! the selected domain, the whole constrained to the selection with
  ! its exterior fixed, is stated afterwards; its affine part R_s(0)
  ! is then the exterior contribution. Stopping rules, the block
  ! width and the component width are the same on the selection.
  !
  ! A minimizer owning metadata indexed by the whole, or children over
  ! domains derived from it, overrides this: it calls this base, maps
  ! its metadata through the selection, and restricts each child by the
  ! selection induced on the child's domain.
  !
  ! Invalid input: an empty selection, an index below one or beyond
  ! the stated unknowns, a repeated index, or a selection that splits
  ! or misaligns a block.
  !===================================================================!

  subroutine restrict(this, selected)

    class(minimizer), intent(inout) :: this
    integer         , intent(in)    :: selected(:)

    logical, allocatable :: chosen(:)
    integer :: i, b, w, nb
    character(len=250) :: message

    if (size(selected) < 1) then
       error stop 'minimization: a restriction must select an unknown at least; selected is empty'
    end if
    if (any(selected < 1)) then
       write(message,'(a,i0)') 'minimization: every selected index must belong to the whole &
            &domain; minval(selected) = ', minval(selected)
       error stop trim(message)
    end if
    if (this % num_unknowns > 0) then
       if (any(selected > this % num_unknowns)) then
          write(message,'(a,i0,a,i0)') 'minimization: every selected index must belong to the &
               &whole domain; maxval(selected) = ', maxval(selected), ', num_unknowns = ', &
               & this % num_unknowns
          error stop trim(message)
       end if
    end if
    allocate(chosen(maxval(selected)), source=.false.)
    do i = 1, size(selected)
       if (chosen(selected(i))) then
          write(message,'(a,i0,a,i0)') 'minimization: a restriction must select each unknown &
               &once; selected(', i, ') repeats the index ', selected(i)
          error stop trim(message)
       end if
       chosen(selected(i)) = .true.
    end do

    w = this % block_width
    if (w > 1) then
       nb = size(selected) / w
       if (nb * w /= size(selected)) then
          write(message,'(a,i0,a,i0)') 'minimization: a restriction must select whole blocks; &
               &size(selected) = ', size(selected), ' is not a multiple of block_width = ', w
          error stop trim(message)
       end if
       do b = 1, nb
          if (mod(selected((b - 1) * w + 1) - 1, w) /= 0) then
             write(message,'(a,i0,a,i0)') 'minimization: a restriction must select whole blocks; &
                  &block ', b, ' starts at unaligned index ', selected((b - 1) * w + 1)
             error stop trim(message)
          end if
          do i = 2, w
             if (selected((b - 1) * w + i) /= selected((b - 1) * w + 1) + i - 1) then
                write(message,'(a,i0,a,i0,a,i0)') 'minimization: a restriction must select whole &
                     &blocks; block ', b, ' entry ', i, ' is not contiguous with the block''s &
                     &first index ', selected((b - 1) * w + 1)
                error stop trim(message)
             end if
          end do
       end do
    end if

    if (allocated(this % action))   deallocate(this % action)
    if (allocated(this % graph))    deallocate(this % graph)
    if (allocated(this % coupling)) deallocate(this % coupling)
    if (allocated(this % stored))   deallocate(this % stored)
    if (allocated(this % affine))   deallocate(this % affine)
    this % num_unknowns   = 0
    this % num_residuals  = 0
    this % diagonal_valid = .false.
    call this % initialize_residual_history()

  end subroutine restrict

  !===================================================================!
  ! The entries a minimizer retains for a statement over num_unknowns
  ! unknowns beyond a fixed number of vectors over them: a
  ! factorisation, a Krylov basis, a block diagonal, a complement, and
  ! its children's. An enclosing minimizer with a storage limit reads
  ! this before stating its inner minimizer, so a requirement beyond
  ! the limit is refused before the allocation. The base retains none.
  !===================================================================!

  pure integer(int64) function storage_entries(this, num_unknowns) result(entries)

    class(minimizer), intent(in) :: this
    integer         , intent(in) :: num_unknowns

    associate (u1 => this, u2 => num_unknowns); end associate
    entries = 0_int64

  end function storage_entries

  !===================================================================!
  ! a + b in 64-bit integers, the largest representable value where
  ! the sum exceeds it: a storage requirement beyond every limit is
  ! reported as such, never as a wrapped count.
  !===================================================================!

  pure integer(int64) function saturated_sum(a, b) result(total)

    integer(int64), intent(in) :: a, b

    if (a < 0_int64 .or. b < 0_int64) then
       error stop 'minimization: a storage requirement must be a nonnegative count, but a &
            &negative term was passed to saturated_sum'
    end if
    if (a > huge(a) - b) then
       total = huge(a)
    else
       total = a + b
    end if

  end function saturated_sum

  !===================================================================!
  ! Labels renumbered compactly in order of first appearance:
  ! mapped(i) is the compact label of label(i), and representative(j)
  ! the original label of the compact label j. A metadata
  ! restriction maps its labels through this, and the representatives
  ! are the selection induced on the domain the labels index.
  !===================================================================!

  subroutine compact_labels(label, mapped, representative)

    integer, intent(in) :: label(:)
    integer, allocatable, intent(out) :: mapped(:), representative(:)

    integer, allocatable :: distinct(:)
    integer :: i, at, n

    allocate(mapped(size(label)), distinct(size(label)))
    n = 0
    do i = 1, size(label)
       at = 0
       if (n > 0) at = findloc(distinct(1:n), label(i), dim=1)
       if (at == 0) then
          n = n + 1
          distinct(n) = label(i)
          at = n
       end if
       mapped(i) = at
    end do
    representative = distinct(1:n)

  end subroutine compact_labels

  !===================================================================!
  ! The inputs a statement is evaluated on at the state x: x stored
  ! on the unknown domain, n members of the given width, then the
  ! fixed inputs where given. The residual and every tangent taken of
  ! it are built on this one tuple, so they linearize the same function.
  !===================================================================!

  function state_tuple(domain, n, components, x, stored) result(inputs)

    type(graph)       , intent(in)           :: domain
    integer           , intent(in)           :: n, components
    real(dp)          , intent(in)           :: x(:)
    type(stored_field), intent(in), optional :: stored(:)
    type(stored_field), allocatable          :: inputs(:)

    type(stored_field) :: state
    type(typed_field_domain) :: unknowns

    unknowns = typed_field_domain(domain, n, components)
    state    = unknowns % state(x)

    if (present(stored)) then
       inputs = [state, stored]
    else
       inputs = [state]
    end if

  end function state_tuple

  !===================================================================!
  ! The operation applied unmodified, affine part included. The
  ! input tuple the operation was applied on is returned when
  ! requested, so that a tangent frozen on that tuple linearizes the
  ! function that was evaluated.
  !===================================================================!

  subroutine evaluate(this, x, y, inputs)

    class(minimizer), intent(in)   :: this
    real(dp), intent(in)               :: x(:)
    real(dp), allocatable, intent(out) :: y(:)
    type(stored_field), allocatable, intent(out), optional :: inputs(:)

    type(stored_field), allocatable :: tuple(:)
    class(field), allocatable :: image
    character(len=250) :: message

    if (.not. allocated(this % action)) then
       error stop 'minimization: evaluate was called before the operator was stated on the &
            &solver domain; this % action is not allocated'
    end if
    tuple = state_tuple(this % unknown_domain, this % num_unknowns, this % num_components, x, this % stored)
    call this % action % apply(this % graph, this % action % bind(tuple), image)

    if (.not. image % defined_on(this % residual_domain)) then
       error stop 'minimization: the action returned a field that is not defined on the &
            &statement''s residual domain'
    end if

    call image % real_vector(y)

    ! A statement with no constant part was not evaluated when it was
    ! stated, so the family's squareness is checked here instead.
    if (size(this % affine) == 0 .and. this % num_unknowns > 0) then
       if (size(y) /= this % num_unknowns * this % num_components) then
          write(message,'(a,i0,a,i0)') 'minimization: the current solver family requires equal &
               &unknown and residual value dimensions; size(y) = ', size(y), ', num_unknowns * &
               &num_components = ', this % num_unknowns * this % num_components
          error stop trim(message)
       end if
    end if

    if (present(inputs)) call move_alloc(tuple, inputs)

  end subroutine evaluate

  !===================================================================!
  ! The solver's operations, each a delegation.
  !===================================================================!

  subroutine matvec(this, x, y)

    class(minimizer), intent(in)   :: this
    real(dp), intent(in)               :: x(:)
    real(dp), allocatable, intent(out) :: y(:)
    character(len=250) :: message

    call evaluate(this, x, y)

    ! A x is the statement less its value at the zero state, which a
    ! statement whose domain excludes that state does not have.
    if (size(this % affine) /= size(y)) then
       write(message,'(a,i0,a,i0)') 'minimization: a matrix-vector product subtracts the &
            &statement at the zero state, which lies outside this statement''s domain; &
            &size(affine) = ', size(this % affine), ', size(y) = ', size(y)
       error stop trim(message)
    end if

    y = y - this % affine

  end subroutine matvec

  !===================================================================!
  ! The linear solver imbalance: rhs - A x, where A is the action
  ! with its affine part removed by matvec. Iterative solvers call
  ! this one procedure instead of each assembling the residual.
  !===================================================================!

  subroutine imbalance(this, rhs, x, r)

    class(minimizer), intent(in)   :: this
    real(dp), intent(in)               :: rhs(:), x(:)
    real(dp), allocatable, intent(out) :: r(:)

    call this % matvec(x, r)
    r = rhs - r

  end subroutine imbalance

  real(dp) function inner_product(this, u, v) result(prod)

    class(minimizer), intent(in) :: this
    real(dp), intent(in) :: u(:), v(:)

    ! a sum reduction of u weighted by v is the sum of the products,
    ! taken here without a field allocated to store one number
    prod = sum(u * v)

  end function inner_product

  real(dp) function norm(this, u) result(length)

    class(minimizer), intent(in) :: this
    real(dp), intent(in) :: u(:)

    length = euclidean_norm(u)

  end function norm

  subroutine sweep_order(this, colours)

    class(minimizer), intent(in)  :: this
    integer, allocatable, intent(out) :: colours(:)

    type(traversal) :: colouring
    class(field), allocatable :: image

    ! THE COLOURING IS OF THE UNKNOWNS' COUPLING, never of the
    ! execution context. Two unknowns may share a colour only when
    ! nothing couples them, and the graph an action executes over does
    ! not record that.
    if (.not. allocated(this % coupling)) then
       error stop 'minimization: sweep_order requires the dependent-variable coupling, but this &
            &% coupling is not allocated - state it with coupling='
    end if

    colouring = traversal(TRAVERSAL_COLOURING)
    call colouring % apply(this % coupling, output=image)
    call image % integer_vector(colours)

  end subroutine sweep_order

  !===================================================================!
  ! THE BLOCK DIAGONAL by coloured indicators. The coupling stated is over
  ! the blocks; blocks of one colour do not couple, so an indicator of
  ! one on the k-th component of every block of a colour, applied
  ! through the matvec, returns the k-th column of every one of those
  ! blocks' square submatrices at once. width indicators per colour read the whole
  ! block diagonal, and with a width of one this is the diagonal.
  !===================================================================!

  subroutine block_diagonal(this, d)

    class(minimizer), intent(in)   :: this
    real(dp), allocatable, intent(out) :: d(:,:,:)

    integer , allocatable :: colours(:)
    real(dp), allocatable :: indicator(:), y(:)
    integer :: n, w, nb, col, b, k, i
    character(len=250) :: message

    if (this % num_components > 1) then
       write(message,'(a,i0)') 'diagonal: the coloured indicator must evaluate one number per cell; &
            &num_components = ', this % num_components
       error stop trim(message)
    end if

    n  = this % num_unknowns * this % num_components
    w  = this % block_width
    nb = n / w
    if (nb * w /= n) then
       write(message,'(a,i0,a,i0)') 'diagonal: the unknowns must come in whole blocks; n = ', n, &
            & ', block_width = ', w
       error stop trim(message)
    end if

    ! an explicit operator states its diagonal blocks on its edges;
    ! any other is read by coloured indicator products
    select type (a => this % action)
    type is (stencil)
       call a % diagonal_blocks(w, d)
       return
    end select

    allocate(d(w, w, nb), indicator(n))
    d = 0.0_dp

    call this % sweep_order(colours)
    if (size(colours) /= nb) then
       write(message,'(a,i0,a,i0)') 'diagonal: the coupling must be stated over the blocks, one &
            &colour each; size(colours) = ', size(colours), ', nb = ', nb
       error stop trim(message)
    end if

    do col = 1, maxval(colours)
       do k = 1, w
          indicator = 0.0_dp
          do b = 1, nb
             if (colours(b) == col) indicator((b - 1) * w + k) = 1.0_dp
          end do
          call this % matvec(indicator, y)
          do b = 1, nb
             if (colours(b) /= col) cycle
             do i = 1, w
                d(i, k, b) = y((b - 1) * w + i)
             end do
          end do
       end do
    end do

  end subroutine block_diagonal

  !===================================================================!
  ! The operation interface.
  !===================================================================!

  subroutine solver_domain(this, input_graph, domain, num_entries)

    class(minimizer), intent(in)       :: this
    class(directed_graph), intent(in)               :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries

    associate (u1 => input_graph); end associate

    ! The solver's output is a solution on U.
    domain   = this % unknown_domain
    num_entries = this % num_unknowns

  end subroutine solver_domain

  subroutine solver_apply(this, input_graph, inputs, output)

    class(minimizer), intent(in)                   :: this
    class(directed_graph), intent(in)                       :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    class(minimizer), allocatable :: copy
    class(field), allocatable :: right_hand_side
    type(stored_field) :: out
    type(typed_field_domain) :: unknowns
    real(dp), allocatable :: rhs(:), x(:)
    real(dp) :: achieved

    associate (u1 => input_graph); end associate

    ! x IS a state on the unknown domain, allocated at that extent.
    allocate(x(this % num_unknowns * this % num_components))
    x = 0.0_dp

    if (present(inputs)) then
       call bound_value(inputs, this % argument(1), right_hand_side)
       if (.not. right_hand_side % defined_on(this % residual_domain)) then
          error stop 'minimization: the bound right-hand side is not defined on the stated &
               &residual domain'
       end if
       call right_hand_side % real_vector(rhs)
       allocate(copy, source=this)
       call copy % solve(rhs, x, achieved)
    end if

    unknowns = typed_field_domain(this % unknown_domain, this % num_unknowns, this % num_components)
    out      = unknowns % solution(x)

    call emit(out, output)

  end subroutine solver_apply

end module operation_minimization
