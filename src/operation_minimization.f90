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
! attach a statement, drive its residual to zero. Linear solvers,
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

  use util_precision  , only : dp, half_digits
  use operation_action  , only : operation, contract
  use operation_action, only : binding, bound_value
  use operation_action, only : emit
  use view_directed   , only : directed_graph
  use field_calculus  , only : field, FIELD_REAL
  use graph_fractal      , only : graph
  use field_calculus        , only : functional
  use field_stored     , only : stored_field
  use operation_reduction , only : reduction, REDUCE_SUM, REDUCE_NORM
  use operation_traversal      , only : traversal, TRAVERSAL_COLOURING

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

  ! How many of the last imbalances a rate is fitted over. Three would
  ! leave one degree of freedom for the scatter and no more, so five is
  ! the smallest window that gives the slope a usable standard error.
  integer, parameter :: window = 5
  public :: minimizer
  public :: attach
  public :: state_tuple

  !===================================================================!
  ! The base: an attached operation, the graph it reads, and the
  ! tolerances every iteration honours.
  !===================================================================!

  type, abstract, extends(operation) :: minimizer

     class(operation), allocatable :: action

     ! THE EXECUTION CONTEXT. The graph passed to the action when the
     ! action is applied, and nothing more. The graph is whatever the
     ! action needs to compute with - a mesh, a compatibility host, an
     ! interface - and the graph determines nothing of the solver's
     ! own structure.
     class(directed_graph)          , allocatable :: on

     ! THE DEPENDENT-VARIABLE COUPLING. Which unknowns depend on which:
     ! the stencil the structural algorithms need, and the ONLY graph
     ! sweep_order is permitted to colour.
     !
     ! The coupling is OPTIONAL AT ATTACH AND HAS NO FALLBACK. A
     ! structure-free minimizer - gmres, conjugate gradient, newton -
     ! never reads the coupling and need never supply one. A structured
     ! minimizer - jacobi, gauss-seidel - receives the coupling from a
     ! caller that records which object owns the dependent axis, and
     ! stops the program if the coupling was not supplied.
     !
     ! `coupling := on` is the error this component exists to prevent:
     ! the graph an action executes over does not determine which
     ! unknowns are coupled. Where the two are the same graph, the
     ! CALLER states so, at its own call site.
     class(directed_graph)          , allocatable :: coupling

     ! The unknown domain U: the domain of the solution, explicit at
     ! attach, identity preserved - never inferred from the host.
     type(graph) :: unknown_domain
     integer         :: num_unknowns = 0

     ! The residual domain Y: the codomain of the action, read from
     ! the action itself at attach.
     type(graph) :: residual_domain
     integer         :: num_residuals = 0

     ! How wide an entry is. One number per cell is the common case
     ! and the default; a state with several numbers per cell - a
     ! complex plane point, a species vector, a whole spatial field
     ! at one instant - declares the width at attach, and every
     ! procedure below then measures the entire vector instead of its
     ! first component.
     integer :: num_components = 1

     ! THE BLOCK WIDTH. Unknowns that come in consecutive blocks of
     ! this width - a point's components, for example - and are coupled
     ! within the block far more strongly than across it, are smoothed
     ! one block at a time: the coupling attached is then over blocks,
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
     ! recorded. Written by note_imbalance and by nothing else.
     real(dp), private :: began_at = 0.0_dp
     real(dp), private :: recent(window) = 0.0_dp
     integer , private :: noted = 0

     ! Whether any window has yet had a slope significantly below
     ! zero. An iteration that has never descended has not flattened
     ! either, whatever the slope over a window of its early iterations.
     logical , private :: descended = .false.

     ! false whenever the operator is attached: a cached block diagonal
     ! is valid only for the operator it was evaluated from
     logical :: diagonal_valid = .false.

   contains

     procedure :: begin_imbalance
     procedure :: note_imbalance
     procedure :: converged
     procedure :: flattened
     procedure :: diverging
     procedure :: began
     procedure, private :: fitted
     procedure :: exhausted
     procedure :: halted

     procedure :: attach
     procedure :: raw_apply
     procedure :: matvec
     procedure :: imbalance
     procedure :: inner_product
     procedure :: norm
     procedure :: sweep_order
     procedure :: block_diagonal

     ! The operation interface: a solver IS an operation - the one
     ! that solves the attached statement. apply solves from zero, so a
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
  ! The imbalance an iteration begins at, against which a relative
  ! tolerance is measured, and the last few recorded, from which a
  ! rate is fitted. Written here and nowhere else: begin clears them,
  ! and the first imbalance noted is the one the iteration began at.
  !===================================================================!

  subroutine begin_imbalance(this)

    class(minimizer), intent(inout) :: this

    this % began_at  = 0.0_dp
    this % recent    = 0.0_dp
    this % noted     = 0
    this % descended = .false.

  end subroutine begin_imbalance

  subroutine note_imbalance(this, imbalance)

    class(minimizer), intent(inout) :: this
    real(dp)        , intent(in)    :: imbalance

    real(dp) :: slope, error
    logical  :: usable
    integer  :: i

    if (this % noted == 0 .and. this % began_at <= 0.0_dp) then
       this % began_at = imbalance
    end if

    do i = 1, window - 1
       this % recent(i) = this % recent(i + 1)
    end do
    this % recent(window) = imbalance

    this % noted = this % noted + 1

    call this % fitted(slope, error, usable)
    if (usable .and. slope < -error) this % descended = .true.

  end subroutine note_imbalance

  !===================================================================!
  ! Whether the imbalance meets the tolerance. Relative divides
  ! by the imbalance the iteration began at; absolute does not divide
  ! at all. A criterion that is neither stops the program.
  !===================================================================!

  logical function converged(this, imbalance) result(done)

    class(minimizer), intent(in) :: this
    real(dp)        , intent(in) :: imbalance

    select case (this % criterion)
    case (relative)
       done = imbalance <= this % tolerance * max(this % began_at, tiny(1.0_dp))
    case (absolute)
       done = imbalance <= this % tolerance
    case default
       error stop 'minimizer: a tolerance is measured relative or absolute'
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

  logical function flattened(this) result(flat)

    class(minimizer), intent(in) :: this

    real(dp) :: slope, error
    logical  :: usable

    call this % fitted(slope, error, usable)
    flat = usable .and. this % descended .and. abs(slope) < error

  end function flattened

  !===================================================================!
  ! The slope of the logarithm of the last window of imbalances, and
  ! the standard error of that slope. Not usable until the window is
  ! full, and not usable where an imbalance has reached zero, there
  ! being no logarithm of it.
  !===================================================================!

  subroutine fitted(this, slope, error, usable)

    class(minimizer), intent(in)  :: this
    real(dp)        , intent(out) :: slope, error
    logical         , intent(out) :: usable

    real(dp) :: x(window), y(window), mx, my, sxx, sxy, scatter
    integer  :: i

    slope  = 0.0_dp
    error  = 0.0_dp
    usable = .false.

    if (this % noted < window) return

    do i = 1, window
       if (this % recent(i) <= 0.0_dp) return
       x(i) = real(i, dp)
       y(i) = log(this % recent(i))
    end do

    mx  = sum(x) / real(window, dp)
    my  = sum(y) / real(window, dp)
    sxx = sum((x - mx) ** 2)
    sxy = sum((x - mx) * (y - my))

    slope   = sxy / sxx
    scatter = sum((y - (my + slope * (x - mx))) ** 2)
    error   = sqrt(scatter / real(window - 2, dp) / sxx)

    usable = .true.

  end subroutine fitted

  !===================================================================!
  ! Whether the imbalance is diverging: above the imbalance the
  ! iteration began at, growing, and at the least rate consistent with
  ! the window - the slope less its error - projected to exceed the
  ! largest representable number before the iteration limit is
  ! reached. An early rise that later descends is not divergence: its
  ! projection over the remaining iterations stays finite. The one
  ! limit named is the arithmetic's own.
  !===================================================================!

  logical function diverging(this, imbalance) result(yes)

    class(minimizer), intent(in) :: this
    real(dp)        , intent(in) :: imbalance

    real(dp) :: slope, error, remaining
    logical  :: usable

    yes = .false.
    if (imbalance <= this % began_at .or. imbalance <= 0.0_dp) return

    call this % fitted(slope, error, usable)
    if (.not. usable .or. slope <= error) return

    remaining = real(max(this % max_iterations - this % noted, 0), dp)
    yes = log(imbalance) + (slope - error) * remaining > log(huge(1.0_dp))

  end function diverging

  pure real(dp) function began(this) result(imbalance)

    class(minimizer), intent(in) :: this

    imbalance = this % began_at

  end function began

  !===================================================================!
  ! Whether the iteration has reached its limit. A limit setting
  ! that is neither by_count nor by_rate stops the program.
  !===================================================================!

  logical function exhausted(this, iteration) result(done)

    class(minimizer), intent(in) :: this
    integer         , intent(in) :: iteration

    done = iteration >= this % max_iterations

    select case (this % limit_kind)
    case (by_count)
       continue
    case (by_rate)
       done = done .or. this % flattened()
    case default
       error stop 'minimizer: an iteration limit is counted or taken from the rate'
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

  logical function halted(this, imbalance, iteration) result(stop_here)

    class(minimizer), intent(inout) :: this
    real(dp)        , intent(in)    :: imbalance
    integer         , intent(in)    :: iteration

    call this % note_imbalance(imbalance)

    stop_here = .true.
    if (this % converged(imbalance)) return
    if (imbalance /= imbalance) return
    if (imbalance > huge(1.0_dp) / 2.0_dp) return
    if (this % diverging(imbalance)) return
    if (this % exhausted(iteration)) return
    stop_here = .false.

  end function halted

  !===================================================================!
  ! Store the operation and the graph it reads. The affine part is
  ! evaluated here, once: the operation applied to the zero state is
  ! the contribution of its boundary values and sources alone.
  !===================================================================!

  subroutine attach(this, action, on, unknown_domain, num_unknowns, &
       & num_components, coupling, stored_inputs)

    class(minimizer)  , intent(inout) :: this
    class(operation), intent(in)    :: action
    class(directed_graph)          , intent(in)    :: on
    type(graph)       , intent(in)    :: unknown_domain
    integer               , intent(in)    :: num_unknowns
    integer, intent(in), optional         :: num_components
    class(directed_graph)  , intent(in), optional  :: coupling
    type(stored_field), intent(in), optional :: stored_inputs(:)

    real(dp), allocatable :: zero(:)
    integer :: n

    if (allocated(this % action)) deallocate(this % action)
    allocate(this % action, source=action)
    if (allocated(this % on)) deallocate(this % on)
    allocate(this % on, source=on)

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
    call action % domain(on, this % residual_domain, this % num_residuals)

    n = this % num_unknowns

    allocate(zero(n * this % num_components))
    zero = 0.0_dp
    call raw_apply(this, zero, this % affine)

    ! Classification is not admissibility: U and Y stay distinct
    ! identities, but THIS solver family is square - the scalar
    ! dimensions must agree. A rectangular least-squares family may
    ! support R^n -> R^m later; none exists yet.
    if (size(this % affine) /= n * this % num_components) then
       error stop 'minimization: the current solver family requires equal &
            &unknown and residual value dimensions'
    end if

    this % diagonal_valid = .false.

  end subroutine attach

  !===================================================================!
  ! The inputs a statement is evaluated on at the state x: x stored
  ! on the unknown domain, n members of the given width, then the
  ! fixed inputs where given. The residual and every tangent taken of
  ! it are built on this one tuple, so they linearize the same function.
  !===================================================================!

  function state_tuple(on, n, components, x, stored) result(inputs)

    type(graph)       , intent(in)           :: on
    integer           , intent(in)           :: n, components
    real(dp)          , intent(in)           :: x(:)
    type(stored_field), intent(in), optional :: stored(:)
    type(stored_field), allocatable          :: inputs(:)

    type(stored_field) :: state

    state = stored_field('state', on, n, num_components=components)
    call state % set_real_vector(x)

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

  subroutine raw_apply(this, x, y, inputs)

    class(minimizer), intent(in)   :: this
    real(dp), intent(in)               :: x(:)
    real(dp), allocatable, intent(out) :: y(:)
    type(stored_field), allocatable, intent(out), optional :: inputs(:)

    type(stored_field), allocatable :: tuple(:)
    class(field), allocatable :: image

    tuple = state_tuple(this % unknown_domain, this % num_unknowns, this % num_components, x, this % stored)
    call this % action % apply(this % on, this % action % bind(tuple), image)

    if (.not. image % defined_on(this % residual_domain)) then
       error stop 'minimization: the action must return a field on its stated residual domain'
    end if

    call image % real_vector(y)
    if (present(inputs)) call move_alloc(tuple, inputs)

  end subroutine raw_apply

  !===================================================================!
  ! The solver's operations, each a delegation.
  !===================================================================!

  subroutine matvec(this, x, y)

    class(minimizer), intent(in)   :: this
    real(dp), intent(in)               :: x(:)
    real(dp), allocatable, intent(out) :: y(:)

    call raw_apply(this, x, y)
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

    ! the reduction's power is two by default, so the norm is the
    ! euclidean length
    length = sqrt(sum(u * u))

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
       error stop 'minimization: a sweep needs the dependent-variable &
            &coupling - attach it with coupling='
    end if

    colouring = traversal(TRAVERSAL_COLOURING)
    call colouring % apply(this % coupling, output=image)
    call image % integer_vector(colours)

  end subroutine sweep_order

  !===================================================================!
  ! THE BLOCK DIAGONAL by coloured indicators. The coupling attached is over
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

    if (this % num_components > 1) then
       error stop 'diagonal: the coloured indicator evaluates one number per cell'
    end if

    n  = size(this % affine)
    w  = this % block_width
    nb = n / w
    if (nb * w /= n) then
       error stop 'diagonal: the unknowns come in whole blocks'
    end if

    allocate(d(w, w, nb), indicator(n))
    d = 0.0_dp

    call this % sweep_order(colours)
    if (size(colours) /= nb) then
       error stop 'diagonal: the coupling attached is over the blocks, one colour each'
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
    real(dp), allocatable :: rhs(:), x(:)
    real(dp) :: achieved

    associate (u1 => input_graph); end associate

    ! x IS a state on the unknown domain, allocated at that extent.
    allocate(x(this % num_unknowns * this % num_components))
    x = 0.0_dp

    if (present(inputs)) then
       call bound_value(inputs, this % argument(1), right_hand_side)
       if (.not. right_hand_side % defined_on(this % residual_domain)) then
          error stop 'minimization: a right-hand side is defined on the residual domain'
       end if
       call right_hand_side % real_vector(rhs)
       allocate(copy, source=this)
       call copy % solve(rhs, x, achieved)
    end if

    out = stored_field('solution', this % unknown_domain, this % num_unknowns, num_components=this % num_components)
    call out % set_real_vector(x)

    call emit(out, output)

  end subroutine solver_apply

end module operation_minimization
