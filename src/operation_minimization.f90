!=====================================================================!
! LEVEL 7 OF THE NEW TOWER . THE MINIMIZATION
!
! The first level with a goal, and its law in one line: given a
! residual map R : U -> Y, vary the values on the UNKNOWN domain U
! to drive the values on the RESIDUAL domain Y toward zero. U and Y
! are set graph identities - never assumed to be anyone's
! vertices; the graph argument survives only as the legacy
! operation host the compatibility apply() signature still wants. This
! module holds the minimizer base - ONE family for one story:
! attach a statement, drive its residual to zero. Linear solvers,
! newton, and whatever else minimizes a residual are its
! concretions; their differences are governance inside the family,
! never a second taxonomy. The solver vocabulary is defined here as
! thin delegations to engine entries, so the engine never learns a
! solver word and the solver never says apply or measure:
!
!      matvec ········· the operation applied, minus its constant
!      inner_product ·· a sum reduction with the second field as
!                       the measure
!      norm ··········· the norm reduction
!      sweep_order ···· the colouring walk
!      diagonal ······· matvec probed by colour: applied to one
!                       colour's indicator, the answer at a member
!                       IS its diagonal entry, because no two
!                       neighbours share a colour
!      constant ······· the affine part of the attached operation -
!                       what boundary values and sources contribute
!                       at zero state; the assembled right hand side
!                       is its negative
!
! A concrete solver works in plain arrays - fetched once, worked,
! written back once, as the field banner orders - and states only
! its iteration. Everything else is inherited from here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_minimization

  use iso_fortran_env       , only : dp => REAL64
  use operation_action  , only : operation
  use view_directed   , only : directed_graph
  use field_calculus  , only : field
  use graph_fractal      , only : graph
  use field_calculus        , only : functional
  use field_stored     , only : stored_field
  use operation_reduction , only : reduction, REDUCE_SUM, REDUCE_NORM
  use operation_walk      , only : walk, WALK_COLOURING

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
  ! leave one degree of freedom for the scatter and no more, so the
  ! smallest window that gives a slope an error worth comparing to.
  integer, parameter :: window = 5
  public :: minimizer

  !===================================================================!
  ! The base: an attached operation, the graph it reads, and the
  ! tolerances every iteration honours.
  !===================================================================!

  type, abstract, extends(operation) :: minimizer

     class(operation), allocatable :: action

     ! THE EXECUTION CONTEXT. The graph handed to the action when it
     ! is applied, and nothing more. It is whatever the action needs
     ! to compute with - a mesh, a compatibility host, a conduit -
     ! and it carries no authority over the solver's own structure.
     class(directed_graph)          , allocatable :: on

     ! THE DEPENDENT-VARIABLE COUPLING. Which unknowns feed which:
     ! the stencil the structural algorithms need, and the ONLY thing
     ! sweep_order is permitted to colour.
     !
     ! It is OPTIONAL AT ATTACH AND HAS NO FALLBACK. A structure-free
     ! minimizer - gmres, conjugate gradient, newton - never asks for
     ! it and need never supply one. A structured one - jacobi,
     ! gauss-seidel - is handed it by a caller that knows which object
     ! owns the dependent axis, and fails loudly if it was not.
     !
     ! `coupling := on` would be exactly the mistake this component exists
     ! to prevent: the graph an action happens to execute over is not
     ! evidence about which unknowns are coupled. Where the two really
     ! are the same graph, the CALLER says so, at its own call site.
     class(directed_graph)          , allocatable :: coupling

     ! The unknown domain U: where the answer lives, explicit at
     ! attach, identity preserved - never inferred from the host.
     type(graph) :: unknown_domain
     integer         :: num_unknowns = 0

     ! The residual domain Y: what the action answers on, asked of
     ! the action itself at attach.
     type(graph) :: residual_domain
     integer         :: num_residuals = 0

     ! A second domain, one member per NUMBER rather than per cell.
     ! The pairings live here: a measure carries one weight per
     ! entry, so a dot product over wide entries must be taken on
     ! the values themselves - the calculus says as much in its own
     ! banner. With one number per cell the two domains are the same
     ! set, and nothing changes.
     type(graph) :: numbers
     integer         :: num_numbers = 0

     ! How wide an entry is. One number per cell is the common case
     ! and the default; a state with several numbers per cell - a
     ! complex plane point, a species vector, a whole spatial field
     ! standing at one instant - says so at attach, and every word
     ! below then measures the entire vector instead of its first
     ! stripe.
     integer :: num_components = 1

     real(dp), allocatable :: affine(:)

     ! THE HELD INPUTS. Inputs held fixed during this solve - a
     ! scheme's history, an action's parameters - applied after the
     ! unknown on every evaluation. Solve-context data: the solver
     ! varies the unknown and nothing else.
     type(stored_field), allocatable :: held(:)

     integer  :: max_iterations = 1000
     real(dp) :: tolerance      = 1.0d-10

     ! WHAT THE TOLERANCE IS MEASURED AGAINST. relative divides the
     ! imbalance by the one the iteration began with, which is the
     ! question asked whenever the target is a reduction. absolute
     ! compares the imbalance itself, which is the question asked only
     ! when the number has a meaning of its own - a functional driven
     ! under a stated value, and nothing else.
     integer :: criterion = relative

     ! WHERE THE ITERATION STOPS WHEN IT HAS NOT CONVERGED. by_count
     ! stops at max_iterations. by_rate also stops where the imbalance
     ! has flattened, which is where the slope of its logarithm over
     ! the last few iterations is no longer distinguishable from zero
     ! against its own scatter.
     integer :: budget = by_count

     ! The imbalance the iteration began at, and the last few it has
     ! seen. Written by note_imbalance and by nothing else.
     real(dp), private :: began_at = 0.0_dp
     real(dp), private :: recent(window) = 0.0_dp
     integer , private :: noted = 0

     ! Whether a window has yet shown a slope significantly below
     ! zero. An iteration that has never descended has not flattened
     ! either, whatever a window of its early wandering looks like.
     logical , private :: descended = .false.

   contains

     procedure :: begin_imbalance
     procedure :: note_imbalance
     procedure :: converged
     procedure :: flattened
     procedure, private :: fitted
     procedure :: exhausted

     procedure :: attach
     procedure :: evaluation_inputs
     procedure :: matvec
     procedure :: inner_product
     procedure :: norm
     procedure :: sweep_order
     procedure :: diagonal
     procedure :: constant

     ! The operation face: a solver IS an operation - the one that
     ! answers the attached statement. apply solves from zero, so a
     ! solver composes wherever operations go; a preconditioner is
     ! exactly this face of an inner solver.
     procedure :: domain => solver_domain
     procedure :: apply  => solver_apply

     procedure(solve_interface), deferred :: solve

  end type minimizer

  abstract interface

     !----------------------------------------------------------------!
     ! Drive || rhs - matvec(x) || under the tolerance, within the
     ! iteration budget. The achieved norm reports the truth either
     ! way.
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
  ! tolerance is measured, and the last few it has seen, from which a
  ! rate is fitted. Written here and nowhere else.
  !===================================================================!

  subroutine begin_imbalance(this, imbalance)

    class(minimizer), intent(inout) :: this
    real(dp)        , intent(in)    :: imbalance

    this % began_at  = imbalance
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
  ! Whether the imbalance meets what was asked of it. Relative divides
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
  ! error; scatter about a floor has a slope inside it. Nothing is
  ! chosen here except the width of the window, and no scale enters,
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
  ! Whether the iteration has run out of what it was given. A budget
  ! that is neither counted nor taken from the rate stops the program.
  !===================================================================!

  logical function exhausted(this, iteration) result(done)

    class(minimizer), intent(in) :: this
    integer         , intent(in) :: iteration

    done = iteration >= this % max_iterations

    select case (this % budget)
    case (by_count)
       continue
    case (by_rate)
       done = done .or. this % flattened()
    case default
       error stop 'minimizer: a budget is counted or taken from the rate'
    end select

  end function exhausted

  !===================================================================!
  ! Take the operation and the graph it reads. The affine part is
  ! measured here, once: the operation applied to nothing is what
  ! its boundary values and sources say by themselves.
  !===================================================================!

  subroutine attach(this, action, on, unknown_domain, num_unknowns, &
       & num_components, coupling, held_inputs)

    class(minimizer)  , intent(inout) :: this
    class(operation), intent(in)    :: action
    class(directed_graph)          , intent(in)    :: on
    type(graph)       , intent(in)    :: unknown_domain
    integer               , intent(in)    :: num_unknowns
    integer, intent(in), optional         :: num_components
    class(directed_graph)  , intent(in), optional  :: coupling
    type(stored_field), intent(in), optional :: held_inputs(:)

    real(dp), allocatable :: zero(:)
    integer :: n

    if (allocated(this % action)) deallocate(this % action)
    allocate(this % action, source=action)
    if (allocated(this % on)) deallocate(this % on)
    allocate(this % on, source=on)

    ! the inputs held fixed during this solve follow the unknown in
    ! every evaluation, the affine part's included
    if (allocated(this % held)) deallocate(this % held)
    if (present(held_inputs)) allocate(this % held, source=held_inputs)

    ! the solver's own operation face reads one input, the right-hand side
    call this % declare_arguments(1)

    ! The dependent-variable coupling arrives EXPLICIT or not at all.
    ! No fallback to the execution context: a solver that needs
    ! structure and was given none says so when it reaches for it,
    ! rather than colouring whatever graph happened to be nearby.
    if (allocated(this % coupling)) deallocate(this % coupling)
    if (present(coupling)) allocate(this % coupling, source=coupling)

    this % num_components = 1
    if (present(num_components)) this % num_components = max(num_components, 1)

    ! The unknown domain arrives EXPLICIT and identity-preserving.
    ! No hidden fallback to the host's vertices: a caller that
    ! means vertices says so at its own call site.
    this % unknown_domain   = unknown_domain
    this % num_unknowns = num_unknowns

    ! The residual domain is the action's own answer.
    call action % domain(on, this % residual_domain, this % num_residuals)

    !----------------------------------------------------------------!
    ! attach is re-enterable - Newton calls it once per iteration - and
    ! a graph signs ONCE. The old counted_set constructor minted a
    ! fresh number domain on every attach, so a fresh one is minted
    ! here too, by resetting the component to an unsigned graph before
    ! declaring it. Signing the same variable twice is refused, and
    ! rightly; this says which of the two meanings was intended.
    !----------------------------------------------------------------!

    n = this % num_unknowns

    block
      type(graph) :: unsigned
      this % numbers = unsigned
    end block
    call this % numbers % declare()

    this % num_numbers = n * this % num_components

    allocate(zero(n * this % num_components))
    zero = 0.0_dp
    call raw_apply(this, zero, this % affine)

    ! Classification is not admissibility: U and Y stay distinct
    ! identities, but THIS solver family is square - the scalar
    ! dimensions must agree. A rectangular least-squares family may
    ! earn R^n -> R^m later; it has not yet.
    if (size(this % affine) /= n * this % num_components) then
       error stop 'minimization: the current solver family requires equal &
            &unknown and residual value dimensions'
    end if

  end subroutine attach

  !===================================================================!
  ! The inputs the statement is evaluated on at the unknown x: the
  ! state on the unknown domain, then the held inputs. The residual
  ! and every tangent taken of it are built on this one tuple, so
  ! they linearize the same function.
  !===================================================================!

  subroutine evaluation_inputs(this, x, inputs)

    class(minimizer), intent(in) :: this
    real(dp), intent(in)         :: x(:)
    type(stored_field), allocatable, intent(out) :: inputs(:)

    type(stored_field) :: state

    state = stored_field('state', this % unknown_domain, this % num_unknowns, num_components=this % num_components)
    call state % set_real_vector(x)

    if (allocated(this % held)) then
       inputs = [state, this % held]
    else
       inputs = [state]
    end if

  end subroutine evaluation_inputs

  !===================================================================!
  ! The operation applied as it stands, affine part and all.
  !===================================================================!

  subroutine raw_apply(this, x, y)

    class(minimizer), intent(in)   :: this
    real(dp), intent(in)               :: x(:)
    real(dp), allocatable, intent(out) :: y(:)

    type(stored_field), allocatable :: inputs(:)
    class(field), allocatable :: answer

    call this % evaluation_inputs(x, inputs)
    call this % action % apply(this % on, inputs, answer)

    if (.not. answer % defined_on(this % residual_domain)) then
       error stop 'minimization: the action must answer on its stated residual domain'
    end if

    call answer % real_vector(y)

  end subroutine raw_apply

  !===================================================================!
  ! The solver's words, each a delegation.
  !===================================================================!

  subroutine matvec(this, x, y)

    class(minimizer), intent(in)   :: this
    real(dp), intent(in)               :: x(:)
    real(dp), allocatable, intent(out) :: y(:)

    call raw_apply(this, x, y)
    y = y - this % affine

  end subroutine matvec

  real(dp) function inner_product(this, u, v) result(prod)

    class(minimizer), intent(in) :: this
    real(dp), intent(in) :: u(:), v(:)

    type(reduction) :: total
    type(stored_field) :: uf, vf
    class(functional), allocatable :: answer
    real(dp), allocatable :: got(:)

    uf = stored_field('u', this % numbers, this % num_numbers)
    call uf % set_real_vector(u)
    vf = stored_field('v', this % numbers, this % num_numbers)
    call vf % set_real_vector(v)

    total = reduction(REDUCE_SUM)
    call total % reduce(uf, answer, measure=vf)

    call answer % real_vector(got)
    prod = got(1)

  end function inner_product

  real(dp) function norm(this, u) result(length)

    class(minimizer), intent(in) :: this
    real(dp), intent(in) :: u(:)

    type(reduction) :: measure_of
    type(stored_field) :: uf
    class(functional), allocatable :: answer
    real(dp), allocatable :: got(:)

    uf = stored_field('u', this % numbers, this % num_numbers)
    call uf % set_real_vector(u)

    measure_of = reduction(REDUCE_NORM)
    call measure_of % reduce(uf, answer)

    call answer % real_vector(got)
    length = got(1)

  end function norm

  subroutine sweep_order(this, colours)

    class(minimizer), intent(in)  :: this
    integer, allocatable, intent(out) :: colours(:)

    type(walk) :: colouring
    class(field), allocatable :: answer

    ! THE COLOURING IS OF THE UNKNOWNS' COUPLING, never of the
    ! execution context. Two unknowns may share a colour only when
    ! nothing couples them, and the graph an action runs over knows
    ! nothing about that.
    if (.not. allocated(this % coupling)) then
       error stop 'minimization: a sweep needs the dependent-variable &
            &coupling - attach it with coupling='
    end if

    colouring = walk(WALK_COLOURING)
    call colouring % apply(this % coupling, output=answer)
    call answer % integer_vector(colours)

  end subroutine sweep_order

  !===================================================================!
  ! The diagonal, probed by colour. Applied to the indicator of one
  ! colour class, the answer at a member is that member's diagonal
  ! entry, because none of its neighbours is in the class. One
  ! matvec per colour, a handful in all.
  !===================================================================!

  subroutine diagonal(this, d)

    class(minimizer), intent(in)   :: this
    real(dp), allocatable, intent(out) :: d(:)

    integer , allocatable :: colours(:)
    real(dp), allocatable :: indicator(:), y(:)
    integer :: nv, col, v

    if (this % num_components > 1) then
       ! The probe reads one answer per cell, and a wide entry has
       ! several. A block probe is the honest generalization and no
       ! caller has asked for one yet.
       error stop 'diagonal: the coloured probe answers one number per cell'
    end if

    nv = size(this % affine)
    allocate(d(nv), indicator(nv))
    d = 0.0_dp

    call this % sweep_order(colours)

    do col = 1, maxval(colours)

       indicator = 0.0_dp
       do v = 1, nv
          if (colours(v) == col) indicator(v) = 1.0_dp
       end do

       call this % matvec(indicator, y)

       do v = 1, nv
          if (colours(v) == col) d(v) = y(v)
       end do

    end do

  end subroutine diagonal

  !===================================================================!
  ! The affine part, for the caller assembling an equation: the
  ! statement action(q) = 0 reads matvec(q) = -constant.
  !===================================================================!

  subroutine constant(this, g)

    class(minimizer), intent(in)   :: this
    real(dp), allocatable, intent(out) :: g(:)

    g = this % affine

  end subroutine constant

  !===================================================================!
  ! The operation face.
  !===================================================================!

  subroutine solver_domain(this, input_graph, domain, num_entries)

    class(minimizer), intent(in)       :: this
    class(directed_graph), intent(in)               :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries

    associate (u1 => input_graph); end associate

    ! The solver's answer is a solution on U.
    domain   = this % unknown_domain
    num_entries = this % num_unknowns

  end subroutine solver_domain

  subroutine solver_apply(this, input_graph, input_data, output)

    class(minimizer), intent(in)                   :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    class(minimizer), allocatable :: worker
    type(stored_field) :: out
    real(dp), allocatable :: rhs(:), x(:)
    real(dp) :: achieved

    associate (u1 => input_graph); end associate

    ! x IS a state on the unknown domain; say so.
    allocate(x(this % num_unknowns * this % num_components))
    x = 0.0_dp

    if (present(input_data)) then
       if (.not. input_data(1) % defined_on(this % residual_domain)) then
          error stop 'minimization: a right-hand side lives on the residual domain'
       end if
       call input_data(1) % real_vector(rhs)
       allocate(worker, source=this)
       call worker % solve(rhs, x, achieved)
    end if

    out = stored_field('solution', this % unknown_domain, this % num_unknowns, num_components=this % num_components)
    call out % set_real_vector(x)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine solver_apply

end module operation_minimization
