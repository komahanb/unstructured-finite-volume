!=====================================================================!
! The difference linearization: a derivative bought with two
! residuals.
!
! LEVEL 1 OF THE STRATIFICATION. The first concretion of the
! linearization operator: the tangent of a statement S at the
! standing state q, computed as a difference,
!
!      J v  ~  ( S(q + eps v) - S(q) ) / eps
!
! wrapped as an operation, so whoever governs it - newton, one rank
! up - sees an ordinary linear question and asks no questions. The
! difference buys generality and pays in precision: the residual
! floor is the machine epsilon over the step times the residual
! scale, about eight digits. A statement that can linearize itself
! exactly joins this family as a second concretion when it arrives,
! and the governor never notices the change.
!
! freeze moves the standing state between the governor's steps; the
! base residual rides along when the caller already holds it, and is
! measured here once when it does not.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_linearization

  use iso_fortran_env    , only : dp => REAL64
  use graph_operation_view, only : graph_operation
  use graph_directed_view, only : directed_graph
  use graph_field_calculus, only : graph_field
  use fractal_graph      , only : set_graph => graph
  use graph_discretization     , only : linearization_operator, &
       & differentiable_operation
  use class_graph_field  , only : field

  implicit none

  private
  public :: difference_linearization

  type, extends(linearization_operator) :: difference_linearization

     class(graph_operation), allocatable :: of

     real(dp), allocatable :: at(:)
     real(dp), allocatable :: base(:)
     real(dp) :: step = 1.0d-7

   contains

     procedure :: name   => derivative_name
     procedure :: domain => derivative_domain
     procedure :: apply  => derivative_apply
     procedure :: freeze => derivative_freeze

  end type difference_linearization

  interface difference_linearization
     module procedure create_difference
  end interface difference_linearization

!=====================================================================!
! Exact linearization: the second concrete member of the
! linearization family. The tangent of a differentiable statement
! S at the frozen state q is one partial action,
!
!      J v = D S(q) [v],
!
! exact to the statement's declared max_degree - no finite
! difference. freeze moves the frozen state; a base residual
! passed to freeze is accepted and ignored, since an exact tangent
! does not use the residual value. Like the difference member,
! this linearizes S : U -> U in its first input slot only.
!
!=====================================================================!

  public :: exact_linearization
  public :: tangent_of

  type, extends(linearization_operator) :: exact_linearization

     class(differentiable_operation), allocatable :: of

     real(dp), allocatable :: at(:)

   contains

     procedure :: name   => exact_name
     procedure :: domain => exact_domain
     procedure :: apply  => exact_apply
     procedure :: freeze => exact_freeze

  end type exact_linearization

  interface exact_linearization
     module procedure create_exact
  end interface exact_linearization

contains

  !===================================================================!
  ! Built about a statement; the state and base may arrive now or
  ! through freeze.
  !===================================================================!

  function create_difference(of, at, base) result(this)

    class(graph_operation), intent(in) :: of
    real(dp), intent(in), optional     :: at(:)
    real(dp), intent(in), optional     :: base(:)
    type(difference_linearization)     :: this

    allocate(this % of, source=of)
    if (present(at))   allocate(this % at  , source=at)
    if (present(base)) allocate(this % base, source=base)

  end function create_difference

  !===================================================================!
  ! Move the standing state. The base residual is stored when handed
  ! over, and forgotten otherwise so apply measures it fresh.
  !===================================================================!

  subroutine derivative_freeze(this, at, base)

    class(difference_linearization), intent(inout) :: this
    real(dp), intent(in)           :: at(:)
    real(dp), intent(in), optional :: base(:)

    this % at = at

    if (present(base)) then
       this % base = base
    else
       if (allocated(this % base)) deallocate(this % base)
    end if

  end subroutine derivative_freeze

  pure function derivative_name(this) result(name)

    class(difference_linearization), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'derivative of ' // this % of % name()

  end function derivative_name

  subroutine derivative_domain(this, input_graph, domain, nentries)

    class(difference_linearization), intent(in) :: this
    class(directed_graph), intent(in)                    :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries

    call this % of % domain(input_graph, domain, nentries)

  end subroutine derivative_domain

  !===================================================================!
  ! J v as two residuals and a difference. The direction arrives as
  ! the input field; the answer leaves on the same cells.
  !===================================================================!

  !===================================================================!
  ! THE DIFFERENCE IS TAKEN WHERE THE OPERATION LIVES.
  !
  ! derivative_domain has always delegated to the underlying
  ! operation; the states this routine builds now do the same. A
  ! finite difference of an operation is a statement about that
  ! operation's own domain, and the graph it is reached through was
  ! never the seat of it.
  !
  ! For every operation that reads its domain off the graph - which
  ! is all of them on the directed-graph road - the domain asked for
  ! here is exactly the vertex set that was being used before, so
  ! nothing a graph-based caller sees changes.
  !
  ! THIS CITIZEN REMAINS SAME-DOMAIN. It differences L : U -> U and
  ! nothing wider. No rectangular U -> Y, no transpose, no reverse
  ! action: those would be a different mathematical object, and this
  ! correction deliberately does not become one.
  !===================================================================!

  subroutine derivative_apply(this, input_graph, input_data, output)

    class(difference_linearization), intent(in)    :: this
    class(directed_graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(field)   :: state
    class(graph_field), allocatable :: pushed
    type(set_graph) :: on, given
    integer         :: n_on, n_given
    real(dp), allocatable :: v(:), y(:), base(:)
    integer :: n, ncomp

    call this % of % domain(input_graph, on, n_on)

    n = n_on
    if (n <= 0) then
       error stop 'linearization: the operation''s domain is empty'
    end if
    if (mod(size(this % at), n) /= 0) then
       error stop 'linearization: the frozen state must carry a whole number &
            &per member of the operation''s domain'
    end if

    ! The width the frozen state carries. A statement of several
    ! numbers per member is differenced whole, like any other.
    ncomp = max(size(this % at) / n, 1)

    if (present(input_data)) then
       given   = input_data(1) % domain()
       n_given = input_data(1) % num_entries()
       if (.not. given % same_as(on)) then
          error stop 'linearization: the direction must live on the operation''s domain'
       end if
       call input_data(1) % get_real_vector(v)
       if (size(v) /= size(this % at)) then
          error stop 'linearization: the direction must match the frozen state''s width'
       end if
    else
       allocate(v(n * ncomp))
       v = 0.0_dp
    end if

    ! The base residual: taken from the freeze when it rode along,
    ! measured here once when it did not.
    if (allocated(this % base)) then
       base = this % base
    else
       state = field('state', on, n_on, ncomp=ncomp)
       call state % set_real_vector(this % at)
       call this % of % apply(input_graph, [state], pushed)
       call answered_on(pushed, on)
       call pushed % get_real_vector(base)
    end if

    state = field('state', on, n_on, ncomp=ncomp)
    call state % set_real_vector(this % at + this % step * v)

    call this % of % apply(input_graph, [state], pushed)
    call answered_on(pushed, on)
    call pushed % get_real_vector(y)

    y = (y - base) / this % step

    state = field('J v', on, n_on, ncomp=ncomp)
    call state % set_real_vector(y)
    if (allocated(output)) deallocate(output)
    allocate(output, source=state)

  end subroutine derivative_apply

  !===================================================================!
  ! A same-domain difference subtracts two answers, so both must
  ! have come from the same place. Equal length is not that claim.
  !===================================================================!

  subroutine answered_on(answer, expected)

    class(graph_field), intent(in) :: answer
    type(set_graph)   , intent(in) :: expected

    type(set_graph) :: got

    got = answer % domain()
    if (.not. got % same_as(expected)) then
       error stop 'linearization: the operation must answer on its stated domain'
    end if

  end subroutine answered_on


  !===================================================================!
  ! Select the linearization by the statement's type: exact for a
  ! differentiable_operation, difference for anything else.
  ! Callers (newton, the marcher) use this instead of dispatching
  ! themselves, so the dispatch rule exists in one place.
  !===================================================================!

  function tangent_of(action) result(tangent)

    class(graph_operation), intent(in)         :: action
    class(linearization_operator), allocatable :: tangent

    select type (action)
    class is (differentiable_operation)
       allocate(tangent, source=exact_linearization(action))
    class default
       allocate(tangent, source=difference_linearization(action))
    end select

  end function tangent_of

  !===================================================================!
  ! Construct from a differentiable statement; the state may
  ! arrive now or later through freeze.
  !===================================================================!

  function create_exact(of, at) result(this)

    class(differentiable_operation), intent(in) :: of
    real(dp), intent(in), optional              :: at(:)
    type(exact_linearization)                   :: this

    allocate(this % of, source=of)
    if (present(at)) allocate(this % at, source=at)

  end function create_exact

  !===================================================================!
  ! Move the frozen state. base is part of the family's freeze
  ! signature and is ignored here.
  !===================================================================!

  subroutine exact_freeze(this, at, base)

    class(exact_linearization), intent(inout) :: this
    real(dp), intent(in)           :: at(:)
    real(dp), intent(in), optional :: base(:)

    if (present(base)) then
       associate(unread => base)
       end associate
    end if

    this % at = at

  end subroutine exact_freeze

  pure function exact_name(this) result(name)

    class(exact_linearization), intent(in) :: this
    character(len=:), allocatable :: name

    name = 'exact derivative of ' // this % of % name()

  end function exact_name

  subroutine exact_domain(this, input_graph, domain, nentries)

    class(exact_linearization), intent(in) :: this
    class(directed_graph), intent(in)      :: input_graph
    type(set_graph), intent(out) :: domain
    integer        , intent(out) :: nentries

    call this % of % domain(input_graph, domain, nentries)

  end subroutine exact_domain

  !===================================================================!
  ! Apply J v as one partial action in input slot 1 at the frozen
  ! state. Checks, each stopping the program: the domain must be
  ! nonempty; a state must have been frozen; the frozen state must
  ! hold a whole number of components per domain member; the
  ! direction must live on the statement's domain and match the
  ! frozen state's size. Without input data the result is zero,
  ! returned without calling the statement.
  !===================================================================!

  subroutine exact_apply(this, input_graph, input_data, output)

    class(exact_linearization), intent(in)          :: this
    class(directed_graph), intent(in)               :: input_graph
    class(graph_field), intent(in), optional        :: input_data(:)
    class(graph_field), allocatable, intent(inout)  :: output

    type(field)     :: state, direction, zero
    type(set_graph) :: on, given
    real(dp), allocatable :: v(:)
    integer         :: n_on
    integer         :: n, ncomp

    call this % of % domain(input_graph, on, n_on)

    n = n_on
    if (n <= 0) then
       error stop 'linearization: the operation''s domain is empty'
    end if
    if (.not. allocated(this % at)) then
       error stop 'linearization: an exact tangent is taken at a frozen state'
    end if
    if (mod(size(this % at), n) /= 0) then
       error stop 'linearization: the frozen state must carry a whole number &
            &per member of the operation''s domain'
    end if

    ncomp = max(size(this % at) / n, 1)

    state = field('state', on, n_on, ncomp=ncomp)
    call state % set_real_vector(this % at)

    if (present(input_data)) then

       given = input_data(1) % domain()
       if (.not. given % same_as(on)) then
          error stop 'linearization: the direction must live on the operation''s domain'
       end if
       call input_data(1) % get_real_vector(v)
       if (size(v) /= size(this % at)) then
          error stop 'linearization: the direction must match the frozen state''s width'
       end if

       ! copy into a concrete field: a polymorphic entity cannot
       ! appear in an array constructor
       direction = field('direction', on, n_on, ncomp=ncomp)
       call direction % set_real_vector(v)

       call this % of % partial_action(input_graph, [state], [1], &
            & [direction], output)

    else

       ! no direction means the zero direction, and the tangent of
       ! zero is zero
       zero = field('J v', on, n_on, ncomp=ncomp)
       call zero % set_real_vector(spread(0.0_dp, 1, size(this % at)))
       if (allocated(output)) deallocate(output)
       allocate(output, source=zero)

    end if

  end subroutine exact_apply


end module class_graph_linearization
