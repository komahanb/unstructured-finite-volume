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
  use graph_ordinary_view, only : graph
  use graph_field_calculus, only : graph_field
  use fractal_graph      , only : set_graph => graph
  use graph_calculus     , only : linearization_operator
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
     module procedure create
  end interface difference_linearization

contains

  !===================================================================!
  ! Built about a statement; the state and base may arrive now or
  ! through freeze.
  !===================================================================!

  function create(of, at, base) result(this)

    class(graph_operation), intent(in) :: of
    real(dp), intent(in), optional     :: at(:)
    real(dp), intent(in), optional     :: base(:)
    type(difference_linearization)     :: this

    allocate(this % of, source=of)
    if (present(at))   allocate(this % at  , source=at)
    if (present(base)) allocate(this % base, source=base)

  end function create

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
    class(graph), intent(in)                    :: input_graph
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
  ! is all of them on the ordinary-graph road - the domain asked for
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
    class(graph), intent(in)                       :: input_graph
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

end module class_graph_linearization
