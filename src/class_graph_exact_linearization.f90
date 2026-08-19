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
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_exact_linearization

  use iso_fortran_env     , only : dp => REAL64
  use graph_directed_view , only : directed_graph
  use graph_field_calculus, only : graph_field
  use fractal_graph       , only : set_graph => graph
  use graph_calculus      , only : linearization_operator, &
       & differentiable_operation
  use graph_operation_view, only : graph_operation
  use class_graph_linearization, only : difference_linearization
  use class_graph_field   , only : field

  implicit none

  private
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
     module procedure create
  end interface exact_linearization

contains

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

  function create(of, at) result(this)

    class(differentiable_operation), intent(in) :: of
    real(dp), intent(in), optional              :: at(:)
    type(exact_linearization)                   :: this

    allocate(this % of, source=of)
    if (present(at)) allocate(this % at, source=at)

  end function create

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

end module class_graph_exact_linearization
