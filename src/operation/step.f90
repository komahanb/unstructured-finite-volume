!=====================================================================!
! The step operator: time's discretization stencil.
!
! LEVEL 1 OF THE STRATIFICATION. One step of a time recurrence, as
! an operation - the temporal concretion of the discretization
! operator, sibling of the spatial stencil:
!
!      a0 q  +  a1 qold  +  a2 qolder  +  h S(q)
!
! The scheme is a MOTIF stamped along the chain: backward euler
! reaches one instant back, bdf-2 reaches two, and the dependency
! pattern this type answers by contract is exactly that motif - a
! little chain of reach + 1 instants. The triangularity a causal
! solver exploits is made here, by the maker, not assumed there, by
! the solver.
!
! The shelf of schemes is this module's constructor roster - the
! module is the library and the compiler is its index:
!
!      backward_euler(action, h)          reach 1, [1, -1]
!      bdf(k, action, h)                  reach k, the k-step table
!
! One name per family, the order riding as an argument, never as a
! type. The statement it holds returns MINUS the velocity, matching
! the house convention that a balance measures what a cell has left
! over.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_step

  use iso_fortran_env    , only : dp => REAL64
  use structure_graph, only : graph
  use data_graph_field, only : graph_field
  use operation_graph_operation, only : graph_operation
  use operation_graph_operation, only : graph_discretization
  use structure_graph, only : GRAPH_SIDE_VERTEX
  use structure_support, only : support
  use data_field  , only : field
  use structure_stored_graph        , only : stored_graph

  implicit none

  private
  public :: step, chain
  public :: backward_euler, forward_euler, bdf
  public :: CHAIN_FORWARD, CHAIN_BACKWARD, CHAIN_BDF2

  !-------------------------------------------------------------------!
  ! The three rules a chain may follow, absorbed as answers.
  !-------------------------------------------------------------------!

  integer, parameter :: CHAIN_FORWARD  = 1
  integer, parameter :: CHAIN_BACKWARD = 2
  integer, parameter :: CHAIN_BDF2     = 3

  type, extends(graph_discretization) :: step

     class(graph_operation), allocatable :: action

     real(dp), allocatable :: qold(:)
     real(dp), allocatable :: qolder(:)

     real(dp) :: a0 = 1.0_dp
     real(dp) :: a1 = -1.0_dp
     real(dp) :: a2 = 0.0_dp
     real(dp) :: hs = 0.0_dp

     integer :: reach = 1

     !----------------------------------------------------------------!
     ! Where the velocity is read. An explicit rule reads it at the
     ! state it came FROM, which is what makes its row linear in the
     ! state it is solving for - and what makes the whole chain
     ! walkable in one pass.
     !----------------------------------------------------------------!

     logical :: explicit = .false.

   contains

     procedure :: name         => step_name
     procedure :: domain       => step_domain
     procedure :: apply        => step_apply
     procedure :: dependencies => step_dependencies

  end type step

  !===================================================================!
  ! THE CHAIN OPERATOR. The whole recurrence, as ONE statement on the
  ! time graph: a trajectory arrives, one residual per instant
  ! leaves.
  !
  !      instant 1 ····· q_1 - initial          the held instant
  !      instant n ····· a0 q_n + a1 q_(n-1) + a2 q_(n-2) + h S
  !
  ! Every row reads its own instant and instants already behind it,
  ! and NOTHING ahead. That is block lower triangularity, and it is
  ! made here rather than assumed by whoever solves it - which is
  ! why one causal sweep answers the whole trajectory exactly, with
  ! no iteration at all.
  !
  ! The spatial state rides as the components of an instant: the
  ! chain's vertices are moments, and each moment carries a whole
  ! field.
  !===================================================================!

  type, extends(graph_discretization) :: chain

     class(graph_operation), allocatable :: action
     class(graph)          , allocatable :: space

     real(dp), allocatable :: initial(:)

     real(dp) :: hs   = 1.0_dp
     integer  :: rule = 1

   contains

     procedure :: name         => chain_name
     procedure :: domain       => chain_domain
     procedure :: apply        => chain_apply
     procedure :: dependencies => chain_dependencies

     procedure :: row

  end type chain

  interface chain
     module procedure create_chain
  end interface chain

contains

  !===================================================================!
  ! The shelf. Backward euler: one instant back.
  !===================================================================!

  function backward_euler(action, h) result(this)

    class(graph_operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(step)                :: this

    allocate(this % action, source=action)
    this % a0    = 1.0_dp
    this % a1    = -1.0_dp
    this % a2    = 0.0_dp
    this % hs    = h
    this % reach = 1

  end function backward_euler

  !===================================================================!
  ! Forward euler: the velocity read where the step began.
  !===================================================================!

  function forward_euler(action, h) result(this)

    class(graph_operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(step)                :: this

    allocate(this % action, source=action)
    this % a0       = 1.0_dp
    this % a1       = -1.0_dp
    this % a2       = 0.0_dp
    this % hs       = h
    this % reach    = 1
    this % explicit = .true.

  end function forward_euler

  !===================================================================!
  ! The bdf family: k instants back, one name, k absorbed. Order one
  ! is backward euler by another road; order two is the second-order
  ! table; higher orders join here when their tables are ruled.
  !===================================================================!

  function bdf(k, action, h) result(this)

    integer, intent(in)                :: k
    class(graph_operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(step)                :: this

    allocate(this % action, source=action)
    this % hs    = h
    this % reach = k

    select case (k)
    case (1)
       this % a0 = 1.0_dp
       this % a1 = -1.0_dp
       this % a2 = 0.0_dp
    case (2)
       this % a0 = 1.5_dp
       this % a1 = -2.0_dp
       this % a2 = 0.5_dp
    case default
       error stop 'bdf: only orders one and two carry tables so far'
    end select

  end function bdf

  !===================================================================!
  ! The contract's answer: the motif - a chain of reach + 1
  ! instants, each edge one look further back.
  !===================================================================!

  subroutine step_dependencies(this, pattern)

    class(step), intent(in)       :: this
    class(graph), allocatable, intent(out) :: pattern

    integer :: n

    allocate(pattern, source=stored_graph(this % reach + 1, &
         & tails=[(n, n = 1, this % reach)], &
         & heads=[(n + 1, n = 1, this % reach)]))

  end subroutine step_dependencies

  !===================================================================!
  ! The step's three answers as an operation.
  !===================================================================!

  pure function step_name(this) result(name)

    class(step), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'step'

  end function step_name

  subroutine step_domain(this, input_graph, domain)

    class(step), intent(in)       :: this
    class(graph), intent(in)               :: input_graph
    class(graph), allocatable, intent(out) :: domain

    associate (u1 => this); end associate

    call input_graph % all_vertices(domain)

  end subroutine step_domain

  subroutine step_apply(this, input_graph, input_data, output)

    class(step), intent(in)               :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(support) :: cells
    type(field)   :: out
    class(graph_field), allocatable :: velocity
    real(dp), allocatable :: q(:), s(:), y(:)
    integer :: nv, ncomp, v

    nv = input_graph % num_vertices()

    if (present(input_data)) then

       call input_data(1) % get_real_vector(q)

       if (this % explicit) then
          ! The velocity where the step began, not where it lands.
          call read_at(this % action, input_graph, this % qold, s)
       else
          call this % action % apply(input_graph, input_data, velocity)
          call velocity % get_real_vector(s)
       end if

       y = this % a0 * q + this % a1 * this % qold + this % hs * s
       if (allocated(this % qolder) .and. abs(this % a2) > 0.0_dp) then
          y = y + this % a2 * this % qolder
       end if

    else
       allocate(y(nv))
       y = 0.0_dp
    end if

    ncomp = size(y) / max(nv, 1)
    cells = support(GRAPH_SIDE_VERTEX, [(v, v = 1, nv)])
    out = field('step residual', cells, ncomp=ncomp)
    call out % set_real_vector(y)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine step_apply

  !===================================================================!
  ! One statement read at a given state.
  !===================================================================!

  subroutine read_at(action, on, q, s)

    class(graph_operation), intent(in) :: action
    class(graph), intent(in)           :: on
    real(dp), intent(in)               :: q(:)
    real(dp), allocatable, intent(out) :: s(:)

    type(support) :: cells
    type(field)   :: state
    class(graph_field), allocatable :: answer
    integer :: nv, ncomp, v

    nv    = on % num_vertices()
    ncomp = size(q) / max(nv, 1)

    cells = support(GRAPH_SIDE_VERTEX, [(v, v = 1, nv)])
    state = field('state', cells, ncomp=ncomp)
    call state % set_real_vector(q)

    call action % apply(on, [state], answer)
    call answer % get_real_vector(s)

  end subroutine read_at

  !===================================================================!
  ! Build the recurrence: a statement, the space it reads, the step,
  ! the rule, and the instant that is held.
  !===================================================================!

  function create_chain(action, space, h, rule, initial) result(this)

    class(graph_operation), intent(in) :: action
    class(graph), intent(in)           :: space
    real(dp), intent(in)               :: h
    integer , intent(in)               :: rule
    real(dp), intent(in)               :: initial(:)
    type(chain)               :: this

    allocate(this % action , source=action)
    allocate(this % space  , source=space)
    allocate(this % initial, source=initial)
    this % hs   = h
    this % rule = rule

  end function create_chain

  !===================================================================!
  ! One row of the recurrence, as a step operator standing at the
  ! given instant. The rule picks the table; the first step of a
  ! two-step rule is taken by the one-step rule, as it must be.
  !===================================================================!

  type(step) function row(this, n, qold, qolder) result(block)

    class(chain), intent(in) :: this
    integer , intent(in)              :: n
    real(dp), intent(in)              :: qold(:)
    real(dp), intent(in), optional    :: qolder(:)

    select case (this % rule)
    case (CHAIN_FORWARD)
       block = forward_euler(this % action, this % hs)
    case (CHAIN_BDF2)
       if (n > 2 .and. present(qolder)) then
          block = bdf(2, this % action, this % hs)
          block % qolder = qolder
       else
          block = bdf(1, this % action, this % hs)
       end if
    case default
       block = backward_euler(this % action, this % hs)
    end select

    block % qold = qold

  end function row

  pure function chain_name(this) result(name)

    class(chain), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'chain'

  end function chain_name

  subroutine chain_domain(this, input_graph, domain)

    class(chain), intent(in)      :: this
    class(graph), intent(in)               :: input_graph
    class(graph), allocatable, intent(out) :: domain

    associate (u1 => this); end associate

    call input_graph % all_vertices(domain)

  end subroutine chain_domain

  !===================================================================!
  ! The dependency pattern IS the chain the caller marches on: each
  ! instant leans on the one before it.
  !===================================================================!

  subroutine chain_dependencies(this, pattern)

    class(chain), intent(in)      :: this
    class(graph), allocatable, intent(out) :: pattern

    integer :: reach, n

    reach = 1
    if (this % rule == CHAIN_BDF2) reach = 2

    allocate(pattern, source=stored_graph(reach + 1, &
         & tails=[(n, n = 1, reach)], heads=[(n + 1, n = 1, reach)]))

  end subroutine chain_dependencies

  !===================================================================!
  ! A trajectory in, one residual per instant out.
  !===================================================================!

  subroutine chain_apply(this, input_graph, input_data, output)

    class(chain), intent(in)              :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(support) :: instants
    type(field)   :: out
    type(step) :: block
    class(graph_field), allocatable :: answer
    type(field) :: at
    type(support) :: cells
    real(dp), allocatable :: q(:), y(:), rn(:)
    integer :: ninstants, width, n, v, lo, hi

    ninstants = input_graph % num_vertices()
    width     = size(this % initial)

    allocate(y(ninstants * width))
    y = 0.0_dp

    if (present(input_data)) then

       call input_data(1) % get_real_vector(q)

       ! The held instant: its row says only that it stands where it
       ! was put.
       y(1 : width) = q(1 : width) - this % initial

       cells = support(GRAPH_SIDE_VERTEX, &
            & [(v, v = 1, this % space % num_vertices())])

       do n = 2, ninstants

          lo = (n - 1) * width + 1
          hi = n * width

          if (n > 2) then
             block = this % row(n, q(lo - width : lo - 1), &
                  &             qolder=q(lo - 2 * width : lo - width - 1))
          else
             block = this % row(n, q(lo - width : lo - 1))
          end if

          at = field('instant', cells, &
               & ncomp = width / max(this % space % num_vertices(), 1))
          call at % set_real_vector(q(lo : hi))

          call block % apply(this % space, [at], answer)
          call answer % get_real_vector(rn)

          y(lo : hi) = rn

       end do

    end if

    instants = support(GRAPH_SIDE_VERTEX, [(v, v = 1, ninstants)])
    out = field('chain residual', instants, ncomp=width)
    call out % set_real_vector(y)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine chain_apply

end module operation_step
