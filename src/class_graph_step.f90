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

module class_graph_step

  use iso_fortran_env    , only : dp => REAL64
  use graph_grammar      , only : graph, graph_field, graph_operation
  use graph_calculus     , only : discretization_operator
  use graph_calculus     , only : GRAPH_SIDE_VERTEX
  use class_graph_support, only : support
  use class_graph_field  , only : field
  use class_graph        , only : stored_graph

  implicit none

  private
  public :: step_operator
  public :: backward_euler, bdf

  type, extends(discretization_operator) :: step_operator

     class(graph_operation), allocatable :: action

     real(dp), allocatable :: qold(:)
     real(dp), allocatable :: qolder(:)

     real(dp) :: a0 = 1.0_dp
     real(dp) :: a1 = -1.0_dp
     real(dp) :: a2 = 0.0_dp
     real(dp) :: hs = 0.0_dp

     integer :: reach = 1

   contains

     procedure :: name         => step_name
     procedure :: domain       => step_domain
     procedure :: apply        => step_apply
     procedure :: dependencies => step_dependencies

  end type step_operator

contains

  !===================================================================!
  ! The shelf. Backward euler: one instant back.
  !===================================================================!

  function backward_euler(action, h) result(this)

    class(graph_operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(step_operator)                :: this

    allocate(this % action, source=action)
    this % a0    = 1.0_dp
    this % a1    = -1.0_dp
    this % a2    = 0.0_dp
    this % hs    = h
    this % reach = 1

  end function backward_euler

  !===================================================================!
  ! The bdf family: k instants back, one name, k absorbed. Order one
  ! is backward euler by another road; order two is the second-order
  ! table; higher orders join here when their tables are ruled.
  !===================================================================!

  function bdf(k, action, h) result(this)

    integer, intent(in)                :: k
    class(graph_operation), intent(in) :: action
    real(dp), intent(in)               :: h
    type(step_operator)                :: this

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

    class(step_operator), intent(in)       :: this
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

    class(step_operator), intent(in) :: this
    character(len=:), allocatable :: name

    associate (u1 => this); end associate

    name = 'step'

  end function step_name

  subroutine step_domain(this, input_graph, domain)

    class(step_operator), intent(in)       :: this
    class(graph), intent(in)               :: input_graph
    class(graph), allocatable, intent(out) :: domain

    associate (u1 => this); end associate

    call input_graph % all_vertices(domain)

  end subroutine step_domain

  subroutine step_apply(this, input_graph, input_data, output)

    class(step_operator), intent(in)               :: this
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

       call this % action % apply(input_graph, input_data, velocity)
       call velocity % get_real_vector(s)

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

end module class_graph_step
