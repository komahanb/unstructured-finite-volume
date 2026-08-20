!=====================================================================!
! Walks over a graph, written as operations.
!
! A walk reads the structure of a graph and assigns a whole number to
! every cell. That makes every walk a vertex field operation: the
! graph supplies the structure, the walk supplies the rule, and an
! integer field comes back.
!
!      graph structure  ---> walk --->  an integer per cell
!
! The graph holds no algorithms. It answers structural queries;
! algorithms act on it from outside. This separation keeps the graph
! contract small, and a new algorithm arrives without touching the
! graph.
!
! One type holds the rule, so a caller can hold several walks in a
! plain array and a new walk costs a case rather than a class.
!
!=====================================================================!
!
!                        WHAT EACH ONE ANSWERS
!
! COLOURING gives every cell a colour such that no face has the same
! colour at both ends:
!
!            (1)---(2)---(3)---(4)
!             1     2     1     2
!
! This property makes a Gauss-Seidel sweep safe to run in parallel:
! all cells of one colour update at the same time, because no two of
! them are neighbours.
!
! VISIT ORDER numbers the cells in the order a breadth-first walk
! reaches them, starting from cell one. Solvers that want to march
! through a mesh in a sensible order read this.
!
! COMPONENT gives all the cells that can reach each other the same
! number. Two cells share a number exactly when a path joins them, so
! a mesh in two pieces reports itself as two pieces.
!
! DEPTH counts faces from the seed cell. Cells the seed cannot reach
! are marked minus one, because no distance to them exists.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_walk

  use operation_action, only : operation
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use field_stored  , only : stored_field

  implicit none

  private
  public :: walk
  public :: WALK_COLOURING, WALK_VISIT_ORDER, WALK_COMPONENT, WALK_DEPTH

  integer, parameter :: WALK_COLOURING   = 1
  integer, parameter :: WALK_VISIT_ORDER = 2
  integer, parameter :: WALK_COMPONENT   = 3
  integer, parameter :: WALK_DEPTH       = 4

  !===================================================================!
  ! One walk, holding which question it answers and where it starts.
  !===================================================================!

  type, extends(operation) :: walk

     integer :: rule = WALK_COLOURING
     integer :: seed = 1

   contains

     procedure :: name   => walk_name
     procedure :: domain => walk_domain
     procedure :: apply  => walk_apply

  end type walk

  interface walk
     module procedure create
  end interface walk

contains

  !===================================================================!
  ! Build a walk that follows one rule. The seed names the vertex a
  ! depth walk starts from; the other rules need no seed.
  !===================================================================!

  type(walk) function create(rule, seed) result(this)

    integer, intent(in)           :: rule
    integer, intent(in), optional :: seed

    this % rule = rule

    if (present(seed)) this % seed = seed

    ! a walk reads no input: every answer comes from the graph
    call this % declare_arguments(0)

  end function create

  !===================================================================!
  ! The walk's name is its rule's name.
  !===================================================================!

  pure function walk_name(this) result(name)

    class(walk), intent(in)       :: this
    character(len=:), allocatable :: name

    select case (this % rule)
    case (WALK_VISIT_ORDER)
       name = 'visit order'
    case (WALK_COMPONENT)
       name = 'component'
    case (WALK_DEPTH)
       name = 'depth'
    case default
       name = 'colouring'
    end select

  end function walk_name

  !===================================================================!
  ! Where the answer lives: one value per vertex, whatever the rule.
  !
  !      colouring      a colour per vertex
  !      visit order    a position per vertex
  !      component      a component number per vertex
  !      depth          a distance per vertex
  !===================================================================!

  subroutine walk_domain(this, input_graph, domain, num_entries)

    class(walk) , intent(in)               :: this
    class(directed_graph), intent(in)               :: input_graph
    type(graph), intent(out) :: domain
    integer        , intent(out) :: num_entries

    associate (u1 => this); end associate

    domain   = input_graph % all_vertices()
    num_entries = input_graph % num_vertices()

  end subroutine walk_domain

  !===================================================================!
  ! Walk the graph and return a whole number per cell.
  !
  ! Nothing here reads input_data. Every answer comes from the shape
  ! of the graph alone: structure in, rule applied, integers out.
  !===================================================================!

  subroutine walk_apply(this, input_graph, input_data, output)

    class(walk)       , intent(in)                 :: this
    class(directed_graph)      , intent(in)                 :: input_graph
    class(field), intent(in), optional       :: input_data(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)           :: out
    integer , allocatable :: mark(:)
    integer :: nv

    associate (u1 => present(input_data)); end associate

    nv = input_graph % num_vertices()

    out = stored_field(this % name(), input_graph % vertex_set(), input_graph % num_vertices())

    select case (this % rule)
    case (WALK_VISIT_ORDER)
       call breadth_first(input_graph, this % seed, mark, want_depth=.false.)
    case (WALK_DEPTH)
       call breadth_first(input_graph, this % seed, mark, want_depth=.true.)
    case (WALK_COMPONENT)
       call components(input_graph, mark)
    case default
       call colour(input_graph, mark)
    end select

    call out % set_integer_vector(mark)

    ! A supplied buffer is overwritten, never added to.
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine walk_apply

  !===================================================================!
  ! Give every cell the lowest colour none of its neighbours has taken.
  !
  ! Greedy, so the colour count is not minimal. The guarantee that
  ! matters holds: no face has one colour at both ends, which is what
  ! makes a colour safe to sweep in parallel.
  !===================================================================!

  subroutine colour(input_graph, mark)

    class(directed_graph)        , intent(in)  :: input_graph
    integer, allocatable, intent(out) :: mark(:)

    integer, allocatable :: nbrs(:)
    logical, allocatable :: taken(:)
    integer :: nv, v, i, c

    nv = input_graph % num_vertices()
    allocate(mark(nv))
    mark = 0

    allocate(taken(nv + 1))

    do v = 1, nv

       taken = .false.
       call input_graph % adjacent_vertices(v, nbrs)
       do i = 1, size(nbrs)
          if (mark(nbrs(i)) >= 1) taken(mark(nbrs(i))) = .true.
       end do

       ! The lowest colour nobody next door is wearing.
       c = 1
       do while (c <= nv .and. taken(c))
          c = c + 1
       end do
       mark(v) = c

    end do

  end subroutine colour

  !===================================================================!
  ! Walk outward from the seed, one ring at a time. Either number the
  ! cells in the order they are reached, or count the faces crossed to
  ! reach them.
  !
  ! A cell the seed cannot reach gets minus one either way. Saying
  ! "not reachable" beats reporting a distance of zero, which would be
  ! indistinguishable from being the seed.
  !===================================================================!

  subroutine breadth_first(input_graph, seed, mark, want_depth)

    class(directed_graph)        , intent(in)  :: input_graph
    integer             , intent(in)  :: seed
    integer, allocatable, intent(out) :: mark(:)
    logical             , intent(in)  :: want_depth

    integer, allocatable :: queue(:), depth(:), nbrs(:)
    integer :: nv, head_of_queue, tail_of_queue, v, i, rank

    nv = input_graph % num_vertices()
    allocate(mark(nv), queue(nv), depth(nv))
    mark  = -1
    depth = -1

    if (seed < 1 .or. seed > nv) return

    queue(1)      = seed
    head_of_queue = 1
    tail_of_queue = 1
    depth(seed)   = 0
    rank          = 1
    if (want_depth) then
       mark(seed) = 0
    else
       mark(seed) = rank
    end if

    do while (head_of_queue <= tail_of_queue)

       v = queue(head_of_queue)
       head_of_queue = head_of_queue + 1

       call input_graph % adjacent_vertices(v, nbrs)
       do i = 1, size(nbrs)
          if (depth(nbrs(i)) >= 0) cycle
          depth(nbrs(i))       = depth(v) + 1
          rank                 = rank + 1
          tail_of_queue        = tail_of_queue + 1
          queue(tail_of_queue) = nbrs(i)
          if (want_depth) then
             mark(nbrs(i)) = depth(nbrs(i))
          else
             mark(nbrs(i)) = rank
          end if
       end do

    end do

  end subroutine breadth_first

  !===================================================================!
  ! Give the same number to every cell that can reach every other.
  !
  ! Start a fresh number at the first unmarked cell, flood as far as
  ! it goes, and repeat. Two cells share a number exactly when a path
  ! joins them.
  !===================================================================!

  subroutine components(input_graph, mark)

    class(directed_graph)        , intent(in)  :: input_graph
    integer, allocatable, intent(out) :: mark(:)

    integer, allocatable :: queue(:), nbrs(:)
    integer :: nv, v, i, which, head_of_queue, tail_of_queue, u

    nv = input_graph % num_vertices()
    allocate(mark(nv), queue(nv))
    mark  = 0
    which = 0

    do v = 1, nv

       if (mark(v) /= 0) cycle

       which         = which + 1
       mark(v)       = which
       queue(1)      = v
       head_of_queue = 1
       tail_of_queue = 1

       do while (head_of_queue <= tail_of_queue)
          u = queue(head_of_queue)
          head_of_queue = head_of_queue + 1
          call input_graph % adjacent_vertices(u, nbrs)
          do i = 1, size(nbrs)
             if (mark(nbrs(i)) /= 0) cycle
             mark(nbrs(i))        = which
             tail_of_queue        = tail_of_queue + 1
             queue(tail_of_queue) = nbrs(i)
          end do
       end do

    end do

  end subroutine components

end module operation_walk
