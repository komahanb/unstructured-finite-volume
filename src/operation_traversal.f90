!=====================================================================!
! Traversals over a graph, written as operations.
!
! A traversal reads the structure of a graph and assigns a whole
! number to every cell. That makes every traversal a vertex field
! operation: the graph supplies the structure, the traversal supplies
! the rule, and an integer field is returned.
!
!      graph structure  ---> traversal --->  an integer per cell
!
! The graph stores no algorithms. The graph evaluates structural
! queries; algorithms are applied to it from outside. This separation
! keeps the graph contract small, and a new algorithm is added without
! changing the graph.
!
! One type stores the rule, so a caller can store several traversals
! in a plain array and a new traversal costs a case rather than a
! class.
!
!=====================================================================!
!
!                        WHAT EACH RULE COMPUTES
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
! VISIT ORDER numbers the cells in the order a breadth-first
! traversal reaches them, starting from cell one. Solvers that require
! a consistent ordering of the mesh read this.
!
! COMPONENT gives all the cells that can reach each other the same
! number. Two cells share a number exactly when a path joins them, so
! a mesh with two connected components reports two components.
!
! DEPTH counts faces from the seed cell. Cells the seed cannot reach
! are marked minus one, because no distance to them exists.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_traversal

  use operation_action, only : operation, binding
  use operation_action, only : emit
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use field_stored  , only : stored_field

  implicit none

  private
  public :: traversal
  public :: TRAVERSAL_COLOURING, TRAVERSAL_VISIT_ORDER, TRAVERSAL_COMPONENT, TRAVERSAL_DEPTH

  integer, parameter :: TRAVERSAL_COLOURING   = 1
  integer, parameter :: TRAVERSAL_VISIT_ORDER = 2
  integer, parameter :: TRAVERSAL_COMPONENT   = 3
  integer, parameter :: TRAVERSAL_DEPTH       = 4

  !===================================================================!
  ! One traversal, storing which rule it computes and where it starts.
  !===================================================================!

  type, extends(operation) :: traversal

     integer :: rule = TRAVERSAL_COLOURING
     integer :: seed = 1

   contains

     procedure :: apply  => traversal_apply

  end type traversal

  interface traversal
     module procedure create
  end interface traversal

contains

  !===================================================================!
  ! Construct a traversal that follows one rule, named after it. The
  ! seed names the vertex a depth traversal starts from; the other
  ! rules need no seed.
  !===================================================================!

  type(traversal) function create(rule, seed) result(this)

    integer, intent(in)           :: rule
    integer, intent(in), optional :: seed

    this % rule = rule

    if (present(seed)) this % seed = seed

    ! a traversal reads no input: every value comes from the graph
    select case (rule)
    case (TRAVERSAL_VISIT_ORDER)
       call this % declare_arguments(0, label='visit order')
    case (TRAVERSAL_COMPONENT)
       call this % declare_arguments(0, label='component')
    case (TRAVERSAL_DEPTH)
       call this % declare_arguments(0, label='depth')
    case default
       call this % declare_arguments(0, label='colouring')
    end select

  end function create

  !===================================================================!
  ! Traverse the graph and return a whole number per cell.
  !
  ! Nothing here reads inputs. Every value comes from the structure
  ! of the graph alone: structure in, rule applied, integers out.
  !===================================================================!

  subroutine traversal_apply(this, input_graph, inputs, output)

    class(traversal)       , intent(in)                 :: this
    class(directed_graph)      , intent(in)                 :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    type(stored_field)           :: out
    integer , allocatable :: mark(:)
    integer :: nv

    associate (u1 => present(inputs)); end associate

    nv = input_graph % num_vertices()

    out = stored_field(this % name(), input_graph % vertex_set(), input_graph % num_vertices())

    select case (this % rule)
    case (TRAVERSAL_VISIT_ORDER)
       call breadth_first(input_graph, this % seed, mark, depth_mode=.false.)
    case (TRAVERSAL_DEPTH)
       call breadth_first(input_graph, this % seed, mark, depth_mode=.true.)
    case (TRAVERSAL_COMPONENT)
       call components(input_graph, mark)
    case default
       call colour(input_graph, mark)
    end select

    call out % set_integer_vector(mark)

    call emit(out, output)

  end subroutine traversal_apply

  !===================================================================!
  ! Assign every cell the lowest colour not assigned to any neighbour.
  !
  ! Greedy, so the colour count is not minimal. The required guarantee
  ! is satisfied: no face has one colour at both ends, which is what makes a
  ! colour safe to sweep in parallel.
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

       ! The lowest colour not assigned to any neighbour.
       c = 1
       do while (c <= nv .and. taken(c))
          c = c + 1
       end do
       mark(v) = c

    end do

  end subroutine colour

  !===================================================================!
  ! Traverse outward from the seed, one level at a time. Either number
  ! the cells in the order they are reached, or count the faces crossed
  ! to reach them.
  !
  ! A cell the seed cannot reach is marked minus one in either mode.
  ! Marking "not reachable" is preferred to a distance of zero, which
  ! would be indistinguishable from the seed.
  !===================================================================!

  subroutine breadth_first(input_graph, seed, mark, depth_mode)

    class(directed_graph)        , intent(in)  :: input_graph
    integer             , intent(in)  :: seed
    integer, allocatable, intent(out) :: mark(:)
    logical             , intent(in)  :: depth_mode

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
    if (depth_mode) then
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
          if (depth_mode) then
             mark(nbrs(i)) = depth(nbrs(i))
          else
             mark(nbrs(i)) = rank
          end if
       end do

    end do

  end subroutine breadth_first

  !===================================================================!
  ! Assign the same number to every cell of one connected component.
  !
  ! Start a new number at the first unmarked cell, traverse the whole
  ! component, and repeat. Two cells share a number exactly when a path
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

end module operation_traversal
