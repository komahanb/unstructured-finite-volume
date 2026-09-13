!=====================================================================!
! The connectivity graph: a directed graph whose edges also store the
! two differential orders a discretization reads across.
!
! A family's block_connectivity and stage_connectivity each return one of these:
! vertices are instants or stages, and an edge from tail to head
! states that the row at head, of degree head_degree(edge), reads a
! value at tail, of degree tail_degree(edge). No arithmetic is
! written here - edge_coefficient (operation_family) and weights_terms
! (operation_coupling) read this incidence and these two degree
! labels to compute an edge's number; this type states only which
! edge exists and at which two degrees, stacked by the same ordering
! degree itself imposes: tail_degree and head_degree never combine
! two edges of the same tail and head into one, since a discretization
! never states the same relation between the same two vertices twice.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_directed_connectivity

  use view_directed_stored, only : stored_directed_graph

  implicit none

  private
  public :: connectivity_graph

  type, extends(stored_directed_graph) :: connectivity_graph

     integer, allocatable, private :: tail_degrees(:)
     integer, allocatable, private :: head_degrees(:)

   contains

     procedure :: tail_degree
     procedure :: head_degree

  end type connectivity_graph

  interface connectivity_graph
     module procedure create
  end interface connectivity_graph

contains

  !===================================================================!
  ! Build the incidence through the parent constructor - Fortran has
  ! no super-constructor call, so the parent part is built and
  ! assigned whole - then attach the two degree labels, one pair per
  ! edge.
  !===================================================================!

  type(connectivity_graph) function create(nv, tails, heads, tail_degrees, head_degrees) &
       & result(this)

    integer, intent(in) :: nv
    integer, intent(in) :: tails(:), heads(:)
    integer, intent(in) :: tail_degrees(:), head_degrees(:)

    character(len=250) :: message

    if (size(tail_degrees) /= size(tails) .or. size(head_degrees) /= size(tails)) then
       write(message,'(a,i0,a,i0,a,i0)') 'view_directed_connectivity: one degree pair is &
            &required per edge; size(tails) = ', size(tails), ', size(tail_degrees) = ', &
            & size(tail_degrees), ', size(head_degrees) = ', size(head_degrees)
       error stop trim(message)
    end if

    this % stored_directed_graph = stored_directed_graph(nv, tails=tails, heads=heads)
    this % tail_degrees = tail_degrees
    this % head_degrees = head_degrees

  end function create

  pure integer function tail_degree(this, edge_index)

    class(connectivity_graph), intent(in) :: this
    integer                  , intent(in) :: edge_index

    tail_degree = this % tail_degrees(edge_index)

  end function tail_degree

  pure integer function head_degree(this, edge_index)

    class(connectivity_graph), intent(in) :: this
    integer                  , intent(in) :: edge_index

    head_degree = this % head_degrees(edge_index)

  end function head_degree

end module view_directed_connectivity
