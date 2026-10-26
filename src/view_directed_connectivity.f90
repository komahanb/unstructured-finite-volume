!=====================================================================!
! Module that states the connectivity graph: a directed graph whose
! edges also store the two differential orders a discretization reads
! across.
!
! A family's block_connectivity and stage_connectivity each return one
! of these. Its vertices are instants or stages, and an edge from tail
! to head states that the row at head, of degree head_degree(edge),
! reads a value at tail, of degree tail_degree(edge):
!
!   tail: instant k-1                         head: instant k
!
!     q'(k-1)  o ------------ edge e ------------> o  q(k)
!
!   tail_degree(e) = 1                        head_degree(e) = 0
!   the row of q(k) reads q'(k-1)
!
! No arithmetic is written here. edge_coefficient (operation_family)
! and weights_terms (operation_coupling) read this incidence and these
! two degree labels to compute an edge's number. This type states only
! which edge exists and at which two degrees. tail_degree and
! head_degree never combine two edges of the same tail and head into
! one, since a discretization never states the same relation between
! the same two vertices twice.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_directed_connectivity

  use view_directed_stored, only : stored_directed_graph

  implicit none

  private

  public :: connectivity_graph

  !-------------------------------------------------------------------!
  ! Type that models a stored directed graph with a pair of degrees on
  ! every edge
  !-------------------------------------------------------------------!

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
  ! Creates the graph of nv vertices and the given edges, with one pair
  ! of degrees per edge.
  !
  ! Fortran has no call to a parent's constructor, so the parent part
  ! is built and assigned whole, and the two degree labels are attached
  ! after it.
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

  !===================================================================!
  ! Returns the degree of the value an edge reads at its tail
  !===================================================================!

  pure integer function tail_degree(this, edge_index)

    class(connectivity_graph), intent(in) :: this
    integer                  , intent(in) :: edge_index

    tail_degree = this % tail_degrees(edge_index)

  end function tail_degree

  !===================================================================!
  ! Returns the degree of the row an edge enters at its head
  !===================================================================!

  pure integer function head_degree(this, edge_index)

    class(connectivity_graph), intent(in) :: this
    integer                  , intent(in) :: edge_index

    head_degree = this % head_degrees(edge_index)

  end function head_degree

end module view_directed_connectivity
