!=====================================================================!
! The coupling a weight action reads: the instants as a directed
! graph, and on it the three fields every such action takes, in this
! order - the step at each instant, and on each edge the degree of
! the source and the degree the edge's condition determines. A
! family, a scheme weight and a step scaling all read exactly this
! tuple, so it is built here and nowhere else; weights_of applies
! any of them to it and reads the answer out as a vector.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_coupling

  use util_precision      , only : dp
  use field_calculus      , only : field
  use field_stored        , only : stored_field
  use view_directed_stored, only : stored_directed_graph
  use operation_action    , only : operation

  implicit none

  private
  public :: coupling_inputs, weights_of

contains

  !===================================================================!
  ! The graph and the three fields, filled. One step per vertex, one
  ! source degree and one determined degree per edge.
  !===================================================================!

  subroutine coupling_inputs(num_vertices, tails, heads, steps, source_degree, &
       & determines, edges, inputs)

    integer , intent(in) :: num_vertices
    integer , intent(in) :: tails(:), heads(:)
    real(dp), intent(in) :: steps(:)
    integer , intent(in) :: source_degree(:), determines(:)
    type(stored_directed_graph)    , intent(out) :: edges
    type(stored_field), allocatable, intent(out) :: inputs(:)

    edges = stored_directed_graph(num_vertices, tails=tails, heads=heads)

    allocate(inputs(3))
    inputs(1) = stored_field('dt'           , edges % vertex_set(), num_vertices)
    inputs(2) = stored_field('source degree', edges % edge_set()  , size(tails))
    inputs(3) = stored_field('determines'   , edges % edge_set()  , size(tails))
    call inputs(1) % set_real_vector(steps)
    call inputs(2) % set_integer_vector(source_degree)
    call inputs(3) % set_integer_vector(determines)

  end subroutine coupling_inputs

  !===================================================================!
  ! An action applied to that coupling, its answer read out.
  !===================================================================!

  subroutine weights_of(action, num_vertices, tails, heads, steps, source_degree, &
       & determines, w)

    class(operation), intent(in) :: action
    integer         , intent(in) :: num_vertices
    integer         , intent(in) :: tails(:), heads(:)
    real(dp)        , intent(in) :: steps(:)
    integer         , intent(in) :: source_degree(:), determines(:)
    real(dp), allocatable, intent(out) :: w(:)

    type(stored_directed_graph)     :: edges
    type(stored_field), allocatable :: inputs(:)
    class(field)      , allocatable :: out

    call coupling_inputs(num_vertices, tails, heads, steps, source_degree, determines, &
         & edges, inputs)
    call action % apply(edges, inputs, out)
    call out % real_vector(w)

  end subroutine weights_of

end module operation_coupling
