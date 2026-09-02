!=====================================================================!
! The coupling a weight action reads: the instants as a directed
! graph, and on it the three fields every such action takes, in this
! order - the step at each instant, and on each edge the degree of
! the source and the degree the edge's condition determines. A
! family, a scheme weight and a step scaling all read exactly this
! tuple, so it is built here and nowhere else; weights_of applies
! any of them to it and reads the result out as a vector.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_coupling

  use util_precision      , only : dp
  use field_calculus      , only : field
  use field_stored        , only : stored_field
  use view_directed_stored, only : stored_directed_graph
  use operation_action    , only : operation
  use operation_edge_function, only : edge_function
  use util_derivative_terms, only : derivative_terms, coefficient

  implicit none

  private
  public :: coupling_inputs, weights_of, weights_terms

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
  ! An action applied to that coupling, its result read out.
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
    call action % apply(edges, action % bind(inputs), out)
    call out % real_vector(w)

  end subroutine weights_of

  !===================================================================!
  ! Every coefficient of an action's weights over n directions in the
  ! steps, the steps seeded subset by subset: seeds(k, m) is the total
  ! derivative of step k along the subset with mask m, for m = 1 to
  ! 2^n - 1, and table(e, m) is the total derivative of edge e's
  ! weight along that subset, the value at m = 0. Singleton seeds
  ! alone give the mixed partials; the total derivatives of the steps
  ! along subsets of designs give the weights' total derivatives along
  ! the same, every set partition included by the product rule. The
  ! action computes the weights over derivative terms, so every number
  ! is exact. Invalid input: a seed table of other than one row per
  ! vertex and 2^n - 1 columns, or an action that is not an edge
  ! function.
  !===================================================================!

  subroutine weights_terms(action, num_vertices, tails, heads, steps, seeds, &
       & source_degree, determines, table)

    class(operation), intent(in) :: action
    integer         , intent(in) :: num_vertices
    integer         , intent(in) :: tails(:), heads(:)
    real(dp)        , intent(in) :: steps(:), seeds(:,:)
    integer         , intent(in) :: source_degree(:), determines(:)
    real(dp), allocatable, intent(out) :: table(:,:)

    type(derivative_terms), allocatable :: dt(:)
    type(derivative_terms) :: c
    integer :: n, k, m, e

    n = 0
    do while (2**n - 1 < size(seeds, 2))
       n = n + 1
    end do
    if (size(seeds, 1) /= num_vertices .or. size(seeds, 2) /= 2**n - 1) then
       error stop 'operation_coupling: one seed row per vertex, one column per nonempty subset'
    end if

    allocate(dt(num_vertices))
    do k = 1, num_vertices
       dt(k) = derivative_terms(steps(k), n)
       do m = 1, 2**n - 1
          call dt(k) % set_coefficient(m, seeds(k, m))
       end do
    end do

    allocate(table(size(tails), 0:2**n - 1))
    select type (action)
    class is (edge_function)
       do e = 1, size(tails)
          c = action % edge_coefficient(dt, tails(e), heads(e), source_degree(e), determines(e))
          do m = 0, 2**n - 1
             table(e, m) = coefficient(c, m)
          end do
       end do
    class default
       error stop 'operation_coupling: the weights are an edge function of the steps'
    end select

  end subroutine weights_terms

end module operation_coupling
