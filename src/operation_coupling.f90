!=====================================================================!
! The coupling a weight action reads: a connectivity_graph, and on it
! the two fields every such action takes - the step at each vertex,
! read off the graph's own edge degree labels. A family, a scheme
! weight and a step scaling all read exactly this tuple, so it is
! built here and nowhere else; weights_of applies any of them to it
! and reads the result out as a vector.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_coupling

  use util_precision      , only : dp
  use field_stored        , only : stored_field
  use view_directed_connectivity, only : connectivity_graph
  use operation_action    , only : operation, applied
  use operation_edge_function, only : edge_function
  use util_derivative_terms, only : derivative_terms, coefficient

  implicit none

  private
  public :: coupling_inputs, weights_of, weights_terms

contains

  !===================================================================!
  ! The step at every vertex and the two degree labels every edge
  ! already stores, read out as the tuple an edge_function's
  ! contract expects.
  !===================================================================!

  subroutine coupling_inputs(reach, steps, inputs)

    type(connectivity_graph)       , intent(in)  :: reach
    real(dp)                       , intent(in)  :: steps(:)
    type(stored_field), allocatable, intent(out) :: inputs(:)

    integer :: e, ne

    ne = reach % num_edges()

    allocate(inputs(3))
    inputs(1) = stored_field('dt'          , reach % vertex_set(), reach % num_vertices())
    inputs(2) = stored_field('tail degree' , reach % edge_set()  , ne)
    inputs(3) = stored_field('head degree' , reach % edge_set()  , ne)
    call inputs(1) % set_real_vector(steps)
    call inputs(2) % set_integer_vector([(reach % tail_degree(e), e = 1, ne)])
    call inputs(3) % set_integer_vector([(reach % head_degree(e), e = 1, ne)])

  end subroutine coupling_inputs

  !===================================================================!
  ! An action applied to that coupling, its result read out.
  !===================================================================!

  subroutine weights_of(action, reach, steps, w)

    class(operation)         , intent(in)  :: action
    type(connectivity_graph) , intent(in)  :: reach
    real(dp)                 , intent(in)  :: steps(:)
    real(dp), allocatable    , intent(out) :: w(:)

    type(stored_field), allocatable :: inputs(:)

    call coupling_inputs(reach, steps, inputs)
    call applied(action, reach, inputs, w)

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

  subroutine weights_terms(action, reach, steps, seeds, table)

    class(operation)         , intent(in)  :: action
    type(connectivity_graph) , intent(in)  :: reach
    real(dp)                 , intent(in)  :: steps(:), seeds(:,:)
    real(dp), allocatable    , intent(out) :: table(:,:)

    type(derivative_terms), allocatable :: dt(:)
    type(derivative_terms) :: c
    integer :: n, k, m, e, ne, nv

    nv = reach % num_vertices()
    ne = reach % num_edges()

    n = 0
    do while (2**n - 1 < size(seeds, 2))
       n = n + 1
    end do
    if (size(seeds, 1) /= nv .or. size(seeds, 2) /= 2**n - 1) then
       error stop 'operation_coupling: one seed row per vertex, one column per nonempty subset'
    end if

    allocate(dt(nv))
    do k = 1, nv
       dt(k) = derivative_terms(steps(k), n)
       do m = 1, 2**n - 1
          call dt(k) % set_coefficient(m, seeds(k, m))
       end do
    end do

    allocate(table(ne, 0:2**n - 1))
    select type (action)
    class is (edge_function)
       do e = 1, ne
          c = action % edge_coefficient(dt, reach % edge_tail(e), reach % edge_head(e), &
               & reach % tail_degree(e), reach % head_degree(e))
          do m = 0, 2**n - 1
             table(e, m) = coefficient(c, m)
          end do
       end do
    class default
       error stop 'operation_coupling: the weights are an edge function of the steps'
    end select

  end subroutine weights_terms

end module operation_coupling
