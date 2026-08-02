!=====================================================================!
! An edge derivative - a concrete edge field operation.
!
! The operation computes one number per edge from the values at the
! edge's two ends:
!
!                          z_e
!            (i) --------------------> (j)
!            q_i                        q_j
!
! The number is a sample of the vertex data, and the order of the
! derivative it produces is set by which sample is taken:
!
!    order 0     c * (q_i + q_j) / 2         the average of the ends
!    order 1     c * q_i  or  c * q_j        one end, chosen by the
!                                            sign of c
!    order 2     c * (q_j - q_i) / h         the difference quotient,
!                                            h the spacing
!
! The derivative appears when a balance folds these samples through
! incidence. Folding the order-1 sample over all edges of a vertex
! produces c times the first derivative of q there; folding the
! order-2 sample produces c times the second derivative. The order-0
! sample is interpolation and produces no derivative; it supplies the
! undifferentiated term of an equation.
!
! This layer is mathematics only. The names a physical model gives
! these operators - and its sign conventions - belong to classes that
! extend or hold this one, outside the graph layer. A transport model,
! for example, uses order 1 with c set to a velocity and order 2 with
! c set to the negative of a conductivity; both choices are the
! model's, not this type's.
!
! An edge with no head has one end. The stored boundary value stands
! in for the value at the missing end:
!
!            (i) ----------------o
!            q_i                  boundary_value
!
! The spacing and the boundary value are stored on the operation. A
! graph that holds per-edge spacing as data by name supersedes the
! stored constant; that path arrives with the geometry work.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_derivative

  use iso_fortran_env     , only : dp => REAL64
  use abstract_graph_types, only : graph_edge_field_operation, graph
  use abstract_graph_types, only : graph_data, graph_edge_field, graph_edge_support
  use class_graph_support , only : edge_support
  use class_graph_field   , only : vertex_field, edge_field

  implicit none

  private
  public :: derivative

  !===================================================================!
  ! One edge derivative of a given order.
  !===================================================================!

  type, extends(graph_edge_field_operation) :: derivative

     integer  :: order          = 2
     real(dp) :: coefficient    = 1.0_dp
     real(dp) :: spacing        = 1.0_dp
     real(dp) :: boundary_value = 0.0_dp

   contains

     procedure :: name    => derivative_name
     procedure :: support => derivative_support
     procedure :: apply   => derivative_apply

  end type derivative

  interface derivative
     module procedure create
  end interface derivative

contains

  pure type(derivative) function create(order, coefficient, spacing, boundary_value) result(this)

    integer , intent(in)           :: order
    real(dp), intent(in), optional :: coefficient
    real(dp), intent(in), optional :: spacing
    real(dp), intent(in), optional :: boundary_value

    this % order = order

    if (present(coefficient))    this % coefficient    = coefficient
    if (present(spacing))        this % spacing        = spacing
    if (present(boundary_value)) this % boundary_value = boundary_value

  end function create

  pure function derivative_name(this) result(name)

    class(derivative), intent(in) :: this
    character(len=:), allocatable :: name

    select case (this % order)
    case (0)
       name = 'zeroth derivative'
    case (1)
       name = 'first derivative'
    case default
       name = 'second derivative'
    end select

  end function derivative_name

  !===================================================================!
  ! Every edge. A caller that requires a subset points a second
  ! instance at that subset; the type does not change, only the set
  ! it acts on.
  !===================================================================!

  subroutine derivative_support(this, input_graph, support)

    class(derivative), intent(in)                       :: this
    class(graph)     , intent(in)                       :: input_graph
    class(graph_edge_support), allocatable, intent(out) :: support

    call input_graph % all_edges(support)

  end subroutine derivative_support

  !===================================================================!
  ! Compute the sample on every edge.
  !
  ! The vertex data arrives as the first entry of input_data. Its
  ! values are fetched once, before the loop - never per edge.
  !===================================================================!

  subroutine derivative_apply(this, input_graph, input_data, output)

    class(derivative), intent(in)                       :: this
    class(graph)     , intent(in)                       :: input_graph
    class(graph_data), intent(in), optional             :: input_data(:)
    class(graph_edge_field), allocatable, intent(inout) :: output

    type(edge_field)      :: out
    type(edge_support)    :: on
    real(dp), allocatable :: q(:), z(:)
    integer , allocatable :: indices(:)
    real(dp)              :: qt, qh
    integer               :: ne, e, t, h

    ne = input_graph % num_edges()

    allocate(indices(ne))
    do e = 1, ne
       indices(e) = e
    end do

    on  = edge_support(indices)
    out = edge_field(this % name(), on)

    allocate(z(ne))
    z = 0.0_dp

    ! No vertex data, no derivative. Called with nothing, this returns
    ! zeros rather than reading memory it was never given.
    if (.not. present(input_data)) then
       call out % set_real_vector(z)
       if (allocated(output)) deallocate(output)
       allocate(output, source=out)
       return
    end if

    select type (state => input_data(1))
    class is (vertex_field)
       call state % get_real_vector(q)
    class default
       allocate(q(0))
    end select

    do e = 1, ne

       t = input_graph % edge_tail(e)
       if (t < 1 .or. t > size(q)) cycle
       qt = q(t)

       ! At an edge with no head, the stored boundary value stands in
       ! for the value at the missing end.
       if (input_graph % edge_has_head(e)) then
          h = input_graph % edge_head(e)
          if (h < 1 .or. h > size(q)) cycle
          qh = q(h)
       else
          qh = this % boundary_value
       end if

       select case (this % order)

       case (0)
          z(e) = this % coefficient * 0.5_dp * (qt + qh)

       case (1)
          ! One end, chosen by the sign of the coefficient. Folded
          ! through incidence this sample produces a first derivative,
          ! and this choice of end keeps the fold stable for either
          ! sign.
          if (this % coefficient >= 0.0_dp) then
             z(e) = this % coefficient * qt
          else
             z(e) = this % coefficient * qh
          end if

       case default
          ! The difference quotient. Folded through incidence this
          ! sample produces a second derivative.
          z(e) = this % coefficient * (qh - qt) / this % spacing

       end select

    end do

    call out % set_real_vector(z)

    ! A supplied buffer is overwritten, never added to.
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine derivative_apply

end module class_graph_derivative
