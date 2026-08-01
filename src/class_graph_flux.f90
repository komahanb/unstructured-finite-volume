!=====================================================================!
! A concrete face term - an edge field operation.
!
! This is where the physics between two cells lives:
!
!                          F_e
!            (i) --------------------> (j)
!            q_i                        q_j
!
! One type carrying a rule rather than a class per law. Diffusion and
! advection are the same shape of thing - a number worked out per face
! from the values on its two sides - and holding them in one type is
! what lets a balance keep a list of them.
!
!      diffusion    F = -k (q_j - q_i) / d
!      advection    F = u q, taken from whichever side is upwind
!
! At a wall there is no cell on the far side:
!
!            (i) ----------------o
!            q_i                  the boundary value stands in for q_j
!
! so a wall face is not a special kind of operation. It is this
! operation, aimed at the boundary edges, using the boundary value
! where the missing neighbour would have been.
!
! The face distance and the boundary value ride on the graph and are
! fetched by name, once, before the loop starts. Nothing here reaches
! for graph data per face.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_flux

  use iso_fortran_env     , only : dp => REAL64
  use abstract_graph_types, only : graph_edge_field_operation, graph
  use abstract_graph_types, only : graph_data, graph_edge_field, graph_edge_support
  use class_graph_support , only : edge_support
  use class_graph_field   , only : vertex_field, edge_field

  implicit none

  private
  public :: flux, FLUX_DIFFUSION, FLUX_ADVECTION

  integer, parameter :: FLUX_DIFFUSION = 1
  integer, parameter :: FLUX_ADVECTION = 2

  !===================================================================!
  ! One face term, carrying which law it follows.
  !===================================================================!

  type, extends(graph_edge_field_operation) :: flux

     integer  :: rule        = FLUX_DIFFUSION
     real(dp) :: coefficient = 1.0_dp   ! k for diffusion, u for advection
     real(dp) :: distance    = 1.0_dp   ! used when the graph carries none
     real(dp) :: boundary    = 0.0_dp   ! what stands in beyond a wall

   contains

     procedure :: name    => f_name
     procedure :: support => f_support
     procedure :: apply   => f_apply

  end type flux

  interface flux
     module procedure create
  end interface flux

contains

  pure type(flux) function create(rule, coefficient, distance, boundary) result(this)

    integer , intent(in)           :: rule
    real(dp), intent(in), optional :: coefficient
    real(dp), intent(in), optional :: distance
    real(dp), intent(in), optional :: boundary

    this % rule = rule

    if (present(coefficient)) this % coefficient = coefficient
    if (present(distance))    this % distance    = distance
    if (present(boundary))    this % boundary    = boundary

  end function create

  pure function f_name(this) result(name)

    class(flux), intent(in)       :: this
    character(len=:), allocatable :: name

    if (this % rule == FLUX_ADVECTION) then
       name = 'advection flux'
    else
       name = 'diffusion flux'
    end if

  end function f_name

  !===================================================================!
  ! Every face, walls included. A caller wanting only the inside of
  ! the mesh points a second instance at interior_edges instead; the
  ! type does not change, only the set it is aimed at.
  !===================================================================!

  subroutine f_support(this, input_graph, support)

    class(flux) , intent(in)                            :: this
    class(graph), intent(in)                            :: input_graph
    class(graph_edge_support), allocatable, intent(out) :: support

    call input_graph % all_edges(support)

  end subroutine f_support

  !===================================================================!
  ! Work out the number on every face.
  !
  ! The state arrives as the first thing in input_data. Its values are
  ! pulled out once, before the loop - not fetched per face.
  !===================================================================!

  subroutine f_apply(this, input_graph, input_data, output)

    class(flux)      , intent(in)                          :: this
    class(graph)     , intent(in)                          :: input_graph
    class(graph_data), intent(in), optional                :: input_data(:)
    class(graph_edge_field), allocatable, intent(inout)    :: output

    type(edge_field)      :: out
    type(edge_support)    :: on
    real(dp), allocatable :: q(:), fv(:)
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

    allocate(fv(ne))
    fv = 0.0_dp

    ! No state, no flux. An operation asked to work with nothing hands
    ! back zeros rather than reaching into memory it was never given.
    if (.not. present(input_data)) then
       call out % set_real_vector(fv)
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

       ! Beyond a wall there is no cell, so the boundary value stands
       ! in for the neighbour that is not there.
       if (input_graph % edge_has_head(e)) then
          h = input_graph % edge_head(e)
          if (h < 1 .or. h > size(q)) cycle
          qh = q(h)
       else
          qh = this % boundary
       end if

       select case (this % rule)

       case (FLUX_ADVECTION)
          ! Take the value from whichever side the flow is coming from.
          if (this % coefficient >= 0.0_dp) then
             fv(e) = this % coefficient * qt
          else
             fv(e) = this % coefficient * qh
          end if

       case default
          fv(e) = -this % coefficient * (qh - qt) / this % distance

       end select

    end do

    call out % set_real_vector(fv)

    ! Overwrite, by the law. A caller may lend a buffer; it does not
    ! get added to.
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine f_apply

end module class_graph_flux
