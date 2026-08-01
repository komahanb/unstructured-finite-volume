!=====================================================================!
! The balance - a vertex field operation built out of smaller ones.
!
! A balance is what a cell has left over: what its faces bring in,
! plus whatever the cell itself contributes.
!
!                          z_e
!                   (i) ---------> (j)
!
!                   y_i = y_i - z_e
!                   y_j = y_j + z_e
!
! Every face gives its number to the cell it leaves and takes it from
! the cell it enters. That is incidence, and it is the whole of the
! folding.
!
! EXACTLY ONCE. Each face is walked one time and touches its two cells
! one time. Walk a face twice and the balance is wrong by that face;
! miss one and it is wrong by that face the other way. Neither shows
! up as a crash - it shows up as a solution that is quietly not the
! solution, which is why the count is checked in the suite rather than
! trusted.
!
! A wall face has no far cell, so its number lands on the one cell it
! touches and stops there:
!
!                   (i) --------o
!
!                   y_i = y_i - z_b
!
!=====================================================================!
!
!                   WHY THIS IS NOT A NEW KIND OF THING
!
! A balance is a vertex field operation like any other. It happens to
! drive an edge operation on the way, but what it hands back is a
! value per cell, and that is all its type promises.
!
! A solver calls the result a residual. That word names a stage in a
! solve, not a thing in this library, so it does not appear here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_balance

  use iso_fortran_env     , only : dp => REAL64
  use abstract_graph_types, only : graph_vertex_field_operation, graph
  use abstract_graph_types, only : graph_data, graph_vertex_field
  use abstract_graph_types, only : graph_vertex_support, graph_edge_field
  use class_graph_support , only : vertex_support
  use class_graph_field   , only : vertex_field
  use class_graph_flux    , only : flux

  implicit none

  private
  public :: balance

  !===================================================================!
  ! A balance holds its face terms and, if it has one, a number added
  ! to every cell.
  !
  ! The face terms are one concrete type carrying a rule. A Fortran
  ! array holds a single dynamic type, so one concrete flux type is
  ! what lets the balance keep its terms in a plain array.
  !===================================================================!

  type, extends(graph_vertex_field_operation) :: balance

     type(flux), allocatable :: face_terms(:)

     real(dp) :: source = 0.0_dp

   contains

     procedure :: name    => b_name
     procedure :: support => b_support
     procedure :: apply   => b_apply

  end type balance

  interface balance
     module procedure create
  end interface balance

contains

  pure type(balance) function create(face_terms, source) result(this)

    type(flux), intent(in), optional :: face_terms(:)
    real(dp)  , intent(in), optional :: source

    if (present(face_terms)) allocate(this % face_terms, source=face_terms)
    if (present(source))     this % source = source

  end function create

  pure function b_name(this) result(name)

    class(balance), intent(in)    :: this
    character(len=:), allocatable :: name

    name = 'balance'

  end function b_name

  !===================================================================!
  ! Every cell. A balance answers for the whole graph it is given; a
  ! caller wanting only part of one hands it only part of one.
  !===================================================================!

  subroutine b_support(this, input_graph, support)

    class(balance), intent(in)                            :: this
    class(graph)  , intent(in)                            :: input_graph
    class(graph_vertex_support), allocatable, intent(out) :: support

    call input_graph % all_vertices(support)

  end subroutine b_support

  !===================================================================!
  ! Work the balance out.
  !
  !    1. start every cell at its own source term
  !    2. for each face term, work out every face
  !    3. fold each face onto the two cells it touches, once
  !===================================================================!

  subroutine b_apply(this, input_graph, input_data, output)

    class(balance)   , intent(in)                       :: this
    class(graph)     , intent(in)                       :: input_graph
    class(graph_data), intent(in), optional             :: input_data(:)
    class(graph_vertex_field), allocatable, intent(inout) :: output

    class(graph_edge_field), allocatable :: face_values

    type(vertex_field)    :: out
    type(vertex_support)  :: on
    real(dp), allocatable :: y(:), z(:)
    integer , allocatable :: indices(:)
    integer               :: nv, ne, v, e, t, h, k

    nv = input_graph % num_vertices()
    ne = input_graph % num_edges()

    allocate(indices(nv))
    do v = 1, nv
       indices(v) = v
    end do

    on  = vertex_support(indices)
    out = vertex_field('balance', on)

    allocate(y(nv))
    y = this % source

    if (allocated(this % face_terms)) then
       do k = 1, size(this % face_terms)

          ! One face term, worked out for every face at once. This is
          ! the only place the face values are computed.
          call this % face_terms(k) % apply(input_graph, input_data, face_values)
          call face_values % get_real_vector(z)

          ! And folded onto the cells, each face touching its cells
          ! exactly one time.
          do e = 1, ne
             if (e > size(z)) exit

             t = input_graph % edge_tail(e)
             if (t >= 1 .and. t <= nv) y(t) = y(t) - z(e)

             if (input_graph % edge_has_head(e)) then
                h = input_graph % edge_head(e)
                if (h >= 1 .and. h <= nv) y(h) = y(h) + z(e)
             end if

          end do

       end do
    end if

    call out % set_real_vector(y)

    ! Overwrite, by the law. The buffer may be lent; it is not added to.
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine b_apply

end module class_graph_balance
