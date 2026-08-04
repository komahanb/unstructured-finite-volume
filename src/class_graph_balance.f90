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
! reduction through incidence.
!
! EXACTLY ONCE. Each face is walked one time and touches its two cells
! one time. Walk a face twice and the balance is wrong by that face;
! miss one and it is wrong by that face the other way. Neither shows
! up as a crash - it shows up as a solution that is quietly not the
! solution, so the incidence count is verified in the test suite
! rather than assumed.
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
! drive an edge operation on the way, but what it returns is a
! value per cell, which is exactly what its type declares.
!
! A solver calls the result a residual. That word names a stage in a
! solve, not a thing in this library, so it does not appear here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_balance

  use iso_fortran_env    , only : dp => REAL64
  use graph_grammar      , only : graph_operation, graph, graph_field
  use graph_calculus     , only : GRAPH_SIDE_VERTEX
  use class_graph_support, only : support
  use class_graph_field  , only : field
  use class_graph_differential_operator, only : differential_operator

  implicit none

  private
  public :: balance

  !===================================================================!
  ! A balance holds its face terms and, if it has one, a number added
  ! to every cell.
  !
  ! The face terms are one concrete type holding a rule. A Fortran
  ! array holds a single dynamic type, so one concrete edge-derivative
  ! type is what lets the balance keep its terms in a plain array.
  !===================================================================!

  type, extends(graph_operation) :: balance

     type(differential_operator), allocatable :: edge_terms(:)

     real(dp) :: source = 0.0_dp

   contains

     procedure :: name   => balance_name
     procedure :: domain => balance_domain
     procedure :: apply  => balance_apply

  end type balance

  interface balance
     module procedure create
  end interface balance

contains

  !===================================================================!
  ! Build a balance from what a model declares: the face terms, and
  ! one number per cell added as a source. Both optional - an empty
  ! balance answers zero.
  !===================================================================!

  pure type(balance) function create(edge_terms, source) result(this)

    type(differential_operator), intent(in), optional :: edge_terms(:)
    real(dp)  , intent(in), optional :: source

    if (present(edge_terms)) allocate(this % edge_terms, source=edge_terms)
    if (present(source))     this % source = source

  end function create

  !===================================================================!
  ! The operation's name, for reports.
  !===================================================================!

  pure function balance_name(this) result(name)

    class(balance), intent(in)    :: this
    character(len=:), allocatable :: name

    name = 'balance'

  end function balance_name

  !===================================================================!
  ! Every cell. A balance answers for the whole graph it is given; a
  ! caller wanting only part of one hands it only part of one.
  !===================================================================!

  subroutine balance_domain(this, input_graph, domain)

    class(balance), intent(in)             :: this
    class(graph)  , intent(in)             :: input_graph
    class(graph), allocatable, intent(out) :: domain

    call input_graph % all_vertices(domain)

  end subroutine balance_domain

  !===================================================================!
  ! Work the balance out.
  !
  !    1. start every cell at its own source term
  !    2. for each face term, work out every face
  !    3. reduce each edge onto the two vertices it touches, once,
  !       through incidence
  !===================================================================!

  subroutine balance_apply(this, input_graph, input_data, output)

    class(balance)    , intent(in)                 :: this
    class(graph)      , intent(in)                 :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    class(graph_field), allocatable :: edge_values

    type(field)           :: out
    type(support)         :: on
    real(dp), allocatable :: y(:), z(:)
    integer , allocatable :: indices(:)
    integer               :: nv, ne, v, e, t, h, k

    nv = input_graph % num_vertices()
    ne = input_graph % num_edges()

    allocate(indices(nv))
    do v = 1, nv
       indices(v) = v
    end do

    on  = support(GRAPH_SIDE_VERTEX, indices)
    out = field('balance', on)

    allocate(y(nv))
    y = this % source

    if (allocated(this % edge_terms)) then
       do k = 1, size(this % edge_terms)

          ! One edge term, computed for every edge at once. This is
          ! the only place the edge values are computed.
          call this % edge_terms(k) % apply(input_graph, input_data, edge_values)
          call edge_values % get_real_vector(z)

          ! And reduced onto the vertices through incidence, each
          ! edge touching its two
          ! ends exactly one time.
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

    ! A supplied buffer is overwritten, never added to.
    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine balance_apply

end module class_graph_balance
