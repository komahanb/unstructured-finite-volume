!=====================================================================!
! The balance - a vertex field operation composed from smaller ones.
!
! A balance is the net quantity of a cell: the sum of its face
! contributions, plus the cell's own source.
!
!                          z_e
!                   (i) ---------> (j)
!
!                   y_i = y_i - z_e
!                   y_j = y_j + z_e
!
! Every face's value is subtracted from the cell it leaves and added
! to the cell it enters. That is incidence, and it is the whole of the
! reduction through incidence.
!
! EXACTLY ONCE. Each face is visited one time and updates its two
! cells one time. Visit a face twice and the balance is wrong by that
! face; omit one and it is wrong by that face with the opposite sign.
! Neither produces a run-time error - it produces a solution that
! differs from the correct solution without any diagnostic, so the
! incidence count is verified in the test suite rather than assumed.
!
! A boundary face has no far cell, so its value is applied to the one
! cell it is incident to:
!
!                   (i) --------o
!
!                   y_i = y_i - z_b
!
!=====================================================================!
!
!                 WHY THIS IS NOT A NEW KIND OF OPERATION
!
! A balance is a vertex field operation like any other. It applies an
! edge operation in the process, but what it returns is a value per
! cell, which is exactly what its type declares.
!
! A solver calls the result a residual. That word names a stage in a
! solve, not an object in this library, so it does not appear here.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module operation_balance

  use util_precision  , only : dp
  use operation_action, only : operation, contract, binding
  use operation_action, only : emit
  use view_directed, only : directed_graph
  use field_calculus, only : field, FIELD_REAL
  use graph_fractal      , only : graph
  use field_stored  , only : stored_field
  use operation_differential, only : differential_operator

  implicit none

  private
  public :: balance

  !===================================================================!
  ! A balance stores its face terms and, if it has one, a number
  ! added to every cell.
  !
  ! The face terms are one concrete type storing a rule. A Fortran
  ! array stores a single dynamic type, so one concrete edge-derivative
  ! type is what lets the balance store its terms in a plain array.
  !===================================================================!

  type, extends(operation) :: balance

     type(differential_operator), allocatable :: edge_terms(:)

     real(dp) :: source = 0.0_dp

   contains

     procedure :: apply  => balance_apply

  end type balance

  interface balance
     module procedure create
  end interface balance

contains

  !===================================================================!
  ! Build a balance from what a model declares: the face terms, and
  ! one number per cell added as a source. Both optional - an empty
  ! balance returns zero.
  !===================================================================!

  type(balance) function create(edge_terms, source) result(this)

    type(differential_operator), intent(in), optional :: edge_terms(:)
    real(dp)  , intent(in), optional :: source

    if (present(edge_terms)) allocate(this % edge_terms, source=edge_terms)
    if (present(source))     this % source = source

    ! one argument: the state the balance is taken of, of any component count
    call this % declare_arguments(1, [contract(FIELD_REAL)], label='balance')

  end function create

  !===================================================================!
  ! Compute the balance.
  !
  !    1. start every cell at its own source term
  !    2. for each face term, evaluate every face
  !    3. reduce each edge onto the two vertices it touches, once,
  !       through incidence
  !===================================================================!

  subroutine balance_apply(this, input_graph, inputs, output)

    class(balance)    , intent(in)                 :: this
    class(directed_graph)      , intent(in)                 :: input_graph
    type(binding), intent(in), optional       :: inputs(:)
    class(field), allocatable, intent(inout) :: output

    class(field), allocatable :: edge_values

    type(stored_field)           :: out
    real(dp), allocatable :: y(:), z(:)
    integer               :: nv, ne, e, t, h, k

    nv = input_graph % num_vertices()
    ne = input_graph % num_edges()

    out = stored_field('balance', input_graph % vertex_set(), input_graph % num_vertices())

    allocate(y(nv))
    y = this % source

    if (allocated(this % edge_terms)) then
       do k = 1, size(this % edge_terms)

          ! One edge term, computed for every edge at once. This is
          ! the only place the edge values are computed.
          call this % edge_terms(k) % apply(input_graph, &
               & this % edge_terms(k) % bind(inputs), edge_values)
          call edge_values % real_vector(z)

          ! Then reduced onto the vertices through incidence, each
          ! edge updating its two ends exactly one time.
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

    call emit(out, output)

  end subroutine balance_apply

end module operation_balance
