!=====================================================================!
! THE SHIFTED-LAPLACIAN FIXTURE - test-local, deliberately outside
! src/: one chain has not earned a production shifted operator. It
! represents
!
!      A(q) = 2 q - L(q)
!
! where L is the PRODUCTION vertex Laplacian, and it is the whole
! point of this tower that the second term is not written here.
!
! This adapter OWNS NO GRAPH. Its graph arrives through domain()
! and apply(), and it hands that graph straight to the production
! differential operator. It does not inspect edge_tail, does not
! inspect edge_head, implements no incidence and reproduces no
! Laplacian loop: topology belongs entirely to the operator it
! delegates to.
!
! The call chain IS the experiment:
!
!      caller / minimizer
!          |  supplies a graph
!          v
!      shifted_laplacian          (this file - no topology)
!          |
!          v
!      production laplacian       (traverses the graph)
!          |
!          v
!      input_graph incidence
!
! Because the adapter stores no graph, applying it on a different
! host yields a different mathematical action - which is how this
! tower proves behaviourally that the graph a minimizer carries is
! load-bearing rather than scenery.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module shifted_laplacian_fixture

  use iso_fortran_env  , only : dp => REAL64
  use graph_carrier    , only : member_set
  use graph_grammar    , only : graph, graph_field, graph_operation
  use class_graph_field, only : field
  use class_graph_differential_operator, only : differential_operator, &
       &                                        laplacian

  implicit none

  private
  public :: shifted_laplacian

  type, extends(graph_operation) :: shifted_laplacian
   contains
     procedure :: name   => shifted_name
     procedure :: domain => shifted_domain
     procedure :: apply  => shifted_apply
  end type shifted_laplacian

contains

  pure function shifted_name(this) result(name)
    class(shifted_laplacian), intent(in) :: this
    character(len=:), allocatable :: name
    name = 'shifted laplacian'
  end function shifted_name

  !===================================================================!
  ! The operation answers on whatever graph it is handed - it holds
  ! no domain of its own.
  !===================================================================!

  subroutine shifted_domain(this, input_graph, domain)

    class(shifted_laplacian), intent(in) :: this
    class(graph), intent(in) :: input_graph
    class(member_set), allocatable, intent(out) :: domain

    allocate(domain, source=input_graph % vertex_set())

  end subroutine shifted_domain

  !===================================================================!
  ! A(q) = 2q - L(q). The state is checked against THIS graph's
  ! vertex carrier by identity - a field of the right size on a
  ! foreign carrier is refused - and the Laplacian is the production
  ! one, applied to the same graph.
  !===================================================================!

  subroutine shifted_apply(this, input_graph, input_data, output)

    class(shifted_laplacian), intent(in)           :: this
    class(graph), intent(in)                       :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(differential_operator)     :: lap
    type(field)                     :: out
    class(graph_field), allocatable :: lq
    class(member_set), allocatable  :: dom
    real(dp), allocatable           :: q(:), l(:)

    if (.not. present(input_data)) then
       error stop 'shifted laplacian: the action needs a state to read'
    end if
    if (size(input_data) /= 1) then
       error stop 'shifted laplacian: the action reads exactly one state'
    end if

    call input_data(1) % domain(dom)
    if (.not. dom % same_as(input_graph % vertex_set())) then
       error stop 'shifted laplacian: the state must live on this graph''s vertex carrier'
    end if

    ! The topology is consumed HERE, by production, on THIS graph.
    lap = laplacian(coefficient=1.0_dp, spacing=1.0_dp, measure=1.0_dp)
    call lap % apply(input_graph, input_data, lq)

    call input_data(1) % get_real_vector(q)
    call lq % get_real_vector(l)

    out = field('shifted laplacian', input_graph % vertex_set())
    call out % set_real_vector(2.0_dp * q - l)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine shifted_apply

end module shifted_laplacian_fixture
