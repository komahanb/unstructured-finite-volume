!=====================================================================!
! THE TRIANGULAR DECAY FIXTURE - earned at LEVEL 6.
!
! The tower's dynamical action, and the first operation it has ever
! built:
!
!      S : Q -> Q          S([x, y]) = [ x, y - x ]
!
! which under the house convention qdot = -S(q) is the coupled
! decay
!
!      dx/dt = -x
!      dy/dt =  x - y
!
! The name is the matrix: S is lower triangular, so x decays alone
! and y is driven by it.
!
!                    WHAT THIS OPERATION STORES
!
! IT STORES Q, AND IT STORES NO GRAPH.
!
! That single sentence is the Level-6 experiment. The
! graph_operation contract hands every operation an input_graph,
! and the established habit - visible in the marching fixtures, in
! operation_fitting, in reduction and broadcast - is to answer
!
!      call input_graph % all_vertices(domain)
!
! and let the host decide what the mathematics is about. This
! operation does not. Its domain is a thing it was constructed
! with, it answers that thing, and it validates its input against
! it. The input_graph is accepted, unread, and ignored.
!
! THAT IS NOT A CRITICISM OF THE CONTRACT. The partitioned tower
! settled seam A1 on production evidence: the host is a real
! conduit for actions that consume topology, and cannot be removed
! generically. This action simply is not one of those. For a
! triangular 2x2 decay there is no topology to traverse, so the
! host is context and nothing more - which is precisely the
! configuration that makes it possible to ask whether anything
! DOWNSTREAM quietly reads the host when it should be reading the
! action.
!
! Nothing here is derivative machinery. There is no tangent, no
! adjoint, no transpose and no linearization: S is the ODE action
! itself, written once, in the plainest arithmetic available.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module triangular_decay_fixture

  use iso_fortran_env  , only : dp => REAL64
  use operation_action, only : graph_operation
  use view_directed, only : directed_graph
  use field_calculus, only : graph_field
  use graph_fractal    , only : set_graph => graph
  use field_stored, only : field

  implicit none

  private
  public :: triangular_decay

  !===================================================================!
  ! The action carries its own domain. One component, and it is the
  ! state carrier Q - never a graph, and never a vertex set.
  !===================================================================!

  type, extends(graph_operation) :: triangular_decay

     ! WHICH domain, and how many coordinates stand in it. Both of
     ! the action's questions are answered from these two, so no map
     ! enters this type - an operation that asks nothing about
     ! membership carries nothing that answers it.
     type(set_graph) :: state
     integer         :: n_state = 0

   contains

     procedure :: name   => decay_name
     procedure :: domain => decay_domain
     procedure :: apply  => decay_apply

  end type triangular_decay

  interface triangular_decay
     module procedure create_decay
  end interface triangular_decay

contains

  !===================================================================!
  ! Built FROM the state domain, and from nothing else.
  !===================================================================!

  type(triangular_decay) function create_decay(state, n_state) result(this)

    type(set_graph), intent(in) :: state
    integer        , intent(in) :: n_state

    if (.not. state % same_as(state)) then
       error stop 'triangular_decay: an action needs a declared state domain'
    end if

    this % state   = state
    this % n_state = n_state

  end function create_decay

  pure function decay_name(this) result(name)

    class(triangular_decay), intent(in) :: this
    character(len=:), allocatable       :: name

    associate (u1 => this); end associate

    name = 'triangular decay'

  end function decay_name

  !===================================================================!
  ! THE answer that makes this fixture the experiment: the domain is
  ! the STORED state carrier. input_graph is present because the
  ! contract says so, and is not consulted.
  !===================================================================!

  subroutine decay_domain(this, input_graph, domain, nentries)

    class(triangular_decay), intent(in)  :: this
    class(directed_graph)           , intent(in)  :: input_graph
    type(set_graph)        , intent(out) :: domain
    integer                , intent(out) :: nentries

    associate (u1 => input_graph); end associate

    domain   = this % state
    nentries = this % n_state

  end subroutine decay_domain

  !===================================================================!
  ! S([x, y]) = [x, y - x], with the input checked against the
  ! stored domain by IDENTITY - not by size. A field of the right
  ! width on the wrong carrier is a different mathematical object,
  ! and this operation says so rather than computing something
  ! plausible.
  !===================================================================!

  subroutine decay_apply(this, input_graph, input_data, output)

    class(triangular_decay), intent(in)            :: this
    class(directed_graph)           , intent(in)            :: input_graph
    class(graph_field), intent(in), optional       :: input_data(:)
    class(graph_field), allocatable, intent(inout) :: output

    type(set_graph)       :: given
    type(field)           :: out
    real(dp), allocatable          :: q(:), s(:)

    associate (u1 => input_graph); end associate

    if (.not. present(input_data)) then
       error stop 'triangular_decay: the action needs a state to act on'
    end if

    given = input_data(1) % domain()
    if (.not. given % same_as(this % state)) then
       error stop 'triangular_decay: the state must live on the action''s own domain'
    end if

    call input_data(1) % get_real_vector(q)
    if (size(q) /= this % n_state) then
       error stop 'triangular_decay: one number per state coordinate'
    end if

    allocate(s(size(q)))
    s(1) = q(1)
    s(2) = q(2) - q(1)

    out = field('triangular decay', this % state, this % n_state, ncomp=1)
    call out % set_real_vector(s)

    if (allocated(output)) deallocate(output)
    allocate(output, source=out)

  end subroutine decay_apply

end module triangular_decay_fixture
