!=====================================================================!
! FRACTAL GRAPH
!
!     G = (B1, B2)
!     B in { NULL, UNKNOWN, KNOWN -> G }
!
! LAWS
!
!     graph -> branch -> graph, and no type stands between them
!     status == KNOWN iff associated(known)
!     NULL and UNKNOWN imply .not. associated(known)
!     NULL and UNKNOWN are distinct by status, not by association
!     (NULL, NULL) is a graph
!     graph identity is independent of branch state
!     branch references do not own their targets
!     identity is assigned once, and is not chosen
!
! The status and the reference are private, so no caller can set one
! without the other. A branch enters existence only through the three
! constructors below, each of which establishes the iff.
!
! branch(2) stays a public component: assigning a whole branch value
! cannot break the iff, and the recursion stays visible in the
! declarations and in every navigation.
!
! The kernel carries shape, status, reference and identity. Numbers,
! symbols and indices are bound in graph_views.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module fractal_graph

  use graph_identity, only : token, mint_token

  implicit none

  private
  public :: graph, graph_branch
  public :: GRAPH_NULL, GRAPH_UNKNOWN, GRAPH_KNOWN
  public :: null_branch, unknown_branch, known_branch

  ! The three branch states are values, not subtypes.

  integer, parameter :: GRAPH_NULL    = 0
  integer, parameter :: GRAPH_UNKNOWN = 1
  integer, parameter :: GRAPH_KNOWN   = 2

  ! graph_branch precedes graph: a forward reference to an undefined
  ! derived type is admitted only as a pointer component.

  type :: graph_branch

     private
     integer              :: status_ = GRAPH_NULL
     type(graph), pointer :: known_  => null()

   contains

     procedure :: status
     procedure :: known

  end type graph_branch

  type :: graph

     type(graph_branch)   :: branch(2)
     type(token), private :: identity

   contains

     procedure :: declare
     procedure :: id
     procedure :: same_as

  end type graph

contains

  !===================================================================!
  ! Branch queries. known answers a disassociated pointer whenever
  ! the status is NULL or UNKNOWN, so the iff is directly observable.
  !
  ! known is not pure: a pure function result may not be associated
  ! with a target reached through an INTENT(IN) dummy (F2018 C1594).
  !===================================================================!

  pure integer function status(this)

    class(graph_branch), intent(in) :: this

    status = this % status_

  end function status

  function known(this) result(g)

    class(graph_branch), intent(in) :: this
    type(graph), pointer            :: g

    g => this % known_

  end function known

  !===================================================================!
  ! Branch constructors. The only admitted introductions of a branch
  ! value, and hence the only place the iff can be established. The
  ! structure constructor is unavailable outside this module because
  ! both components are private.
  !===================================================================!

  pure type(graph_branch) function null_branch() result(this)

    this % status_ = GRAPH_NULL

  end function null_branch

  pure type(graph_branch) function unknown_branch() result(this)

    this % status_ = GRAPH_UNKNOWN

  end function unknown_branch

  type(graph_branch) function known_branch(that) result(this)

    type(graph), target, intent(in) :: that

    if (.not. that % identity % declared()) then
       error stop 'fractal_graph: KNOWN requires a graph with assigned identity'
    end if

    this % status_ = GRAPH_KNOWN
    this % known_  => that

  end function known_branch

  !===================================================================!
  ! Identity. Minted, never chosen; assigned once; the sole equality.
  !===================================================================!

  subroutine declare(this)

    class(graph), intent(inout) :: this

    if (this % identity % declared()) then
       error stop 'fractal_graph: graph identity is assigned once'
    end if

    this % identity = mint_token()

  end subroutine declare

  pure type(token) function id(this)

    class(graph), intent(in) :: this

    id = this % identity

  end function id

  pure logical function same_as(this, other)

    class(graph), intent(in) :: this
    type(graph) , intent(in) :: other

    same_as = this % identity % matches(other % identity)

  end function same_as

end module fractal_graph
