!=====================================================================!
! GRAPH KERNEL
!
!     G = (B1, B2)
!     B in { NULL, UNKNOWN, KNOWN -> G }
!
! LAWS
!
!     graph -> branch -> graph, and no type stands between them
!     KNOWN implies associated(known)
!     NULL and UNKNOWN imply .not. associated(known)
!     NULL and UNKNOWN are distinct by status, not by association
!     (NULL, NULL) is a graph
!     graph identity is independent of branch state
!     branch references do not own their targets
!     identity is assigned once, and is not chosen
!
! The kernel carries structure and identity only. Numbers, symbols,
! indices and every other attribute are bound in graph_views.
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

     integer              :: status = GRAPH_NULL
     type(graph), pointer :: known  => null()

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
  ! Branch constructors. The only admitted constructions, and hence
  ! the point at which the association laws hold by construction.
  !
  ! known_branch is not pure: pointer assignment to an INTENT(IN)
  ! target is prohibited in a pure procedure (F2018 C1594).
  !===================================================================!

  pure type(graph_branch) function null_branch() result(this)

    this % status = GRAPH_NULL

  end function null_branch

  pure type(graph_branch) function unknown_branch() result(this)

    this % status = GRAPH_UNKNOWN

  end function unknown_branch

  type(graph_branch) function known_branch(that) result(this)

    type(graph), target, intent(in) :: that

    if (.not. that % identity % declared()) then
       error stop 'fractal_graph: KNOWN requires a graph with assigned identity'
    end if

    this % status = GRAPH_KNOWN
    this % known  => that

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
