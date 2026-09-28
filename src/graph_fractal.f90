!=====================================================================!
! FRACTAL GRAPH
!
! A graph as a recursive binary composition:
!
!     G = (B1, B2)
!     B in { NULL, UNKNOWN, KNOWN -> G }
!
! The two branches are symmetric in structure, independent in semantics.
! They may represent complementary views of the same domain (structure
! and parameters, skeleton and learning, model and correction, or any
! other dual decomposition a caller defines). A caller interprets what
! each branch means; the graph enforces only that both exist and are
! recursively well-formed.
!
! LAWS
!
!     graph -> branch -> graph, and no type is defined between them
!     status == KNOWN iff associated(known)
!     NULL and UNKNOWN imply .not. associated(known)
!     NULL and UNKNOWN are distinct by status, not by association
!     (NULL, NULL) is a graph
!     graph identity is independent of branch state
!     each branch may store independent data or interpretation
!     branch references do not own their targets
!     identity is assigned once, and is not chosen
!
! The status and the reference are private, so no caller can set one
! without the other. A branch is created only through the three
! constructors below, each of which establishes the iff.
!
! branch(2) stays a public component: assigning a whole branch value
! cannot break the iff, and the recursion stays visible in the
! declarations and in every navigation. The binary structure is exposed
! so callers can exploit the natural pairing of the two branches for
! their own semantics.
!
! The kernel stores shape, status, reference and identity. Numbers,
! symbols and indices are bound in graph_views. Interpretation belongs
! to the caller: this module defines only the scaffold.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_fractal

  use token_identity, only : token, next_token

  implicit none

  private
  public :: graph, branch
  public :: BRANCH_NULL, BRANCH_UNKNOWN, BRANCH_KNOWN
  public :: null_branch, unknown_branch, known_branch

  ! The three branch states are values, not subtypes.

  integer, parameter :: BRANCH_NULL    = 0
  integer, parameter :: BRANCH_UNKNOWN = 1
  integer, parameter :: BRANCH_KNOWN   = 2

  ! branch precedes graph: a forward reference to an undefined
  ! derived type is admitted only as a pointer component.

  type :: branch

     private
     integer              :: status_ = BRANCH_NULL
     type(graph), pointer :: known_  => null()

   contains

     procedure :: status
     procedure :: known

  end type branch

  type :: graph

     ! Two branches, independent in semantics but symmetric in structure.
     ! A caller may assign different roles to branch(1) and branch(2),
     ! such as representation and correction, constraint and freedom,
     ! or any complementary decomposition of the domain. The graph
     ! enforces only well-formedness; interpretation is the caller's.
     type(branch)   :: branch(2)
     type(token), private :: identity

   contains

     procedure :: declare
     procedure :: id
     procedure :: same_as

  end type graph

contains

  !===================================================================!
  ! Branch queries. known returns a disassociated pointer whenever
  ! the status is NULL or UNKNOWN, so the iff is directly observable.
  !
  ! known is not pure: a pure function result may not be associated
  ! with a target reached through an INTENT(IN) dummy (F2018 C1594).
  !===================================================================!

  pure integer function status(this)

    class(branch), intent(in) :: this

    status = this % status_

  end function status

  function known(this) result(g)

    class(branch), intent(in) :: this
    type(graph), pointer            :: g

    g => this % known_

  end function known

  !===================================================================!
  ! Branch constructors. The only admitted introductions of a branch
  ! value, and hence the only place the iff can be established. The
  ! structure constructor is unavailable outside this module because
  ! both components are private.
  !===================================================================!

  pure type(branch) function null_branch() result(this)

    this % status_ = BRANCH_NULL

  end function null_branch

  pure type(branch) function unknown_branch() result(this)

    this % status_ = BRANCH_UNKNOWN

  end function unknown_branch

  type(branch) function known_branch(that) result(this)

    type(graph), target, intent(in) :: that

    if (.not. that % identity % declared()) then
       error stop 'graph_fractal: known_branch requires that to have an assigned identity, &
            &but that % identity % declared() is false'
    end if

    this % status_ = BRANCH_KNOWN
    this % known_  => that

  end function known_branch

  !===================================================================!
  ! Identity. Generated, never chosen; assigned once; the sole equality.
  !===================================================================!

  subroutine declare(this)

    class(graph), intent(inout) :: this

    if (this % identity % declared()) then
       error stop 'graph_fractal: declare was called on a graph whose identity is already &
            &assigned; identity is assigned once'
    end if

    this % identity = next_token()

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

end module graph_fractal
