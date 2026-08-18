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
     procedure :: similar_as

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

  !===================================================================!
  ! STRUCTURAL COMPARISON: similar_as
  !
  ! Two graphs are similar if their branch states match structurally,
  ! regardless of identity. For KNOWN branches, similarity is tested
  ! recursively. Distinguishes from same_as (identity) to enable
  ! structural change detection in evolving systems.
  !
  ! SAME_AS vs SIMILAR_AS:
  !
  !   g1.same_as(g2)     => g1 and g2 are the identical graph
  !   g1.similar_as(g2)  => g1 and g2 have the same structure
  !
  ! Use cases:
  !
  !   1. CONVERGENCE DETECTION (adaptive refinement):
  !      do
  !        refined = refiner.apply(mesh)
  !        if (refined.similar_as(mesh)) exit  ! No finer structure
  !        mesh = refined
  !      end do
  !
  !   2. COARSENING TERMINATION (multigrid levels):
  !      do level = 1, max_levels
  !        coarse = coarsener.apply(current)
  !        if (coarse.similar_as(current)) exit  ! Already coarsest
  !        current = coarse
  !      end do
  !
  !   3. STRUCTURE VALIDATION (post-transformation):
  !      partition = partitioner.apply(mesh)
  !      call assert(partition.similar_as(expected_structure), &
  !           & "partition structure mismatch")
  !
  !   4. EQUILIBRIUM DETECTION (evolving systems):
  !      state_old = state
  !      call marcher.step(state)
  !      if (state.similar_as(state_old)) converged = .true.
  !
  ! STRUCTURE COMPARISON LAW:
  !
  !   g1 and g2 are similar iff:
  !   - branch(1) states match (NULL, UNKNOWN, or KNOWN)
  !   - branch(2) states match
  !   - for each KNOWN branch, the referenced graphs are recursively similar
  !
  ! Example (sequences [a,b] vs [c,d]):
  !
  !   g1: (a, (b, NULL))      g2: (c, (d, NULL))
  !   both KNOWN -> KNOWN     both KNOWN -> KNOWN
  !   -> similar = .true.     (structure same, identity different)
  !
  !===================================================================!

  recursive logical function similar_as(this, other)

    class(graph), intent(in) :: this
    type(graph) , intent(in) :: other

    type(graph), pointer :: this_known_1, this_known_2
    type(graph), pointer :: other_known_1, other_known_2

    ! Branch(1) states must match
    if (this % branch(1) % status() /= other % branch(1) % status()) then
       similar_as = .false.
       return
    end if

    ! Branch(2) states must match
    if (this % branch(2) % status() /= other % branch(2) % status()) then
       similar_as = .false.
       return
    end if

    ! If branch(1) is KNOWN, recursively check the referenced graphs
    if (this % branch(1) % status() == GRAPH_KNOWN) then
       this_known_1  => this % branch(1) % known()
       other_known_1 => other % branch(1) % known()
       if (.not. this_known_1 % similar_as(other_known_1)) then
          similar_as = .false.
          return
       end if
    end if

    ! If branch(2) is KNOWN, recursively check the referenced graphs
    if (this % branch(2) % status() == GRAPH_KNOWN) then
       this_known_2  => this % branch(2) % known()
       other_known_2 => other % branch(2) % known()
       if (.not. this_known_2 % similar_as(other_known_2)) then
          similar_as = .false.
          return
       end if
    end if

    similar_as = .true.

  end function similar_as

end module fractal_graph
