!=====================================================================!
! THE FRACTAL KERNEL . ONE PRIMITIVE                        (spike)
!
! A graph is two branches; each branch is in exactly one epistemic
! status; a known branch realizes another graph:
!
!      graph
!         branch(2)         status  in  { NULL, UNKNOWN, KNOWN }
!                           KNOWN  ->  realization is a graph
!
! and NOTHING ELSE. Atoms are (NULL, NULL) nodes told apart by
! identity alone. Tuples, relations, operators, fields, meshes and
! the computational pair are INTERPRETATIONS of node shapes, owned
! by views above this module - the kernel holds generators only.
!
!                    TWO ABSENCES, ONE GROWTH LAW
!
! NULL is definite absence - the boundary face's missing far side.
! UNKNOWN is epistemic bottom - the seat not yet realized. They
! are different absences (COMPUTATIONAL-GRAPH.md), and here they
! are different STATUSES, never collapsed into unallocated storage.
! Knowledge only grows: UNKNOWN may become KNOWN, exactly once;
! NULL and KNOWN are final. Nine status pairs are nine VALUES of
! one type - a finite classification is never a hierarchy.
!
!                  IDENTITY IS A SEAT IN A UNIVERSE
!
! Storage is an arena: a universe signs the identity roll once
! (graph_identity), and every node it mints is identified by its
! seat. Two (NULL, NULL) atoms are two atoms because two seats.
! A handle is (universe token, seat) - literally: every handle
! carries the stamp of the signing that minted it, and a universe
! reassigned over is a DIFFERENT universe whose old citizens are
! refused, not misread. So a realization may sit in many branches
! (sharing), a branch may point back up the road it came by
! (cycles, tied with realize), and growth reallocates arrays
! without invalidating anyone, because a seat is an index, never
! an address.
!
! Three lifetime laws no guard can check for you: the universe
! must OUTLIVE its citizens' handles (awake guards seats and
! generations, never lifetimes - a dangled pointer is undefined
! before any guard runs); every arena declaration and every arena
! dummy on any road to node/citizen carries TARGET; and an arena
! never moves while it hosts - not an element of a growable
! collection, not a value passed by copy, not the object of an
! assignment.
!
!                     VALUES ARE COMPRESSED ATOMS
!
! An atom may carry one number. That is not a second ontology: it
! is the compressed representation of a numeral the structure
! could spell out and never should. Only atoms carry values, and
! meaning still arrives by interpretation, never by payload.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module fractal_graph

  use iso_fortran_env, only : dp => REAL64
  use graph_identity , only : token, mint_token

  implicit none

  private
  public :: BRANCH_NULL, BRANCH_UNKNOWN, BRANCH_KNOWN
  public :: graph_arena, graph, branch_spec
  public :: null_branch, unknown_branch, known_branch

  integer, parameter :: BRANCH_NULL    = 0
  integer, parameter :: BRANCH_UNKNOWN = 1
  integer, parameter :: BRANCH_KNOWN   = 2

  !===================================================================!
  ! The universe: flat storage for every citizen, one signature on
  ! the roll. Append-only, and SEAT ORDER IS MINT ORDER - the
  ! contiguity promise the storage laws above ride on. Structure is
  ! immutable; the only transition anywhere is UNKNOWN -> KNOWN
  ! inside these arrays.
  !===================================================================!

  type :: graph_arena

     type(token), private :: identity
     integer    , private :: n = 0

     integer , allocatable, private :: status(:,:)   ! (2, capacity)
     integer , allocatable, private :: target(:,:)   ! (2, capacity)
     logical , allocatable, private :: carries(:)
     real(dp), allocatable, private :: payload(:)

   contains

     procedure :: node       => mint_node
     procedure :: value_atom => mint_value_atom
     procedure :: citizen
     procedure :: population

  end type graph_arena

  interface graph_arena
     module procedure create_arena
  end interface graph_arena

  !===================================================================!
  ! The handle: a universe and a seat. Copying a handle is free and
  ! changes nothing; identity is the coordinate, never the copy.
  !===================================================================!

  type :: graph

     class(graph_arena), pointer, private :: universe => null()
     type(token)                , private :: stamp
     integer                    , private :: seat_at = 0

   contains

     procedure :: status  => status_of
     procedure :: known
     procedure :: realize
     procedure :: same_as
     procedure :: carries_value
     procedure :: value   => value_of
     procedure :: seat
     procedure :: universe_size

  end type graph

  !===================================================================!
  ! The branch as a TYPE was tried and evicted: it is status() plus
  ! a guarded known() - a composition, owned by call sites, per
  ! GENERATION. The ontology's 'branch(2)' is the pair of
  ! observations, not a bundle the kernel must sell.
  !===================================================================!

  type :: branch_spec

     integer    , private :: status = BRANCH_NULL
     type(graph), private :: realization

  end type branch_spec

contains

  type(graph_arena) function create_arena() result(this)

    this % identity = mint_token()
    allocate(this % status(2, 64), this % target(2, 64))
    allocate(this % carries(64), this % payload(64))
    this % status = BRANCH_NULL
    this % target = 0
    this % carries = .false.
    this % payload = 0.0_dp

  end function create_arena

  type(branch_spec) function null_branch() result(s)
    s % status = BRANCH_NULL
  end function null_branch

  type(branch_spec) function unknown_branch() result(s)
    s % status = BRANCH_UNKNOWN
  end function unknown_branch

  type(branch_spec) function known_branch(realization) result(s)

    type(graph), intent(in) :: realization

    if (.not. associated(realization % universe)) then
       error stop 'fractal: an unhatched graph realizes nothing'
    end if
    s % status = BRANCH_KNOWN
    s % realization = realization

  end function known_branch

  !===================================================================!
  ! The sole introduction form: one node, two branch specs, both at
  ! once - there is no first and second act to commute.
  !===================================================================!

  type(graph) function mint_node(this, left, right) result(g)

    class(graph_arena), target, intent(inout) :: this
    type(branch_spec)         , intent(in)    :: left, right

    if (.not. this % identity % declared()) then
       error stop 'fractal: a universe signs before it hosts'
    end if
    call admit(this, left)
    call admit(this, right)
    call ensure_room(this)

    this % n = this % n + 1
    this % status(1, this % n) = left  % status
    this % status(2, this % n) = right % status
    this % target(1, this % n) = left  % realization % seat_at
    this % target(2, this % n) = right % realization % seat_at

    g % universe => this
    g % stamp    =  this % identity
    g % seat_at  =  this % n

  end function mint_node

  type(graph) function mint_value_atom(this, v) result(g)

    class(graph_arena), target, intent(inout) :: this
    real(dp)                  , intent(in)    :: v

    g = mint_node(this, null_branch(), null_branch())
    this % carries(g % seat_at) = .true.
    this % payload(g % seat_at) = v

  end function mint_value_atom

  subroutine admit(this, s)

    class(graph_arena), target, intent(in) :: this
    type(branch_spec)         , intent(in) :: s

    if (s % status .eq. BRANCH_KNOWN) then
       if (.not. associated(s % realization % universe, this)) then
          error stop 'fractal: a branch realizes within its own universe'
       end if
    end if

  end subroutine admit

  subroutine ensure_room(this)

    class(graph_arena), intent(inout) :: this

    integer , allocatable :: si(:,:), ti(:,:)
    logical , allocatable :: ci(:)
    real(dp), allocatable :: pi(:)
    integer               :: cap

    cap = size(this % carries)
    if (this % n .lt. cap) return

    allocate(si(2, 2*cap), ti(2, 2*cap), ci(2*cap), pi(2*cap))
    si = BRANCH_NULL;  si(:, 1:cap) = this % status
    ti = 0;            ti(:, 1:cap) = this % target
    ci = .false.;      ci(1:cap)    = this % carries
    pi = 0.0_dp;       pi(1:cap)    = this % payload
    call move_alloc(si, this % status)
    call move_alloc(ti, this % target)
    call move_alloc(ci, this % carries)
    call move_alloc(pi, this % payload)

  end subroutine ensure_room

  !===================================================================!
  ! The enumeration face: the universe is the ground carrier of its
  ! citizens - citizen(k) and seat() are member(k) and local_index,
  ! the same contract the carriers have always answered.
  !===================================================================!

  type(graph) function citizen(this, k) result(g)

    class(graph_arena), target, intent(in) :: this
    integer                   , intent(in) :: k

    if (k .lt. 1 .or. k .gt. this % n) then
       error stop 'fractal: no citizen holds that seat'
    end if
    g % universe => this
    g % stamp    =  this % identity
    g % seat_at  =  k

  end function citizen

  pure integer function population(this)
    class(graph_arena), intent(in) :: this
    population = this % n
  end function population

  pure integer function seat(this)
    class(graph), intent(in) :: this
    seat = this % seat_at
  end function seat

  ! Not a generator: population's answer, re-exported across the
  ! capability wall - a handle may ask how large its world is
  ! without gaining the right to mint into it.
  integer function universe_size(this)
    class(graph), intent(in) :: this
    call awake(this)
    universe_size = this % universe % n
  end function universe_size

  !===================================================================!
  ! The observations: status, realization, identity, value.
  !===================================================================!

  integer function status_of(this, side)

    class(graph), intent(in) :: this
    integer     , intent(in) :: side

    call awake(this)
    call lawful(side)
    status_of = this % universe % status(side, this % seat_at)

  end function status_of

  type(graph) function known(this, side) result(r)

    class(graph), intent(in) :: this
    integer     , intent(in) :: side

    call awake(this)
    call lawful(side)
    select case (this % universe % status(side, this % seat_at))
    case (BRANCH_KNOWN)
       r % universe => this % universe
       r % stamp    =  this % stamp
       r % seat_at  =  this % universe % target(side, this % seat_at)
    case (BRANCH_UNKNOWN)
       error stop 'fractal: an unknown branch answers no realization - bottom is not a value'
    case default
       error stop 'fractal: a null branch answers no realization - absence is not a value'
    end select

  end function known

  !===================================================================!
  ! The sole transition: bottom becomes knowledge, once. Absence
  ! never does, and knowledge never re-arrives.
  !===================================================================!

  subroutine realize(this, side, realization)

    class(graph), intent(in) :: this
    integer     , intent(in) :: side
    type(graph) , intent(in) :: realization

    call awake(this)
    call lawful(side)
    if (.not. associated(realization % universe, this % universe)) then
       error stop 'fractal: a branch realizes within its own universe'
    end if
    select case (this % universe % status(side, this % seat_at))
    case (BRANCH_NULL)
       error stop 'fractal: absence is not ignorance - a null branch never realizes'
    case (BRANCH_KNOWN)
       error stop 'fractal: knowledge grows once - a known branch never realizes twice'
    case (BRANCH_UNKNOWN)
       this % universe % status(side, this % seat_at) = BRANCH_KNOWN
       this % universe % target(side, this % seat_at) = realization % seat_at
    end select

  end subroutine realize

  !===================================================================!
  ! Identity answers, never refuses: an unhatched handle equals
  ! nothing, itself included - the undeclared token's own law.
  !===================================================================!

  logical function same_as(this, other)

    class(graph), intent(in) :: this
    type(graph) , intent(in) :: other

    same_as = .false.
    if (.not. associated(this % universe)) return
    if (.not. associated(this % universe, other % universe)) return
    if (.not. this % stamp % matches(other % stamp)) return
    same_as = this % seat_at .eq. other % seat_at

  end function same_as

  logical function carries_value(this)

    class(graph), intent(in) :: this

    call awake(this)
    carries_value = this % universe % carries(this % seat_at)

  end function carries_value

  real(dp) function value_of(this)

    class(graph), intent(in) :: this

    call awake(this)
    if (.not. this % universe % carries(this % seat_at)) then
       error stop 'fractal: this graph carries no value'
    end if
    value_of = this % universe % payload(this % seat_at)

  end function value_of

  !===================================================================!
  ! The guards, and the canonical names for diagnostics.
  !===================================================================!

  subroutine awake(this)

    class(graph), intent(in) :: this

    if (.not. associated(this % universe)) then
       error stop 'fractal: an unhatched graph answers nothing'
    end if
    if (.not. this % universe % identity % matches(this % stamp)) then
       error stop 'fractal: this universe is not the one that minted you'
    end if
    if (this % seat_at .lt. 1 .or. this % seat_at .gt. this % universe % n) then
       error stop 'fractal: no citizen holds that seat'
    end if

  end subroutine awake

  subroutine lawful(side)

    integer, intent(in) :: side

    if (side .lt. 1 .or. side .gt. 2) then
       error stop 'fractal: a graph has exactly two branches'
    end if

  end subroutine lawful

end module fractal_graph
