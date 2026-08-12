!=====================================================================!
! LEVEL 0 OF THE NEW TOWER . THE CARRIERS
!
! The ground object of the relation-centered tower (AGENTS.md): a
! MEMBER SET,
!
!      A = { a_1, ..., a_n }
!
! It answers what members may exist, and nothing else. No structure,
! no relations, no graph words - vertices, edges, cells and faces
! are what higher levels call the sets they declare here.
!
! The structural contract is five questions: how many (size), which
! one (member), all of them (members), IS THIS ONE OF YOURS (has),
! and WHERE DOES IT STAND (local_index - the inverse of member,
! zero for an outsider). Membership and position are primitives,
! not searches: a relation signature validates tuples against has,
! and an indexed relation finds its row through local_index, so
! neither ever has to enumerate a domain. The two enumeration laws
! bind them:
!
!      member(local_index(a)) = a        for every a in A
!      local_index(member(i)) = i        for i = 1 .. size
!
! which forces enumeration to be injective: a member set holds each
! member once - set semantics on the carriers themselves. The name
! a set was declared with is metadata for the reader - carried,
! printable, and no part of the mathematics.
!
!                       STRUCTURAL IDENTITY
!
! Two member sets are the same only when they are the SAME DECLARED
! DOMAIN - never because they hold equal integers, or count the same,
! or happen to enumerate identical values,
!
!      cells  = { 1 2 3 4 }        four cells
!      faces  = { 1 2 3 4 }        four faces - a different world
!
! The identity machinery - the opaque token, its minting, its one
! law - lives beneath the tower in graph_identity, because relations
! sign the same way and the relational graph will be the third to
! do so. Here the law is only applied: a set signs once, at
! declaration; a second signing is refused loudly; id() answers the
! whole token, so identity is never mistaken for an image-local
! integer.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_carrier

  use graph_identity, only : token, mint_token

  implicit none

  private
  public :: member_set, counted_set

  !===================================================================!
  ! The abstract member set: an identity, a count, its members, and
  ! membership. The contract is deliberately this small (AGENTS.md
  ! 4.2); subsets, tags, parts and every other refinement arrive on
  ! higher levels as predicates and relations, never here.
  !===================================================================!

  type, abstract :: member_set

     type(token)                  , private :: identity
     character(len=:), allocatable, private :: label

   contains

     !----------------------------------------------------------------!
     ! The structural questions, deferred to each concretion.
     !----------------------------------------------------------------!

     procedure(set_size_interface)       , deferred :: size
     procedure(set_member_interface)     , deferred :: member
     procedure(set_members_interface)    , deferred :: members
     procedure(set_has_interface)        , deferred :: has
     procedure(set_local_index_interface), deferred :: local_index

     !----------------------------------------------------------------!
     ! Identity, answered once for every concretion.
     !----------------------------------------------------------------!

     procedure :: declare
     procedure :: id
     procedure :: same_as

     !----------------------------------------------------------------!
     ! Metadata, not mathematics: the declared name, or ''.
     !----------------------------------------------------------------!

     procedure :: name

  end type member_set

  abstract interface

     pure integer function set_size_interface(this)
       import member_set
       class(member_set), intent(in) :: this
     end function set_size_interface

     pure integer function set_member_interface(this, local_index)
       import member_set
       class(member_set), intent(in) :: this
       integer          , intent(in) :: local_index
     end function set_member_interface

     pure subroutine set_members_interface(this, indices)
       import member_set
       class(member_set)   , intent(in)  :: this
       integer, allocatable, intent(out) :: indices(:)
     end subroutine set_members_interface

     pure logical function set_has_interface(this, member)
       import member_set
       class(member_set), intent(in) :: this
       integer          , intent(in) :: member
     end function set_has_interface

     pure integer function set_local_index_interface(this, member)
       import member_set
       class(member_set), intent(in) :: this
       integer          , intent(in) :: member
     end function set_local_index_interface

  end interface

  !===================================================================!
  ! The first concrete citizen: a counted set, members 1 to n. Every
  ! domain the repository names today - cells, faces, points, parts,
  ! instants - enumerates exactly so, which makes one integer its
  ! whole storage and membership one comparison. A domain that must
  ! list its members arrives as a second concretion the day
  ! something needs it.
  !===================================================================!

  type, extends(member_set) :: counted_set

     integer, private :: n = 0

   contains

     procedure :: size        => counted_size
     procedure :: member      => counted_member
     procedure :: members     => counted_members
     procedure :: has         => counted_has
     procedure :: local_index => counted_local_index

  end type counted_set

  interface counted_set
     module procedure create_counted
  end interface counted_set

contains

  !===================================================================!
  ! declare stamps a set with a fresh identity and, if given, its
  ! name. A set signs once; a second signing stops the program,
  ! because a silent second stamp would leave the caller believing a
  ! domain it never made.
  !===================================================================!

  subroutine declare(this, name)

    class(member_set), intent(inout)        :: this
    character(len=*) , intent(in), optional :: name

    if (this % identity % declared()) then
       error stop 'graph_carrier: a domain never signs twice'
    end if

    this % identity = mint_token()
    if (present(name)) this % label = name

  end subroutine declare

  !===================================================================!
  ! id answers the whole opaque token - the identity itself, honest
  ! across images, never a bare local integer.
  !===================================================================!

  pure type(token) function id(this)

    class(member_set), intent(in) :: this

    id = this % identity

  end function id

  !===================================================================!
  ! same_as: one declared domain, or not. Equal contents prove
  ! nothing; an undeclared set equals nothing, itself included.
  !===================================================================!

  pure logical function same_as(this, other)

    class(member_set), intent(in) :: this
    class(member_set), intent(in) :: other

    same_as = this % identity % matches(other % identity)

  end function same_as

  !===================================================================!
  ! name answers the declared label, or '' for a set declared
  ! nameless. Metadata only: no law reads it.
  !===================================================================!

  function name(this)

    class(member_set), intent(in) :: this
    character(len=:), allocatable :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function name

  !===================================================================!
  ! Declare a counted domain: a name for the reader, a count for the
  ! mathematics, a fresh stamp for the identity.
  !===================================================================!

  type(counted_set) function create_counted(name, n) result(this)

    character(len=*), intent(in) :: name
    integer         , intent(in) :: n

    this % n = max(0, n)
    call this % declare(name)

  end function create_counted

  pure integer function counted_size(this)

    class(counted_set), intent(in) :: this

    counted_size = this % n

  end function counted_size

  pure integer function counted_member(this, local_index)

    class(counted_set), intent(in) :: this
    integer           , intent(in) :: local_index

    counted_member = local_index

  end function counted_member

  pure subroutine counted_members(this, indices)

    class(counted_set)  , intent(in)  :: this
    integer, allocatable, intent(out) :: indices(:)

    integer :: k

    allocate(indices(this % n))
    do k = 1, this % n
       indices(k) = k
    end do

  end subroutine counted_members

  !===================================================================!
  ! Membership in one comparison - the primitive a relation
  ! signature leans on, never an enumeration and a search.
  !===================================================================!

  pure logical function counted_has(this, member)

    class(counted_set), intent(in) :: this
    integer           , intent(in) :: member

    counted_has = (member >= 1) .and. (member <= this % n)

  end function counted_has

  !===================================================================!
  ! Where a member stands - for the counted set, where it says. Zero
  ! for an outsider, so an indexed relation reads absence without a
  ! second question.
  !===================================================================!

  pure integer function counted_local_index(this, member)

    class(counted_set), intent(in) :: this
    integer           , intent(in) :: member

    if (member >= 1 .and. member <= this % n) then
       counted_local_index = member
    else
       counted_local_index = 0
    end if

  end function counted_local_index

end module graph_carrier
