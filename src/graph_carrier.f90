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
! The structural contract is four questions: how many (size), which
! one (member), all of them (members), and IS THIS ONE OF YOURS
! (has). Membership is a primitive, not a search: a relation
! signature will validate tuples against it, and must never have to
! enumerate a domain to ask whether one index belongs. The name a
! set was declared with is metadata for the reader - carried,
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
! Identity is an OPAQUE TOKEN. Its parts are private, so no caller
! can compose one with chosen contents; the only ways a token moves
! are minting - fresh, unrepeatable - and whole-object copy, which
! is precisely what "the same declared domain" means. A set signs
! once, at declaration, and a second signing is refused loudly: a
! domain that changed identity mid-life would let one question
! answer two ways. The token today is an (image, serial) pair, so
! two coarray images can never mint the same stamp; the
! representation stays free to grow because nothing outside this
! module can read its parts.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_carrier

  implicit none

  private
  public :: member_set, counted_set, token, mint_token

  !===================================================================!
  ! The stamp roll of this image. Serial zero is reserved for the
  ! undeclared: a default-initialized token is no identity at all,
  ! and matches nothing, itself included.
  !===================================================================!

  integer, save :: last_serial = 0

  !===================================================================!
  ! The opaque token. Parts private; minting and copying are the
  ! only ways one comes to exist. serial() is a read-only diagnostic
  ! for messages and tests - matches is the law.
  !===================================================================!

  type :: token

     integer, private :: image  = 0
     integer, private :: serial = 0

   contains

     procedure :: matches
     procedure :: declared
     procedure :: serial_number

  end type token

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

     procedure(set_size_interface)   , deferred :: size
     procedure(set_member_interface) , deferred :: member
     procedure(set_members_interface), deferred :: members
     procedure(set_has_interface)    , deferred :: has

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

     procedure :: size    => counted_size
     procedure :: member  => counted_member
     procedure :: members => counted_members
     procedure :: has     => counted_has

  end type counted_set

  interface counted_set
     module procedure create_counted
  end interface counted_set

contains

  !===================================================================!
  ! mint_token hands out the next stamp of this image: fresh,
  ! unrepeatable, contents unchoosable. Higher levels with their own
  ! identities (relations, graphs) mint here too, so the whole tower
  ! draws on one roll per image.
  !===================================================================!

  type(token) function mint_token()

    last_serial          = last_serial + 1
    mint_token % serial  = last_serial
    mint_token % image   = this_image()

  end function mint_token

  !===================================================================!
  ! The token's three answers.
  !===================================================================!

  pure logical function matches(this, other)

    class(token), intent(in) :: this
    type(token) , intent(in) :: other

    matches = (this % serial /= 0)              .and. &
         &    (this % serial == other % serial) .and. &
         &    (this % image  == other % image)

  end function matches

  pure logical function declared(this)

    class(token), intent(in) :: this

    declared = this % serial /= 0

  end function declared

  pure integer function serial_number(this)

    class(token), intent(in) :: this

    serial_number = this % serial

  end function serial_number

  !===================================================================!
  ! declare stamps a set with a fresh identity and, if given, its
  ! name. A set signs once; a second signing stops the program,
  ! because a silent refusal would leave the caller believing a
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
  ! id reads the stamp's serial - a diagnostic, local to this image;
  ! zero means never declared. same_as is the law.
  !===================================================================!

  pure integer function id(this)

    class(member_set), intent(in) :: this

    id = this % identity % serial_number()

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

end module graph_carrier
