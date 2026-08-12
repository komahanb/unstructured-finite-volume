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
!                       STRUCTURAL IDENTITY
!
! Two member sets are the same only when they are the SAME DECLARED
! DOMAIN - never because they hold equal integers, or count the same,
! or happen to enumerate identical values,
!
!      cells  = { 1 2 3 4 }        four cells
!      faces  = { 1 2 3 4 }        four faces - a different world
!
! So every construction stamps a fresh identity, copies carry the
! stamp along, and same_as reads the stamp alone. A relation
! signature will refer to these identities, not to array contents.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module graph_carrier

  implicit none

  private
  public :: member_set, counted_set, new_identity

  !===================================================================!
  ! The stamp roll. Every declared domain takes the next number, so
  ! no two declarations ever collide. Zero is reserved for the
  ! undeclared - a default-initialized set that was never constructed
  ! is no domain at all, and same_as says so.
  !===================================================================!

  integer, save :: last_identity = 0

  !===================================================================!
  ! The abstract member set: identity, a name for the reader, a
  ! count, and its members. The contract is deliberately this small
  ! (AGENTS.md 4.2); subsets, tags, parts and every other refinement
  ! arrive on higher levels as predicates and relations, never here.
  !===================================================================!

  type, abstract :: member_set

     integer, private :: identity = 0

   contains

     procedure(set_name_interface)   , deferred :: name
     procedure(set_size_interface)   , deferred :: size
     procedure(set_member_interface) , deferred :: member
     procedure(set_members_interface), deferred :: members

     procedure :: id
     procedure :: same_as
     procedure :: stamp

  end type member_set

  abstract interface

     function set_name_interface(this) result(name)
       import member_set
       class(member_set), intent(in) :: this
       character(len=:), allocatable :: name
     end function set_name_interface

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

  end interface

  !===================================================================!
  ! The first concrete citizen: a counted set, members 1 to n. Every
  ! domain the repository names today - cells, faces, points, parts,
  ! instants - enumerates exactly so, which makes one integer its
  ! whole storage. A domain that must list its members arrives as a
  ! second concretion the day something needs it.
  !===================================================================!

  type, extends(member_set) :: counted_set

     character(len=:), allocatable, private :: label
     integer                      , private :: n = 0

   contains

     procedure :: name    => counted_name
     procedure :: size    => counted_size
     procedure :: member  => counted_member
     procedure :: members => counted_members

  end type counted_set

  interface counted_set
     module procedure create_counted
  end interface counted_set

contains

  !===================================================================!
  ! new_identity hands the next stamp to whoever declares a domain -
  ! this module's own constructors, and any concretion that lives
  ! elsewhere.
  !===================================================================!

  integer function new_identity()

    last_identity = last_identity + 1
    new_identity  = last_identity

  end function new_identity

  !===================================================================!
  ! id reads the stamp; zero means never declared.
  !===================================================================!

  pure integer function id(this)

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

    same_as = (this % identity /= 0) .and. &
         &    (this % identity == other % identity)

  end function same_as

  !===================================================================!
  ! stamp marks a set with an identity already minted. Constructors
  ! outside this module use it in one breath with new_identity; the
  ! stamp is refused if the set already carries one, because a domain
  ! that changes identity mid-life would let one question answer two
  ! ways.
  !===================================================================!

  subroutine stamp(this, identity)

    class(member_set), intent(inout) :: this
    integer          , intent(in)    :: identity

    if (this % identity == 0) this % identity = identity

  end subroutine stamp

  !===================================================================!
  ! Declare a counted domain: a name for the reader, a count for the
  ! mathematics, a fresh stamp for the identity.
  !===================================================================!

  type(counted_set) function create_counted(name, n) result(this)

    character(len=*), intent(in) :: name
    integer         , intent(in) :: n

    this % label = name
    this % n     = max(0, n)
    call this % stamp(new_identity())

  end function create_counted

  function counted_name(this) result(name)

    class(counted_set), intent(in) :: this
    character(len=:), allocatable  :: name

    if (allocated(this % label)) then
       name = this % label
    else
       name = ''
    end if

  end function counted_name

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

end module graph_carrier
