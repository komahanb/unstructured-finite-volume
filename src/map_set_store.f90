!=====================================================================!
! SET STORE
!
! One owner for the data attached to set identities: extension, name
! and declared inclusion. The separate maps remain the storage
! detail. Callers that only read or declare sets need not pass three
! side tables through every signature.
!
! The store does not merge labels, extents and inclusions into one
! attribute system. It keeps the three maps distinct and gives each a
! narrow procedure. What is removed is the public interface to three
! maps, not the meaning.
!
! Lifetime law: every stored key is a copied token owned by the map
! below. No graph pointer is stored here, and no TARGET argument is
! required to recognize a set later.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_set_store

  use graph_fractal          , only : graph
  use map_set                , only : set_map
  use map_label              , only : label_map
  use map_inclusion          , only : inclusion_map, declared_subobject
  use map_set_representation , only : set_representation, listed_set_representation

  implicit none

  private
  public :: set_store

  type :: set_store

     type(set_map)       , private :: extents
     type(label_map)     , private :: labels
     type(inclusion_map) , private :: inclusions

   contains

     procedure :: bind
     procedure :: name
     procedure :: include_in
     procedure :: declare_subobject

     procedure :: describes
     procedure :: labelled
     procedure :: label_of

     procedure :: num_members_of
     procedure :: member_of
     procedure :: members_of
     procedure :: has
     procedure :: index_in
     procedure :: extent_of

     procedure :: included
     procedure :: declared_into
     procedure :: subobject_of

  end type set_store

contains

  subroutine bind(this, element, representation)

    class(set_store)          , intent(inout) :: this
    type(graph)               , intent(in)    :: element
    class(set_representation) , intent(in)    :: representation

    call this % extents % bind(element, representation)

  end subroutine bind

  subroutine name(this, element, text)

    class(set_store), intent(inout) :: this
    type(graph)     , intent(in)    :: element
    character(len=*), intent(in)    :: text

    call this % labels % bind(element, text)

  end subroutine name

  subroutine include_in(this, part, ambient)

    class(set_store), intent(inout) :: this
    type(graph)     , intent(in)    :: part
    type(graph)     , intent(in)    :: ambient

    call this % inclusions % include_in(part, ambient)

  end subroutine include_in

  subroutine declare_subobject(this, members, listed_members, text, ambient)

    class(set_store), intent(inout) :: this
    type(graph)     , intent(out)   :: members
    integer         , intent(in)    :: listed_members(:)
    character(len=*), intent(in)    :: text
    type(graph)     , intent(in)    :: ambient

    call members % declare()
    call this % bind(members, listed_set_representation(listed_members))
    call this % name(members, text)
    call this % include_in(members, ambient)

  end subroutine declare_subobject

  logical function describes(this, element)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: element

    describes = this % extents % describes(element)

  end function describes

  logical function labelled(this, element)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: element

    labelled = this % labels % labelled(element)

  end function labelled

  function label_of(this, element) result(text)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: element
    character(len=:), allocatable :: text

    text = this % labels % label_of(element)

  end function label_of

  integer function num_members_of(this, element)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: element

    num_members_of = this % extents % num_members_of(element)

  end function num_members_of

  integer function member_of(this, element, position)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: element
    integer         , intent(in) :: position

    member_of = this % extents % member_of(element, position)

  end function member_of

  subroutine members_of(this, element, values)

    class(set_store)   , intent(in)  :: this
    type(graph)        , intent(in)  :: element
    integer, allocatable, intent(out) :: values(:)

    call this % extents % members_of(element, values)

  end subroutine members_of

  logical function has(this, element, value)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: element
    integer         , intent(in) :: value

    has = this % extents % has(element, value)

  end function has

  integer function index_in(this, element, value)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: element
    integer         , intent(in) :: value

    index_in = this % extents % index_in(element, value)

  end function index_in

  subroutine extent_of(this, element, extent)

    class(set_store)                     , intent(in)  :: this
    type(graph)                          , intent(in)  :: element
    class(set_representation), allocatable, intent(out) :: extent

    call this % extents % extent_of(element, extent)

  end subroutine extent_of

  logical function included(this, part)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: part

    included = this % inclusions % included(part)

  end function included

  logical function declared_into(this, part, ambient)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: part
    type(graph)     , intent(in) :: ambient

    declared_into = this % inclusions % declared_into(part, ambient)

  end function declared_into

  logical function subobject_of(this, part, ancestor)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: part
    type(graph)     , intent(in) :: ancestor

    subobject_of = declared_subobject(part, ancestor, this % inclusions)

  end function subobject_of

end module map_set_store
