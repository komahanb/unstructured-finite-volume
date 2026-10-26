!=====================================================================!
! Module that keeps, in one owner, the data attached to set
! identities: a set's members, its name and the set it is declared
! inside.
!
!                            set_store
!   +----------------------------------------------------------+
!   |  set_map (extended)    label_map       inclusion_map     |
!   |  set -> its members    set -> name     part -> ambient   |
!   +----------------------------------------------------------+
!      bind, extent, ...     name            include_in
!                            labelled        subobject_of
!                            label_of
!
! The store is the set map, extended by a label map and an inclusion
! map, so that a caller that reads or declares sets need not pass
! three side tables through every signature. The three maps stay
! distinct, each behind its own narrow procedures; only their three
! public interfaces become one.
!
! Every stored key is a copied token, owned by the map that stores
! it. No graph pointer is stored, and no TARGET argument is needed to
! recognize a set later.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_set_store

  use graph_fractal          , only : graph
  use map_set                , only : set_map
  use map_label              , only : label_map
  use map_inclusion          , only : inclusion_map, declared_subobject
  use map_set_representation , only : listed_set_representation

  implicit none

  private

  public :: set_store

  !-------------------------------------------------------------------!
  ! Type that models a set map with a label map and an inclusion map
  ! beside it. Every query and binding of members is the set map's
  ! own; the store adds the name, the declared inclusion and the
  ! subobject order.
  !-------------------------------------------------------------------!

  type, extends(set_map) :: set_store

     type(label_map)     , private :: labels
     type(inclusion_map) , private :: inclusions

   contains

     ! Declaring
     procedure :: name
     procedure :: include_in
     procedure :: declare_subobject

     ! Names
     procedure :: labelled
     procedure :: label_of

     ! The subobject order
     procedure :: subobject_of

  end type set_store

contains

  !===================================================================!
  ! Names a set. A set is named once; the label map rejects a second
  ! name.
  !===================================================================!

  subroutine name(this, element, label)

    class(set_store), intent(inout) :: this
    type(graph)     , intent(in)    :: element
    character(len=*), intent(in)    :: label

    call this % labels % bind(element, label)

  end subroutine name

  !===================================================================!
  ! Declares part inside ambient, part c--> ambient. A set is declared
  ! into one ambient only; the inclusion map rejects a second.
  !===================================================================!

  subroutine include_in(this, part, ambient)

    class(set_store), intent(inout) :: this
    type(graph)     , intent(in)    :: part
    type(graph)     , intent(in)    :: ambient

    call this % inclusions % include_in(part, ambient)

  end subroutine include_in

  !===================================================================!
  ! Declares a new set of listed members, names it and places it
  ! inside ambient:
  !
  !   members = a new set
  !      |
  !      |  bind         its listed members      set map
  !      |  name         its label               label map
  !      v  include_in   the ambient             inclusion map
  !===================================================================!

  subroutine declare_subobject(this, members, listed_members, label, ambient)

    class(set_store), intent(inout) :: this
    type(graph)     , intent(out)   :: members
    integer         , intent(in)    :: listed_members(:)
    character(len=*), intent(in)    :: label
    type(graph)     , intent(in)    :: ambient

    call members % declare()
    call this % bind(members, listed_set_representation(listed_members))
    call this % name(members, label)
    call this % include_in(members, ambient)

  end subroutine declare_subobject

  !===================================================================!
  ! Returns whether a name was bound to the set
  !===================================================================!

  logical function labelled(this, element)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: element

    labelled = this % labels % labelled(element)

  end function labelled

  !===================================================================!
  ! Returns the set's name, or '' when none was bound
  !===================================================================!

  function label_of(this, element) result(label)

    class(set_store), intent(in)  :: this
    type(graph)     , intent(in)  :: element
    character(len=:), allocatable :: label

    label = this % labels % label_of(element)

  end function label_of

  !===================================================================!
  ! Returns whether part lies inside ancestor in the declared
  ! subobject order: part itself, or reached from it by following
  ! declared inclusions,
  !
  !   part c--> A c--> B c--> ... c--> ancestor
  !===================================================================!

  logical function subobject_of(this, part, ancestor)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: part
    type(graph)     , intent(in) :: ancestor

    subobject_of = declared_subobject(part, ancestor, this % inclusions)

  end function subobject_of

end module map_set_store
