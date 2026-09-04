!=====================================================================!
! SET STORE
!
! One owner for the data attached to set identities: extension, name
! and declared inclusion. The store IS the set map, extended by a
! label map and an inclusion map. Callers that only read or declare
! sets need not pass three side tables through every signature.
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
  use map_set_representation , only : listed_set_representation

  implicit none

  private
  public :: set_store

  !===================================================================!
  ! A set map, with a label map and an inclusion map beside it. Every
  ! extent query and bind is the set map's own; the store adds the
  ! name, the declared inclusion and the subobject order.
  !===================================================================!

  type, extends(set_map) :: set_store

     type(label_map)     , private :: labels
     type(inclusion_map) , private :: inclusions

   contains

     procedure :: name
     procedure :: include_in
     procedure :: declare_subobject

     procedure :: labelled
     procedure :: label_of

     procedure :: subobject_of

  end type set_store

contains

  subroutine name(this, element, label)

    class(set_store), intent(inout) :: this
    type(graph)     , intent(in)    :: element
    character(len=*), intent(in)    :: label

    call this % labels % bind(element, label)

  end subroutine name

  subroutine include_in(this, part, ambient)

    class(set_store), intent(inout) :: this
    type(graph)     , intent(in)    :: part
    type(graph)     , intent(in)    :: ambient

    call this % inclusions % include_in(part, ambient)

  end subroutine include_in

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

  logical function labelled(this, element)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: element

    labelled = this % labels % labelled(element)

  end function labelled

  function label_of(this, element) result(label)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: element
    character(len=:), allocatable :: label

    label = this % labels % label_of(element)

  end function label_of

  logical function subobject_of(this, part, ancestor)

    class(set_store), intent(in) :: this
    type(graph)     , intent(in) :: part
    type(graph)     , intent(in) :: ancestor

    subobject_of = declared_subobject(part, ancestor, this % inclusions)

  end function subobject_of

end module map_set_store
