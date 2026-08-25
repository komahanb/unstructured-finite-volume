!=====================================================================!
! Carving: a declared subset minted with its extension, its label
! and its embedding bound together, in one place, so that no half-
! described set escapes. Every module that carves a subobject of a
! carrier - the directed view, the relation algorithms, the
! transports - calls this and states the law nowhere else.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module map_carving

  use graph_fractal         , only : graph
  use map_set               , only : set_map
  use map_label             , only : label_map
  use map_inclusion         , only : inclusion_map
  use map_set_representation, only : listed_set_representation

  implicit none

  private
  public :: carve

contains

  !===================================================================!
  ! CARVE. The one gate every named subset passes through.
  !
  ! A carved set is a NEW set - it signs a fresh identity, exactly as
  ! the subset_set it replaces did - and three things must be said
  ! about it or it is not usable:
  !
  !     its extension     which members, in this order
  !     its label         what the old subset called itself
  !     its embedding     which carrier it was carved from
  !
  ! They are bound together HERE rather than at each of the twelve
  ! call sites, because the third is the one an author forgets: a
  ! missing representation stops the program at the first query, and a
  ! missing label answers '', but a missing inclusion answers FALSE to
  ! is_subobject_of - quietly, and only on a real mesh.
  !
  ! The maps are the caller's. This routine writes into them and keeps
  ! nothing.
  !===================================================================!

  subroutine carve(members, roll, label, ambient, sets, labels, inclusions)

    type(graph)    , intent(out)   :: members
    integer            , intent(in)    :: roll(:)
    character(len=*)   , intent(in)    :: label
    type(graph)    , intent(in)    :: ambient
    type(set_map)      , intent(inout) :: sets
    type(label_map)    , intent(inout) :: labels
    type(inclusion_map), intent(inout) :: inclusions

    call members % declare()

    call sets       % bind(members, listed_set_representation(roll))
    call labels     % bind(members, label)
    call inclusions % include_in(members, ambient)

  end subroutine carve

end module map_carving
