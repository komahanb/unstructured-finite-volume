!=====================================================================!
! Concrete graph partitioners.
!
! P cuts a graph into parts. One concrete type stores the rule, so a
! caller can store partitioners in a plain array and a new rule adds a
! case rather than a class.
!
!         o---o---o---o---o---o
!                     :                cut where few edges cross, so
!         o---o---o   :   o---o---o    few values pass between the
!                part 1     part 2     parts
!
! The output is one part, in its own numbering, and the part is still
! a graph. The part also records its relation to the whole - which
! cells it owns, which it reads as halo copies, and the global index of each
! of its own numbers.
!
!         global graph    1   2   3   4   5   6   7   8
!                                   |   |   |
!         part 2                    1   2   3
!
!                       global_vertex_index(2) = 4
!
!=====================================================================!
!
!                      OWNED, HALO, OVERLAP
!
! A part owns the cells it must produce values for. The part reads as halo copies
! the neighbouring cells that other parts own, because a face term
! needs the value on both sides. Together those are the overlap - the
! cells this part must read to complete what it owns.
!
!            part 1                        part 2
!       +---------------+            +---------------+
!       |  o    o    o  |            |  o    o    o  |
!       |  o    o    o--|------------|--b    o    o  |
!       +---------------+            +---------------+
!                    \______________/
!                       part 1 reads one halo cell from part 2
!
! Every cell is owned by exactly one part. That is what stops a
! conserved quantity being counted twice when the parts are added back
! together.
!
!=====================================================================!
!
!                        WHAT IS NOT HERE
!
! No physics, no geometry, no solver behaviour. A partitioner computes
! which cells belong to which part and restricts the data across the
! same cut. It evaluates nothing.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module transform_partitioner

  use util_precision  , only : dp
  use view_directed , only : directed_graph
  use relation_partition, only : partition_relation
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use map_set_store, only : set_store
  use transform_structure, only : transform
  use view_directed_stored         , only : stored_directed_graph
  use field_stored   , only : stored_field
  use operation_traversal    , only : traversal, TRAVERSAL_VISIT_ORDER

  implicit none

  private
  public :: partitioner
  public :: PARTITION_LINEAR, PARTITION_BREADTH_FIRST, PARTITION_ADOPTED

  !===================================================================!
  ! PARTITIONER. The transform that cuts a whole into parts.
  ! partition_graph returns the part and the relation
  ! r <= S_part x S_whole recording how the part embeds in the
  ! whole; the relation is required output, because a part without
  ! it cannot be reassembled. partition_data restricts a field onto
  ! a part along that same relation, writing the declared domain
  ! into the caller's set store.
  !===================================================================!

  !-------------------------------------------------------------------!
  ! Slice the vertex numbering into equal blocks. Low cost,
  ! deterministic, and independent of how the cells are connected -
  ! useful as a reference the other rules are compared against.
  !-------------------------------------------------------------------!

  integer, parameter :: PARTITION_LINEAR = 1

  !-------------------------------------------------------------------!
  ! Grow each part outward from a seed cell, one ring at a time, until
  ! the part has its share. The rule follows the connections, so it
  ! cuts fewer edges.
  !-------------------------------------------------------------------!

  integer, parameter :: PARTITION_BREADTH_FIRST = 2

  !-------------------------------------------------------------------!
  ! Take a map computed elsewhere - a mesh file, an outside library,
  ! a previous run.
  !-------------------------------------------------------------------!

  integer, parameter :: PARTITION_ADOPTED = 3

  !===================================================================!
  ! One partitioner: how to cut, into how many, and which part to
  ! return.
  !===================================================================!

  type, extends(transform) :: partitioner

     integer :: rule   = PARTITION_LINEAR
     integer :: num_parts = 1
     integer :: part   = 1

     !----------------------------------------------------------------!
     ! The map an adopted partition was given, one owning part per
     ! whole-graph vertex.
     !----------------------------------------------------------------!

     integer, allocatable :: adopted(:)

   contains

     procedure :: defined_on_graph
     procedure :: defined_on_data
     procedure :: partition_graph
     procedure :: partition_data

  end type partitioner

  interface partitioner
     module procedure create
  end interface partitioner

contains

  !===================================================================!
  ! Build a partitioner: the rule, the number of parts, and which
  ! part to return.
  !===================================================================!

  pure type(partitioner) function create(rule, num_parts, part, adopted) result(this)

    integer, intent(in)           :: rule
    integer, intent(in)           :: num_parts
    integer, intent(in), optional :: part
    integer, intent(in), optional :: adopted(:)

    this % rule   = rule
    this % num_parts = num_parts

    if (present(part))    this % part = part
    if (present(adopted)) allocate(this % adopted, source=adopted)

  end function create

  !===================================================================!
  ! A graph can be partitioned when it has vertices, and when an
  ! adopted map - if that is the rule - covers them.
  !===================================================================!

  pure logical function defined_on_graph(this, input_graph)

    class(partitioner), intent(in) :: this
    class(directed_graph)      , intent(in) :: input_graph

    defined_on_graph = input_graph % num_vertices() > 0 .and. this % num_parts >= 1

    if (this % rule == PARTITION_ADOPTED) then
       if (.not. allocated(this % adopted)) then
          defined_on_graph = .false.
       else if (size(this % adopted) < input_graph % num_vertices()) then
          defined_on_graph = .false.
       end if
    end if

  end function defined_on_graph

  !===================================================================!
  ! Data can be restricted when the graph can be partitioned and the
  ! data is defined on that graph.
  !===================================================================!

  logical function defined_on_data(this, input_graph, input_data)

    class(partitioner), intent(in) :: this
    class(directed_graph)      , intent(in) :: input_graph
    class(field) , intent(in) :: input_data

    defined_on_data = this % defined_on_graph(input_graph)

    select type (input_data)
    class is (field)
       defined_on_data = defined_on_data .and. input_data % num_entries() >= 0
    class default
       defined_on_data = .false.
    end select

  end function defined_on_data

  !===================================================================!
  ! P. Assign an owner to every cell, gather this part's cells, and
  ! build the part as a graph in its own numbering.
  !===================================================================!

  subroutine partition_graph(this, global_graph, part_graph, rel)

    class(partitioner), intent(in)               :: this
    class(directed_graph), intent(in)            :: global_graph
    class(directed_graph), allocatable, intent(out) :: part_graph

    ! r <= S_part x S_whole, written here and nowhere else. The
    ! relation is returned beside the part, because the three
    ! procedures that read it receive it as an argument - none of them
    ! may read it from the part, and none may construct its own.
    !
    ! REQUIRED. A cut produces the pair (G_p, r_p); taking the part
    ! and leaving the relation is taking half the output.
    type(partition_relation), intent(out) :: rel

    integer, allocatable :: owner(:), own_part(:), part_index(:)
    integer, allocatable :: ltail(:), lhead(:), eglobal(:), eowner(:), vowner(:)
    integer :: nv, ne, e, t, h, k, num_part_edges

    nv = global_graph % num_vertices()
    ne = global_graph % num_edges()

    call assign_owners(this, global_graph, owner)
    call gather_part(global_graph, owner, this % part, own_part, part_index)

    ! Retain an edge when both its ends are in this part and at least
    ! one of them is owned here. An edge with neither end owned belongs
    ! entirely to another part; retaining it here would add its flux to
    ! the balance twice.
    allocate(ltail(ne), lhead(ne), eglobal(ne), eowner(ne))
    num_part_edges = 0
    do e = 1, ne
       t = global_graph % edge_tail(e)
       h = global_graph % edge_head(e)

       if (part_index(t) == 0) cycle
       if (h >= 1) then
          if (part_index(h) == 0) cycle
       end if
       if (owner(t) /= this % part) then
          if (h < 1) cycle
          if (owner(h) /= this % part) cycle
       end if

       num_part_edges = num_part_edges + 1
       ltail(num_part_edges) = part_index(t)
       if (h >= 1) then
          lhead(num_part_edges) = part_index(h)
       else
          lhead(num_part_edges) = 0
       end if
       eglobal(num_part_edges) = e

       ! An edge is owned by the part that owns its TAIL - always,
       ! whether this part owns that tail or reads it as a halo copy. One global
       ! edge therefore has exactly one owner across all parts, which
       ! is what makes assembly reconstruct a global edge field
       ! exactly once. (The branch below is vestigial: both arms
       ! assign the same value. An earlier design let the head's
       ! owner own the edge for a halo tail; that rule was never
       ! implemented, and the uniqueness property does not need it.)
       if (owner(t) == this % part) then
          eowner(num_part_edges) = owner(t)
       else
          eowner(num_part_edges) = owner(t)
       end if
    end do

    allocate(vowner(size(own_part)))
    do k = 1, size(own_part)
       vowner(k) = owner(own_part(k))
    end do

    ! The part is constructed with r. The tuples are passed to the
    ! constructor, because a graph given its relation after
    ! construction could return two different relations in one
    ! lifetime.
    allocate(part_graph, source = &
         & stored_directed_graph(size(own_part), tails=ltail(1:num_part_edges), heads=lhead(1:num_part_edges), &
         &              number  = this % part,   &
         &              num_parts  = this % num_parts, &
         &              vglobal = own_part,          &
         &              vowner  = vowner,        &
         &              eglobal = eglobal(1:num_part_edges), &
         &              eowner  = eowner(1:num_part_edges),  &
         &              whole_vertices   = global_graph % vertex_set(), &
         &              whole_edges   = global_graph % edge_set(),   &
         &              num_whole_vertices = nv,                          &
         &              num_whole_edges = ne))

    select type (part_graph)
    class is (stored_directed_graph)
       rel = part_graph % whole_relation()
    end select

  end subroutine partition_graph

  !===================================================================!
  ! Assign the owning part of each cell of the whole graph.
  !===================================================================!

  subroutine assign_owners(this, global_graph, owner)

    class(partitioner)  , intent(in)  :: this
    class(directed_graph)        , intent(in)  :: global_graph
    integer, allocatable, intent(out) :: owner(:)

    integer :: nv

    nv = global_graph % num_vertices()
    allocate(owner(nv))

    select case (this % rule)

    case (PARTITION_ADOPTED)
       owner = this % adopted(1:nv)

    case (PARTITION_BREADTH_FIRST)
       call assign_owners_breadth_first(global_graph, this % num_parts, owner)

    case default
       call assign_owners_linear(nv, this % num_parts, owner)

    end select

  end subroutine assign_owners

  !===================================================================!
  ! Equal blocks of the numbering. The first few parts take one extra
  ! cell when the count does not divide evenly.
  !===================================================================!

  pure subroutine assign_owners_linear(nv, num_parts, owner)

    integer, intent(in)    :: nv, num_parts
    integer, intent(inout) :: owner(:)

    integer :: v, base, additional, lo, hi, k

    base  = nv / num_parts
    additional = mod(nv, num_parts)

    hi = 0
    do k = 1, num_parts
       lo = hi + 1
       hi = lo + base - 1
       if (k <= additional) hi = hi + 1
       do v = lo, min(hi, nv)
          owner(v) = k
       end do
    end do

  end subroutine assign_owners_linear

  !===================================================================!
  ! Extend every part from a seed by successive neighbourhoods, so each
  ! part is connected and few edges cross.
  !===================================================================!

  subroutine assign_owners_breadth_first(global_graph, num_parts, owner)

    class(directed_graph), intent(in)    :: global_graph
    integer     , intent(in)    :: num_parts
    integer     , intent(inout) :: owner(:)

    type(stored_directed_graph) :: unassigned_graph
    type(traversal)         :: visit
    class(field), allocatable :: traversal_order
    integer, allocatable :: global_vertices(:), part_index(:), tails(:), heads(:), order(:)
    integer :: nv, ne, max_part_size, k, v, e, t, h, n, m

    nv    = global_graph % num_vertices()
    ne    = global_graph % num_edges()
    owner = 0
    max_part_size = (nv + num_parts - 1) / num_parts

    allocate(global_vertices(nv), part_index(nv), tails(ne), heads(ne))

    do k = 1, num_parts

       ! The induced graph on the unassigned vertices. The traversal implements
       ! breadth-first; this routine only calls it and reads the result.
       n = 0
       part_index = 0
       do v = 1, nv
          if (owner(v) == 0) then
             n = n + 1
             global_vertices(n)  = v
             part_index(v) = n
          end if
       end do
       if (n == 0) exit

       m = 0
       do e = 1, ne
          t = global_graph % edge_tail(e)
          if (.not. global_graph % edge_has_head(e)) cycle
          h = global_graph % edge_head(e)
          if (part_index(t) > 0 .and. part_index(h) > 0) then
             m = m + 1
             tails(m) = part_index(t)
             heads(m) = part_index(h)
          end if
       end do

       unassigned_graph = stored_directed_graph(n, tails=tails(1:m), heads=heads(1:m))

       ! The first unassigned cell seeds the part; the visit order
       ! determines its members, up to max_part_size.
       visit = traversal(TRAVERSAL_VISIT_ORDER, seed=1)
       call visit % apply(unassigned_graph, output=traversal_order)
       call traversal_order % integer_vector(order)

       do v = 1, n
          if (order(v) >= 1 .and. order(v) <= max_part_size) owner(global_vertices(v)) = k
       end do

    end do

    ! Every cell the rings never reached is assigned to the last part,
    ! so every cell is owned exactly once.
    do v = 1, nv
       if (owner(v) == 0) owner(v) = num_parts
    end do

  end subroutine assign_owners_breadth_first

  !===================================================================!
  ! Collect one part's cells: the ones it owns first, then the ones it
  ! must read as halo copies to compute its own values.
  !
  ! The owned-first order has a practical consequence: a part's owned
  ! values are stored at the front of every vector, so the portion a
  ! solver reduces over is a contiguous slice.
  !===================================================================!

  subroutine gather_part(global_graph, owner, part, own_part, part_index)

    class(directed_graph)        , intent(in)  :: global_graph
    integer             , intent(in)  :: owner(:)
    integer             , intent(in)  :: part
    integer, allocatable, intent(out) :: own_part(:)
    integer, allocatable, intent(out) :: part_index(:)

    integer, allocatable :: nbrs(:)
    integer :: nv, v, i, n

    nv = global_graph % num_vertices()

    allocate(part_index(nv))
    part_index = 0
    allocate(own_part(nv))
    n = 0

    do v = 1, nv
       if (owner(v) == part) then
          n = n + 1
          own_part(n)    = v
          part_index(v) = n
       end if
    end do

    do v = 1, nv
       if (owner(v) /= part) cycle
       call global_graph % adjacent_vertices(v, nbrs)
       do i = 1, size(nbrs)
          if (owner(nbrs(i)) /= part .and. part_index(nbrs(i)) == 0) then
             n = n + 1
             own_part(n)          = nbrs(i)
             part_index(nbrs(i)) = n
          end if
       end do
    end do

    own_part = own_part(1:n)

  end subroutine gather_part

  !===================================================================!
  ! Restrict the data across the same cut, by the map the part graph
  ! already stores. Nothing here recomputes the cut, so the values
  ! cannot become inconsistent with the structure.
  !===================================================================!

  subroutine partition_data(this, rel, global_graph, global_data, part_graph, &
       & sets, part_data)

    class(partitioner), intent(in)               :: this
    type(partition_relation), intent(in)         :: rel
    class(directed_graph)      , intent(in)               :: global_graph
    class(field) , intent(in)               :: global_data
    class(directed_graph)      , intent(in)               :: part_graph
    type(set_store)    , intent(inout)            :: sets
    class(field) , allocatable, intent(out) :: part_data

    type(graph) :: dom
    integer         :: n_dom
    character(len=250) :: message

    associate (u1 => this); end associate

    ! ONE RELATION PER PART, checked on the restriction too. This
    ! is the direction where a wrong r is SILENT: values restricted
    ! onto a part by another part's numbering, no abort, wrong values
    ! near a cut - which is the failure the design exists to prevent.
    if (.not. rel % describes(part_graph)) then
       error stop 'partition: rel does not describe part_graph; rel % describes(part_graph) is false'
    end if

    select type (global_data)

    class is (stored_field)
       dom   = global_data % domain()
       n_dom = global_data % num_entries()
       ! Classify by embedding - a DECLARED predicate evaluated through
       ! the set store. Coverage selects the restriction inside.
       if (sets % subobject_of(dom, global_graph % vertex_set())) then
          call restrict_field(global_data, dom, n_dom, &
               & global_graph % vertex_set(), global_graph % num_vertices(), &
               & part_graph, rel, .true., sets, part_data)
       else if (sets % subobject_of(dom, global_graph % edge_set())) then
          call restrict_field(global_data, dom, n_dom, &
               & global_graph % edge_set(), global_graph % num_edges(), &
               & part_graph, rel, .false., sets, part_data)
       else
          write(message,'(a,i0,a)') 'partition: global_data has ', n_dom, &
               & ' entries, which match neither the vertex set nor the edge set of global_graph'
          error stop trim(message)
       end if

    class default
       error stop 'partition: global_data''s dynamic type is not handled by this transform; &
            &only class(stored_field) is'
    end select

  end subroutine partition_data

  !===================================================================!
  ! One restriction for both families and both coverages. A FULL
  ! field - domain same_as the global carrier - is defined on the
  ! part's own carrier, every part member valued through the global
  ! map. A PROPER SUBSET is restricted as a subset: the part-local
  ! members whose global indices the subset contains, each value read
  ! through the GLOBAL DOMAIN'S local_index - never by raw member
  ! arithmetic - and stored on a new subobject of the part's carrier.
  ! A new ambient means a new declared subset: identity is not
  ! preserved across restriction, extension and values are.
  !===================================================================!

  subroutine restrict_field(global_data, dom, n_dom, global_carrier, &
       &                 n_global_carrier, part_graph, rel, on_vertices, &
       &                 sets, part_data)

    type(stored_field)        , intent(in)               :: global_data
    type(graph)    , intent(in)               :: dom
    integer            , intent(in)               :: n_dom
    type(graph)    , intent(in)               :: global_carrier
    integer            , intent(in)               :: n_global_carrier
    class(directed_graph)       , intent(in)               :: part_graph
    type(partition_relation), intent(in)          :: rel
    logical            , intent(in)               :: on_vertices
    type(set_store)    , intent(inout)            :: sets
    class(field) , allocatable, intent(out) :: part_data

    type(stored_field)           :: out
    type(graph)       :: part_carrier
    type(graph)       :: sp
    real(dp), allocatable :: fv(:), lv(:)
    integer , allocatable :: part_members(:)
    integer :: nlocal, num_components, l, c, g, n, at

    if (on_vertices) then
       nlocal       = part_graph % num_vertices()
       part_carrier = part_graph % vertex_set()
    else
       nlocal       = part_graph % num_edges()
       part_carrier = part_graph % edge_set()
    end if
    num_components = global_data % num_components()

    call global_data % real_vector(fv)

    if (dom % same_as(global_carrier)) then

       ! Full coverage: the part field is defined on the part's carrier.
       allocate(lv(nlocal * num_components))
       lv = 0.0_dp
       do l = 1, nlocal
          g = rel % global_index(l, on_vertices)
          at = sets % index_in(dom, g)
          if (at >= 1) then
             do c = 1, num_components
                lv((l - 1) * num_components + c) = fv((at - 1) * num_components + c)
             end do
          end if
       end do
       out = stored_field(global_data % name(), part_carrier, nlocal, num_components=num_components, &
            &      unit_name=global_data % units())
       call out % set_real_vector(lv)

    else

       ! Proper subset: gather the part members the subset names.
       allocate(part_members(nlocal))
       n = 0
       do l = 1, nlocal
          g = rel % global_index(l, on_vertices)
          if (sets % has(dom, g)) then
             n = n + 1
             part_members(n) = l
          end if
       end do
       !-------------------------------------------------------------!
       ! A new ambient means a new declared subset, so this operation
       ! satisfies the subobject invariant: identity, extension, label
       ! and embedding, together. The label is the global domain's
       ! label - restriction renames nothing.
       !-------------------------------------------------------------!

       call sets % declare_subobject(sp, part_members(1:n), sets % label_of(dom), part_carrier)

       allocate(lv(n * num_components))
       do l = 1, n
          g  = rel % global_index(part_members(l), on_vertices)
          at = sets % index_in(dom, g)
          do c = 1, num_components
             lv((l - 1) * num_components + c) = fv((at - 1) * num_components + c)
          end do
       end do
       out = stored_field(global_data % name(), sp, n, num_components=num_components, &
            &      unit_name=global_data % units())
       call out % set_real_vector(lv)

    end if

    allocate(part_data, source=out)

  end subroutine restrict_field

end module transform_partitioner
