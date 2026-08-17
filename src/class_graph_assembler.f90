!=====================================================================!
! The concrete graph assembler.
!
! P inverse, and only that. It puts a piece back into whole-graph
! order and brings its data with it.
!
!         part 2   1   2   3
!                  |   |   |
!         whole    3   4   5           by the relation the cut wrote.
!                                      The assembler is HANDED r; it
!                                      never invents one and never
!                                      keeps one.
!
! The law it has to satisfy:
!
!         assemble( partition( G ) )     ==  G
!         assemble( partition( G, D ) )  ==  ( G, D )
!
!=====================================================================!
!
!                   ONLY OWNED VALUES ARE COLLECTED
!
! A part borrows the cells around its edge so it can work out its own
! answers. A borrowed value is a copy of a value another part owns.
! Collecting both copies counts a conserved quantity twice - mass
! appears from nowhere, and it appears only in parallel, only near a
! partition boundary, where such an error is hardest to locate.
!
!            part 1                        part 2
!       +---------------+            +---------------+
!       |  o    o    o  |            |  o    o    o  |
!       |  o    o    O--|------------|--b    o    o  |
!       +---------------+            +---------------+
!                       \____________/
!                    part 1 borrows this cell.
!                    part 2 owns it and answers for it.
!                    exactly one of them is collected.
!
!=====================================================================!
!
!                     WHAT ONE PART CAN AND CANNOT DO
!
! The contract hands the assembler a single part. So it can restore
! everything that part owns, and it cannot invent what it never saw.
!
! With one part, that is the whole graph and the round trip is exact.
! With several, each call fills in that part's own share and leaves
! the rest alone, so summing the answers over all the parts rebuilds
! the whole. The union of the owned sets is the whole graph and the
! sets do not overlap, which is what makes that sum right.
!
! ASSEMBLER MEANS THIS AND NOTHING ELSE. No physics, no boundary
! conditions, no residual, no matrix, no file, no solver behaviour.
! Any of those added here obscures the one-line law the type exists
! to keep visible.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_graph_assembler

  use iso_fortran_env     , only : dp => REAL64
  use graph_directed_view , only : directed_graph
  use graph_partition_relation, only : partition_relation
  use graph_field_calculus, only : graph_field
  use fractal_graph      , only : set_graph => graph
  use graph_set_map      , only : set_map
  use graph_label_map    , only : label_map
  use graph_inclusion_map, only : inclusion_map, declared_subobject
  use graph_set_representation, only : listed_set_representation
  use graph_calculus      , only : graph_assembler
  use class_graph         , only : directed_stored_graph
  use class_graph_field   , only : field

  implicit none

  private
  public :: assembler

  !===================================================================!
  ! THE ASSEMBLER HOLDS NOTHING.
  !
  ! It briefly held a bound relation, and that was one holder too
  ! many. Partition and assembly are opposite verbs over ONE r, so r
  ! belongs to neither of them: the cut writes it, and every verb
  ! that reads it is handed it. An assembler that carried its own
  ! copy could be paired with a part the copy was never written for,
  ! and the law it exists to keep - one relation per part - would be
  ! a convention rather than an argument.
  !
  ! So: no state, no set_map, no label_map, no inclusion_map. What a
  ! set MEANS is the caller's, and arrives at the semantic boundary -
  ! as arguments to assemble_data, every time. WHERE a member goes
  ! arrives beside it, as r.
  !===================================================================!

  type, extends(graph_assembler) :: assembler

   contains

     procedure :: defined_on_graph
     procedure :: defined_on_data
     procedure :: defined_on_relation
     procedure :: assemble_graph
     procedure :: assemble_data

  end type assembler

contains

  !===================================================================!
  ! THE TRANSFORM'S GENERIC GATE, AND IT IS A WEAK ONE ON PURPOSE.
  !
  ! With no relation in hand there is exactly one thing an assembler
  ! can say about a bare graph: whether it has members to put back.
  ! The real question - is this the part r was written for - needs r,
  ! and r is not an argument here because graph_transform's contract
  ! is not about relations. It is defined_on_relation below, and that
  ! is the gate a caller assembling parts must use.
  !===================================================================!

  pure logical function defined_on_graph(this, input_graph)

    class(assembler)     , intent(in) :: this
    class(directed_graph), intent(in) :: input_graph

    associate (u1 => this); end associate

    defined_on_graph = input_graph % num_vertices() > 0

  end function defined_on_graph

  !===================================================================!
  ! THE REAL GATE: is this r the one written for this part?
  !
  ! A PREDICATE, not an abort. A caller holds one relation per part -
  ! that is the law - so handing over the wrong one is a mistake that
  ! has to be REPORTABLE. Ending the image would leave nothing to
  ! report it to.
  !===================================================================!

  logical function defined_on_relation(this, rel, part_graph)

    class(assembler)        , intent(in) :: this
    type(partition_relation), intent(in) :: rel
    class(directed_graph)   , intent(in) :: part_graph

    associate (u1 => this); end associate

    defined_on_relation = rel % describes(part_graph)

  end function defined_on_relation

  !===================================================================!
  ! Can this assembler say anything about that data? Yes for a field
  ! riding a part that passes the graph gate above.
  !===================================================================!

  ! Reads through defined_on_graph, so it inherits that gate's reach
  ! and no more. The relation question is defined_on_relation.
  logical function defined_on_data(this, input_graph, input_data)

    class(assembler) , intent(in) :: this
    class(directed_graph)     , intent(in) :: input_graph
    class(graph_field), intent(in) :: input_data

    defined_on_data = this % defined_on_graph(input_graph)

    select type (input_data)
    class is (graph_field)
       defined_on_data = defined_on_data .and. input_data % num_entries() >= 0
    class default
       defined_on_data = .false.
    end select

  end function defined_on_data

  !===================================================================!
  ! Put the piece back in whole-graph order.
  !
  ! Every cell of the part is renamed to what the whole graph called
  ! it, and every edge with it. What comes out is a graph again,
  ! and holds no partition record - because a whole graph is not a
  ! part of anything.
  !===================================================================!

  subroutine assemble_graph(this, rel, part_graph, global_graph)

    class(assembler), intent(in)               :: this
    type(partition_relation), intent(in)       :: rel
    class(directed_graph)    , intent(in)               :: part_graph
    class(directed_graph)    , allocatable, intent(out) :: global_graph

    integer, allocatable :: tails(:), heads(:)
    integer :: ne, e, nv_global, l, biggest

    if (.not. this % defined_on_relation(rel, part_graph)) then
       error stop 'assemble: this relation was not written for this part'
    end if

    ne = part_graph % num_edges()

    ! The whole graph is at least as big as the largest whole-graph
    ! index r records.
    biggest = 0
    do l = 1, part_graph % num_vertices()
       biggest = max(biggest, rel % global_vertex_index(l))
    end do
    nv_global = max(biggest, rel % num_whole_vertices())

    allocate(tails(ne), heads(ne))
    do e = 1, ne
       tails(e) = rel % global_vertex_index(part_graph % edge_tail(e))
       if (part_graph % edge_has_head(e)) then
          heads(e) = rel % global_vertex_index(part_graph % edge_head(e))
       else
          heads(e) = 0
       end if
    end do

    allocate(global_graph, source = &
         & directed_stored_graph(nv_global, tails=tails, heads=heads, &
         &                       number=rel % part_id()))

  end subroutine assemble_graph

  !===================================================================!
  ! Assemble the data back onto the whole graph.
  !
  ! The answer is laid out on the whole graph. Only the entries this
  ! part owns are written; everything else is left at zero, so adding
  ! the answers from every part rebuilds the whole field exactly once.
  !===================================================================!

  subroutine assemble_data(this, rel, part_graph, part_data, global_graph, &
       & sets, labels, inclusions, global_data)

    class(assembler) , intent(in)               :: this
    type(partition_relation), intent(in)        :: rel
    class(directed_graph)     , intent(in)               :: part_graph
    class(graph_field), intent(in)               :: part_data
    class(directed_graph)     , intent(in)               :: global_graph
    type(set_map)      , intent(inout)            :: sets
    type(label_map)    , intent(inout)            :: labels
    type(inclusion_map), intent(inout)            :: inclusions
    class(graph_field), allocatable, intent(out) :: global_data

    type(set_graph) :: dom
    integer         :: n_dom

    if (.not. this % defined_on_relation(rel, part_graph)) then
       error stop 'assemble: this relation was not written for this part'
    end if

    select type (part_data)

    class is (field)
       dom   = part_data % domain()
       n_dom = part_data % num_entries()
       ! Classify by embedding - a DECLARED question, so the inclusion
       ! map answers it, never the extension and never the graph.
       if (declared_subobject(dom, part_graph % vertex_set(), inclusions)) then
          call gather_field(part_data, dom, n_dom, part_graph, rel, &
               & part_graph % vertex_set(), part_graph % num_vertices(), &
               & global_graph, .true., sets, labels, inclusions, global_data)
       else if (declared_subobject(dom, part_graph % edge_set(), inclusions)) then
          call gather_field(part_data, dom, n_dom, part_graph, rel, &
               & part_graph % edge_set(), part_graph % num_edges(), &
               & global_graph, .false., sets, labels, inclusions, global_data)
       else
          error stop 'assemble: this field does not live on this part''s domains'
       end if

    class default
       error stop 'assemble: this data does not ride on this transform'
    end select

  end subroutine assemble_data

  !===================================================================!
  ! One gather for both families and both coverages. A FULL part
  ! field lands on the GLOBAL carrier, owned members only, exactly
  ! the established assembly. A PROPER SUBSET maps home through the
  ! part->global map and lands on a new subobject of the global
  ! carrier - its actual mapped subdomain, no manufactured zeros on
  ! members the field never held. A new ambient means a new
  ! declared subset: extension and values return, tokens do not.
  !===================================================================!

  subroutine gather_field(part_data, dom, n_dom, part_graph, rel, part_carrier, &
       &                  n_part_carrier, global_graph, on_vertices, &
       &                  sets, labels, inclusions, global_data)

    type(field)        , intent(in)               :: part_data
    type(set_graph)    , intent(in)               :: dom
    integer            , intent(in)               :: n_dom
    class(directed_graph)       , intent(in)               :: part_graph
    type(partition_relation), intent(in)          :: rel
    type(set_graph)    , intent(in)               :: part_carrier
    integer            , intent(in)               :: n_part_carrier
    class(directed_graph)       , intent(in)               :: global_graph
    logical            , intent(in)               :: on_vertices
    type(set_map)      , intent(inout)            :: sets
    type(label_map)    , intent(inout)            :: labels
    type(inclusion_map), intent(inout)            :: inclusions
    class(graph_field) , allocatable, intent(out) :: global_data

    type(field)           :: out
    type(set_graph)       :: global_carrier
    type(set_graph)       :: sg
    real(dp), allocatable :: lv(:), fv(:)
    integer , allocatable :: kept(:), came(:)
    integer :: nglobal, nlocal, ncomp, l, c, f, me, n, at

    if (on_vertices) then
       nglobal        = global_graph % num_vertices()
       nlocal         = part_graph % num_vertices()
       global_carrier = global_graph % vertex_set()
    else
       nglobal        = global_graph % num_edges()
       nlocal         = part_graph % num_edges()
       global_carrier = global_graph % edge_set()
    end if
    ncomp = part_data % num_components()
    me    = rel % part_id()

    call part_data % get_real_vector(lv)

    if (dom % same_as(part_carrier)) then

       ! Full coverage: the established dense assembly, owned only.
       out = field(part_data % name(), global_carrier, nglobal, ncomp=ncomp, &
            &      unit_name=part_data % units())
       allocate(fv(nglobal * ncomp))
       fv = 0.0_dp

       do l = 1, nlocal
          if (rel % has_part_relation()) then
             if (owner_of(rel, l, on_vertices) /= me) cycle
          end if
          f = global_of(rel, l, on_vertices)
          do c = 1, ncomp
             associate (to => (f - 1) * ncomp + c, from => (l - 1) * ncomp + c)
               if (to >= 1 .and. to <= size(fv) .and. from <= size(lv)) fv(to) = lv(from)
             end associate
          end do
       end do

       call out % set_real_vector(fv)

    else

       ! Proper subset: carry the members home and keep only them.
       allocate(kept(n_dom), came(n_dom))
       n = 0
       do l = 1, n_dom
          at = sets % member_of(dom, l)      ! part-local member
          if (rel % has_part_relation()) then
             if (owner_of(rel, at, on_vertices) /= me) cycle
          end if
          n = n + 1
          kept(n) = global_of(rel, at, on_vertices)
          came(n) = l
       end do
       !-------------------------------------------------------------!
       ! A new ambient means a new declared subset, so this is a carve
       ! and obeys the carve law: identity, extension, label and
       ! embedding, together. Extension and values return home; tokens
       ! do not, and the label does.
       !-------------------------------------------------------------!

       call sg % declare()
       call sets       % bind(sg, listed_set_representation(kept(1:n)))
       call labels     % bind(sg, labels % label_of(dom))
       call inclusions % include_in(sg, global_carrier)

       allocate(fv(n * ncomp))
       do l = 1, n
          do c = 1, ncomp
             fv((l - 1) * ncomp + c) = lv((came(l) - 1) * ncomp + c)
          end do
       end do
       out = field(part_data % name(), sg, n, ncomp=ncomp, &
            &      unit_name=part_data % units())
       call out % set_real_vector(fv)

    end if

    allocate(global_data, source=out)

  end subroutine gather_field

  pure integer function global_of(rel, l, on_vertices)

    type(partition_relation), intent(in) :: rel
    integer                 , intent(in) :: l
    logical                 , intent(in) :: on_vertices

    if (on_vertices) then
       global_of = rel % global_vertex_index(l)
    else
       global_of = rel % global_edge_index(l)
    end if

  end function global_of

  pure integer function owner_of(rel, l, on_vertices)

    type(partition_relation), intent(in) :: rel
    integer                 , intent(in) :: l
    logical                 , intent(in) :: on_vertices

    if (on_vertices) then
       owner_of = rel % vertex_owner_part(l)
    else
       owner_of = rel % edge_owner_part(l)
    end if

  end function owner_of

end module class_graph_assembler
