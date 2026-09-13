!=====================================================================!
! The concrete graph assembler.
!
! P inverse, and only that. It maps a part back into whole-graph
! order and maps its data with it.
!
!         part 2   1   2   3
!                  |   |   |
!         whole    3   4   5           by the relation the partition
!                                      wrote. The assembler is PASSED
!                                      r; it never constructs one and
!                                      never stores one.
!
! The law it must satisfy:
!
!         assemble( partition( G ) )     ==  G
!         assemble( partition( G, D ) )  ==  ( G, D )
!
!=====================================================================!
!
!                   ONLY OWNED VALUES ARE COLLECTED
!
! A part reads halo copies of the cells along its boundary so it can compute its
! own values. A halo value is a copy of a value another part owns.
! Collecting both copies counts a conserved quantity twice - mass is
! created, and only in parallel, only near a partition boundary, where
! such an error is hardest to locate.
!
!            part 1                        part 2
!       +---------------+            +---------------+
!       |  o    o    o  |            |  o    o    o  |
!       |  o    o    O--|------------|--b    o    o  |
!       +---------------+            +---------------+
!                       \____________/
!                    part 1 reads this cell as a halo copy.
!                    part 2 owns it and reports its value.
!                    exactly one of them is collected.
!
!=====================================================================!
!
!                     WHAT ONE PART CAN AND CANNOT DO
!
! The contract passes the assembler a single part. So it can restore
! everything that part owns, and it cannot construct what it was never
! passed.
!
! With one part, that is the whole graph and the round trip is exact.
! With several, each call fills in that part's own subset and leaves
! the rest unchanged, so summing the results over all the parts
! rebuilds the whole. The union of the owned sets is the whole graph
! and the sets do not overlap, which is what makes that sum correct.
!
! ASSEMBLER MEANS THIS AND NOTHING ELSE. No physics, no boundary
! conditions, no residual, no matrix, no file, no solver behaviour.
! Any of those added here obscures the one-line law the type exists
! to keep visible.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module transform_assembler

  use util_precision  , only : dp
  use view_directed , only : directed_graph
  use relation_partition, only : partition_relation
  use field_calculus, only : field
  use graph_fractal      , only : graph
  use map_set_store, only : set_store
  use transform_structure, only : transform
  use view_directed_stored         , only : stored_directed_graph
  use field_stored   , only : stored_field

  implicit none

  private
  public :: assembler

  !===================================================================!
  ! ASSEMBLER. The inverse transform: parts become the whole again,
  ! along the same relation the partitioner wrote -
  ! assemble(partition(G)) = G. defined_on_relation reports
  ! whether a given relation belongs to a given part, so a caller
  ! storing several relations is informed that it passed the wrong
  ! one instead of receiving a wrong assembly. Only owned values are
  ! collected; counting a halo copy twice would violate
  ! conservation.
  !===================================================================!
  ! THE ASSEMBLER STORES NO STATE.
  !
  ! An earlier version stored a bound relation, and that was one
  ! owner too many. Partition and assembly are inverse operations
  ! over ONE r, so r belongs to neither of them: the partition writes
  ! it, and every operation that reads it is passed it. An assembler
  ! that stored its own copy could be paired with a part the copy was
  ! never written for, and the law it exists to enforce - one
  ! relation per part - would be a convention rather than an
  ! argument.
  !
  ! So: no state, no set_map, no label_map, no inclusion_map. What a
  ! set denotes is decided by the caller, and arrives at the semantic
  ! boundary - as arguments to assemble_data, every time. WHERE a
  ! member maps to arrives beside it, as r.
  !===================================================================!

  type, extends(transform) :: assembler

   contains

     procedure :: defined_on_graph
     procedure :: defined_on_data
     procedure :: defined_on_relation
     procedure :: assemble_graph
     procedure :: assemble_data

  end type assembler

contains

  !===================================================================!
  ! THE TRANSFORM'S GENERIC CHECK, AND IT IS A WEAK ONE BY DESIGN.
  !
  ! With no relation passed there is exactly one predicate an
  ! assembler can evaluate on a bare graph: whether it has members to
  ! map back. The complete predicate - is this the part r was written
  ! for - requires r, and r is not an argument here because
  ! transform's contract does not involve relations. That predicate
  ! is defined_on_relation below, and that is the check a caller
  ! assembling parts must use.
  !===================================================================!

  pure logical function defined_on_graph(this, input_graph)

    class(assembler)     , intent(in) :: this
    class(directed_graph), intent(in) :: input_graph

    associate (u1 => this); end associate

    defined_on_graph = input_graph % num_vertices() > 0

  end function defined_on_graph

  !===================================================================!
  ! THE COMPLETE CHECK: is this r the one written for this part?
  !
  ! A PREDICATE, not an error stop. A caller stores one relation per
  ! part - that is the law - so passing the wrong one is an error
  ! that must be REPORTABLE. Stopping the program would leave no
  ! caller to report it to.
  !===================================================================!

  logical function defined_on_relation(this, rel, part_graph)

    class(assembler)        , intent(in) :: this
    type(partition_relation), intent(in) :: rel
    class(directed_graph)   , intent(in) :: part_graph

    associate (u1 => this); end associate

    defined_on_relation = rel % describes(part_graph)

  end function defined_on_relation

  !===================================================================!
  ! Whether this assembler is defined on that data: true for a field
  ! on a part that passes the graph check above.
  !===================================================================!

  ! Evaluates through defined_on_graph, so it inherits that check's
  ! scope and no more. The relation predicate is defined_on_relation.
  logical function defined_on_data(this, input_graph, input_data)

    class(assembler) , intent(in) :: this
    class(directed_graph)     , intent(in) :: input_graph
    class(field), intent(in) :: input_data

    defined_on_data = this % defined_on_graph(input_graph)

    select type (input_data)
    class is (field)
       defined_on_data = defined_on_data .and. input_data % num_entries() >= 0
    class default
       defined_on_data = .false.
    end select

  end function defined_on_data

  !===================================================================!
  ! Map the part back to whole-graph order.
  !
  ! Every cell of the part is renumbered to its whole-graph index,
  ! and every edge with it. The result is a graph again, and stores
  ! no partition record - because a whole graph is not a part of
  ! anything.
  !===================================================================!

  subroutine assemble_graph(this, rel, part_graph, global_graph)

    class(assembler), intent(in)               :: this
    type(partition_relation), intent(in)       :: rel
    class(directed_graph)    , intent(in)               :: part_graph
    class(directed_graph)    , allocatable, intent(out) :: global_graph

    integer, allocatable :: tails(:), heads(:)
    integer :: ne, e, nv_global, l, largest_entry
    character(len=250) :: message

    if (.not. this % defined_on_relation(rel, part_graph)) then
       write(message,'(a,i0)') 'assemble_graph: rel was not written for this part; &
            &rel % part_id() = ', rel % part_id()
       error stop trim(message)
    end if

    ne = part_graph % num_edges()

    ! The whole graph is at least as large as the largest whole-graph
    ! index r records.
    largest_entry = 0
    do l = 1, part_graph % num_vertices()
       largest_entry = max(largest_entry, rel % global_vertex_index(l))
    end do
    nv_global = max(largest_entry, rel % num_whole_vertices())

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
         & stored_directed_graph(nv_global, tails=tails, heads=heads, &
         &                       number=rel % part_id()))

  end subroutine assemble_graph

  !===================================================================!
  ! Assemble the data back onto the whole graph.
  !
  ! The result is laid out on the whole graph. Only the entries this
  ! part owns are written; everything else is left at zero, so adding
  ! the results from every part rebuilds the whole field exactly once.
  !===================================================================!

  subroutine assemble_data(this, rel, part_graph, part_data, global_graph, &
       & sets, global_data)

    class(assembler) , intent(in)               :: this
    type(partition_relation), intent(in)        :: rel
    class(directed_graph)     , intent(in)               :: part_graph
    class(field), intent(in)               :: part_data
    class(directed_graph)     , intent(in)               :: global_graph
    type(set_store)    , intent(inout)            :: sets
    class(field), allocatable, intent(out) :: global_data

    type(graph) :: dom
    integer         :: n_dom
    character(len=250) :: message

    if (.not. this % defined_on_relation(rel, part_graph)) then
       write(message,'(a,i0)') 'assemble_data: rel was not written for this part; &
            &rel % part_id() = ', rel % part_id()
       error stop trim(message)
    end if

    if (.not. rel % describes_whole(global_graph)) then
       write(message,'(a,i0)') 'assemble_data: rel was not written for this whole; &
            &rel % part_id() = ', rel % part_id()
       error stop trim(message)
    end if

    select type (part_data)

    class is (stored_field)
       dom   = part_data % domain()
       n_dom = part_data % num_entries()
       ! Classify by embedding - a DECLARED predicate evaluated through
       ! the set store, never the extension and never the graph.
       if (sets % subobject_of(dom, part_graph % vertex_set())) then
          call gather_field(part_data, dom, n_dom, part_graph, rel, &
               & part_graph % vertex_set(), part_graph % num_vertices(), &
               & global_graph, .true., sets, global_data)
       else if (sets % subobject_of(dom, part_graph % edge_set())) then
          call gather_field(part_data, dom, n_dom, part_graph, rel, &
               & part_graph % edge_set(), part_graph % num_edges(), &
               & global_graph, .false., sets, global_data)
       else
          write(message,'(a)') "assemble_data: field '" // trim(part_data % name()) // &
               & "' is not defined on either of this part's domains"
          error stop trim(message)
       end if

    class default
       error stop 'assemble_data: part_data''s dynamic type is not one this transform handles'
    end select

  end subroutine assemble_data

  !===================================================================!
  ! One gather for both families and both coverages. A FULL part
  ! field maps onto the GLOBAL carrier set, owned members only,
  ! exactly the established assembly. A PROPER SUBSET maps through
  ! the part->global map onto a new subobject of the global carrier
  ! set - its mapped subdomain, no inserted zeros on members the
  ! field never stored. A new ambient set means a new declared
  ! subset: extension and values are mapped, tokens are not.
  !===================================================================!

  subroutine gather_field(part_data, dom, n_dom, part_graph, rel, part_carrier, &
       &                  n_part_carrier, global_graph, on_vertices, &
       &                  sets, global_data)

    type(stored_field)        , intent(in)               :: part_data
    type(graph)    , intent(in)               :: dom
    integer            , intent(in)               :: n_dom
    class(directed_graph)       , intent(in)               :: part_graph
    type(partition_relation), intent(in)          :: rel
    type(graph)    , intent(in)               :: part_carrier
    integer            , intent(in)               :: n_part_carrier
    class(directed_graph)       , intent(in)               :: global_graph
    logical            , intent(in)               :: on_vertices
    type(set_store)    , intent(inout)            :: sets
    class(field) , allocatable, intent(out) :: global_data

    type(stored_field)           :: out
    type(graph)       :: global_carrier
    type(graph)       :: sg
    real(dp), allocatable :: lv(:), fv(:)
    integer , allocatable :: global_members(:), origin(:)
    integer :: nglobal, nlocal, num_components, l, c, f, own_part, n, at
    character(len=250) :: message

    if (on_vertices) then
       nglobal        = global_graph % num_vertices()
       nlocal         = part_graph % num_vertices()
       global_carrier = global_graph % vertex_set()
    else
       nglobal        = global_graph % num_edges()
       nlocal         = part_graph % num_edges()
       global_carrier = global_graph % edge_set()
    end if
    num_components = part_data % num_components()
    own_part = rel % part_id()

    call part_data % real_vector(lv)
    if (size(lv) /= n_dom * num_components) then
       write(message,'(a,i0,a,i0,a,i0)') 'gather_field: the field values must fill its stated &
            &domain; size(lv) = ', size(lv), ', n_dom = ', n_dom, ', num_components = ', num_components
       error stop trim(message)
    end if

    if (dom % same_as(part_carrier)) then

       if (n_dom /= n_part_carrier) then
          write(message,'(a,i0,a,i0)') 'gather_field: a full field must fill the part carrier; &
               &n_dom = ', n_dom, ', n_part_carrier = ', n_part_carrier
          error stop trim(message)
       end if

       ! Full coverage: the established dense assembly, owned only.
       out = stored_field(part_data % name(), global_carrier, nglobal, num_components=num_components, &
            &      unit_name=part_data % units())
       allocate(fv(nglobal * num_components))
       fv = 0.0_dp

       do l = 1, nlocal
          if (rel % has_part_relation()) then
             if (rel % owner_part(l, on_vertices) /= own_part) cycle
          end if
          f = rel % global_index(l, on_vertices)
          if (f < 1 .or. f > nglobal) then
             write(message,'(a,i0,a,i0)') 'gather_field: the full-coverage relation must map into &
                  &the whole carrier 1..nglobal; f = ', f, ', nglobal = ', nglobal
             error stop trim(message)
          end if
          do c = 1, num_components
             associate (to => (f - 1) * num_components + c, from => (l - 1) * num_components + c)
               fv(to) = lv(from)
             end associate
          end do
       end do

       call out % set_real_vector(fv)

    else

       ! Proper subset: map the members to the global set and retain
       ! only the owned ones.
       if (sets % num_members_of(dom) /= n_dom) then
          write(message,'(a,i0,a,i0)') 'gather_field: a subset field must fill its stated domain; &
               &num_members_of(dom) = ', sets % num_members_of(dom), ', n_dom = ', n_dom
          error stop trim(message)
       end if
       allocate(global_members(n_dom), origin(n_dom))
       n = 0
       do l = 1, n_dom
          at = sets % member_of(dom, l)      ! part-local member
          if (at < 1 .or. at > n_part_carrier) then
             write(message,'(a,i0,a,i0)') 'gather_field: a subset must belong to the part carrier &
                  &1..n_part_carrier; at = ', at, ', n_part_carrier = ', n_part_carrier
             error stop trim(message)
          end if
          if (rel % has_part_relation()) then
             if (rel % owner_part(at, on_vertices) /= own_part) cycle
          end if
          f = rel % global_index(at, on_vertices)
          if (f < 1 .or. f > nglobal) then
             write(message,'(a,i0,a,i0)') 'gather_field: the subset relation must map into the &
                  &whole carrier 1..nglobal; f = ', f, ', nglobal = ', nglobal
             error stop trim(message)
          end if
          n = n + 1
          global_members(n) = f
          origin(n) = l
       end do
       !-------------------------------------------------------------!
       ! A new ambient set means a new declared subset, so this
       ! operation follows the subobject law: identity, extension,
       ! label and embedding, together. Extension and values are
       ! mapped to the global set; tokens are not, and the label is.
       !-------------------------------------------------------------!

       call sets % declare_subobject(sg, global_members(1:n), sets % label_of(dom), global_carrier)

       allocate(fv(n * num_components))
       do l = 1, n
          do c = 1, num_components
             fv((l - 1) * num_components + c) = lv((origin(l) - 1) * num_components + c)
          end do
       end do
       out = stored_field(part_data % name(), sg, n, num_components=num_components, &
            &      unit_name=part_data % units())
       call out % set_real_vector(fv)

    end if

    allocate(global_data, source=out)

  end subroutine gather_field

end module transform_assembler
