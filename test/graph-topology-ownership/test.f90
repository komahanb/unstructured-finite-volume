module topology_statistics
  use util_counted_storage, only : counted_storage, counted_reference
  use view_directed_stored, only : stored_directed_graph
  use view_level, only : level_storage
  implicit none
  ! A containing object: every copy of the container copies its graph.
  type :: graph_holder
     type(stored_directed_graph) :: incidence
  end type graph_holder
  ! A cell of the counted primitive with a payload, and the number of
  ! times any cell was cleared.
  type, extends(counted_storage) :: counted_test_cell
     integer, allocatable :: payload(:)
   contains
     procedure :: clear => clear_test_cell
  end type counted_test_cell
  type :: reference_holder
     type(counted_reference) :: reference
  end type reference_holder
  ! A containing object of a hierarchy: every copy of the container
  ! binds one more owner of the same nodes.
  type :: hierarchy_holder
     type(level_storage) :: nodes
  end type hierarchy_holder
  integer :: num_cleared = 0
  type :: fibre_statistics
     integer :: num_members = 0
     integer :: member_sum = 0
     integer :: num_reads = 0
  end type fibre_statistics
contains
  subroutine clear_test_cell(this)
    class(counted_test_cell), intent(inout) :: this
    if (allocated(this % payload)) deallocate(this % payload)
    !$omp atomic
    num_cleared = num_cleared + 1
  end subroutine clear_test_cell

  subroutine sum_context_members(members, context)
    integer, intent(in) :: members(:)
    class(*), intent(inout) :: context
    select type(context)
    type is(fibre_statistics)
       context % num_members = context % num_members + size(members)
       context % member_sum = context % member_sum + sum(members)
       context % num_reads = context % num_reads + 1
    class default
       error stop 'the fibre statistics context is required'
    end select
  end subroutine sum_context_members
end module topology_statistics

program topology_ownership

  use graph_fractal, only : graph
  use map_set, only : set_map
  use map_set_representation, only : listed_set_representation
  use relation_binary, only : csr_relation, transposed_relation, transpose_of, integer_fibre
  use view_directed_stored, only : stored_directed_graph
  use relation_partition, only : partition_relation
  use util_counted_storage, only : counted_storage, counted_reference
  use view_level, only : level_storage, level_num_members, level_member, level_is_leaf
  use topology_statistics, only : fibre_statistics, sum_context_members, graph_holder, &
       & counted_test_cell, reference_holder, num_cleared, hierarchy_holder
  !$ use omp_lib, only : omp_get_num_threads

  implicit none

  integer :: num_failures = 0
  integer :: num_reads = 0
  integer, allocatable :: expected_members(:)
  type(stored_directed_graph) :: active_graph

  call check_fibres()
  call check_incidence()
  call check_orientation()
  call check_counted_storage()
  call check_registry_consistency()
  call check_hierarchy_ownership()
  if (num_failures /= 0) error stop 'topology ownership laws failed'
  print *, 'All topology ownership laws passed.'

contains

  subroutine assert_all(satisfied, description)
    logical, intent(in) :: satisfied
    character(len=*), intent(in) :: description
    if (satisfied) then
       print *, 'PASS : ', description
    else
       print *, 'FAIL : ', description
       num_failures = num_failures + 1
    end if
  end subroutine assert_all

  subroutine check_fibres()
    type(graph) :: source_set, target_set, domain
    type(set_map) :: sets
    type(csr_relation), target :: relation
    type(transposed_relation), target :: transposed_relation_value
    type(transposed_relation) :: twice_transposed
    type(integer_fibre) :: fibre, default_fibre, second_fibre
    type(fibre_statistics) :: image_statistics, preimage_statistics
    integer, parameter :: source_members(4) = [30, 10, 70, 50]
    integer, parameter :: target_members(4) = [900, 100, 400, 200]
    integer, allocatable :: owned(:), values(:), tuples(:,:)
    integer :: i, j, source_member, target_member

    call source_set % declare()
    call target_set % declare()
    call sets % bind(source_set, listed_set_representation(source_members))
    call sets % bind(target_set, listed_set_representation(target_members))
    relation = csr_relation('sparse incidence', source_set, target_set, &
         & reshape([30,900, 10,400, 30,100, 70,400, 30,900], [2,5]), sets)
    transposed_relation_value = transpose_of(relation)
    twice_transposed = transpose_of(transposed_relation_value)
    domain = transposed_relation_value % source()
    call assert_all(domain % same_as(target_set), 'transpose preserves the target domain identity')
    domain = transposed_relation_value % target()
    call assert_all(domain % same_as(source_set), 'transpose preserves the source domain identity')
    call assert_all(relation % num_tuples() == 4, 'duplicate incidence has one tuple')

    fibre = relation % image_view(30)
    second_fibre = relation % preimage_view(400)
    call fibre % read(sum_context_members, image_statistics)
    call second_fibre % read(sum_context_members, preimage_statistics)
    call fibre % read(sum_context_members, image_statistics)
    call default_fibre % read(sum_context_members, preimage_statistics)
    call assert_all(image_statistics % num_members == 4 .and. image_statistics % member_sum == 2000 .and. &
         & image_statistics % num_reads == 2, 'interleaved image reads update only their explicit context')
    call assert_all(preimage_statistics % num_members == 2 .and. preimage_statistics % member_sum == 80 .and. &
         & preimage_statistics % num_reads == 2, 'empty and nonempty preimage reads preserve a separate context')

    do i = 1, size(source_members)
       source_member = source_members(i)
       fibre = relation % image_view(source_member)
       call relation % image(source_member, owned)
       values = fibre % values()
       call assert_all(fibre % num_members() == size(owned) .and. all(values == owned), &
            & 'sparse image agrees with its owning copy')
       expected_members = owned
       num_reads = 0
       call fibre % read(check_members)
       call assert_all(num_reads == 1, 'fibre reader executes exactly once, including empty fibres')
       second_fibre = transposed_relation_value % preimage_view(source_member)
       values = second_fibre % values()
       call assert_all(all(values == owned), 'image equals transpose preimage')
       do j = 1, fibre % num_members()
          target_member = fibre % member(j)
          call assert_all(relation % has([source_member, target_member]), &
               & 'fibre member belongs to the relation')
       end do
       if (size(owned) > 0) then
          owned = -1
          call assert_all(all(fibre % values() /= -1), 'mutating an owned image preserves its source')
          values = fibre % values()
          values = -2
          call assert_all(all(fibre % values() /= -2), 'mutating values preserves its source')
       end if
    end do

    do i = 1, size(target_members)
       target_member = target_members(i)
       fibre = relation % preimage_view(target_member)
       call relation % preimage(target_member, owned)
       second_fibre = transposed_relation_value % image_view(target_member)
       values = second_fibre % values()
       call assert_all(fibre % num_members() == size(owned) .and. all(values == owned), &
            & 'sparse preimage equals transpose image')
       do j = 1, size(source_members)
          source_member = source_members(j)
          call assert_all(relation % has([source_member, target_member]) .eqv. &
               & transposed_relation_value % has([target_member, source_member]), &
               & 'incidence membership is equivalent in both orientations')
       end do
    end do

    fibre = relation % image_view(999)
    call assert_all(fibre % num_members() == 0, 'a source outsider has empty image')
    fibre = relation % preimage_view(999)
    call assert_all(fibre % num_members() == 0, 'a target outsider has empty preimage')
    expected_members = [integer ::]
    num_reads = 0
    call default_fibre % read(check_members)
    call assert_all(default_fibre % num_members() == 0 .and. num_reads == 1, &
         & 'default fibre is empty and readable')
    values = default_fibre % values()
    call assert_all(size(values) == 0, 'default fibre has an empty owning copy')
    call twice_transposed % tuples(tuples)
    call assert_all(size(tuples, 2) == relation % num_tuples(), 'transpose involution preserves cardinality')
    do i = 1, size(tuples, 2)
       call assert_all(relation % has(tuples(:, i)), 'transpose involution preserves each tuple')
    end do
  end subroutine check_fibres

  subroutine check_members(members)
    integer, intent(in) :: members(:)
    num_reads = num_reads + 1
    call assert_all(size(members) == size(expected_members), 'fibre reader preserves cardinality')
    if (size(members) == size(expected_members)) then
       call assert_all(all(members == expected_members), 'fibre reader preserves member order')
    end if
  end subroutine check_members

  subroutine check_incidence()
    type(stored_directed_graph) :: original, transposed_graph, twice_transposed, default_graph
    type(graph) :: original_domain, transposed_domain
    integer, allocatable :: incidence(:), outgoing(:), incoming(:)
    integer :: vertex, edge_index, num_incidences

    active_graph = default_graph
    num_reads = 0
    call active_graph % read_incoming(check_incoming)
    call assert_all(num_reads == 1, 'default empty graph invokes one reader')
    active_graph = stored_directed_graph(0, [integer ::], [integer ::])
    num_reads = 0
    call active_graph % read_incoming(check_incoming)
    call assert_all(num_reads == 1, 'constructed empty graph invokes one reader')
    original = stored_directed_graph(1, [1,1,1], [-1,0,2])
    call assert_all(all([(original % edge_head(edge_index) == 0, edge_index=1,3)]), &
         & 'all heads outside the vertex carrier denote a boundary')

    original = stored_directed_graph(5, [1,1,2,3,3], [2,3,4,4,0], number=23)
    num_incidences = 0
    do vertex = 1, original % num_vertices()
       call original % incident_edges(vertex, incidence)
       num_incidences = num_incidences + size(incidence)
       do edge_index = 1, original % num_edges()
          call assert_all(any(incidence == edge_index) .eqv. &
               & (original % edge_tail(edge_index) == vertex .or. original % edge_head(edge_index) == vertex), &
               & 'incident membership equals endpoint membership')
       end do
    end do
    call assert_all(num_incidences == 9, 'four interior edges and one boundary have nine incidences')
    call assert_all(.not. original % edge_has_head(5) .and. original % edge_head(5) == 0, &
         & 'the boundary has no artificial opposite vertex')
    active_graph = original
    num_reads = 0
    call active_graph % read_incoming(check_incoming)
    call assert_all(num_reads == 1, 'graph traversal invokes one reader')

    original = stored_directed_graph(5, [1,1,2,3,4], [2,3,4,4,4], number=29)
    transposed_graph = original % transpose()
    twice_transposed = transposed_graph % transpose()
    original_domain = original % vertex_set()
    transposed_domain = transposed_graph % vertex_set()
    call assert_all(original_domain % same_as(transposed_domain), 'directed transpose preserves vertex identity')
    original_domain = original % edge_set()
    transposed_domain = transposed_graph % edge_set()
    call assert_all(original_domain % same_as(transposed_domain), 'directed transpose preserves edge identity')
    call assert_all(original % id() == twice_transposed % id(), 'directed transpose involution preserves identity')
    do edge_index = 1, original % num_edges()
       call assert_all(original % edge_tail(edge_index) == transposed_graph % edge_head(edge_index) .and. &
            & original % edge_head(edge_index) == transposed_graph % edge_tail(edge_index), &
            & 'directed transpose exchanges endpoints')
       call assert_all(original % edge_tail(edge_index) == twice_transposed % edge_tail(edge_index) .and. &
            & original % edge_head(edge_index) == twice_transposed % edge_head(edge_index), &
            & 'directed transpose involution restores endpoints')
    end do
    do vertex = 1, original % num_vertices()
       call original % outgoing_edges(vertex, outgoing)
       call transposed_graph % incoming_edges(vertex, incoming)
       call assert_all(size(outgoing) == size(incoming) .and. all(outgoing == incoming), &
            & 'transpose incoming equals original outgoing')
    end do
    active_graph = transposed_graph
    num_reads = 0
    call active_graph % read_incoming(check_incoming)
    call assert_all(num_reads == 1, 'transposed graph invokes one reader')
  end subroutine check_incidence

  !-------------------------------------------------------------------!
  ! Orientation laws of a stored graph value. The owning transpose is
  ! valid after its source is finalized or replaced; reverse changes
  ! the orientation of one value in place; every copy of a value or of
  ! an object containing it is independent under every copy mechanism.
  !-------------------------------------------------------------------!
  subroutine check_orientation()
    integer, parameter :: tails(5) = [1,1,2,3,4], heads(5) = [2,3,4,4,4]
    integer, parameter :: global_vertices(5) = [5,4,3,2,1]
    type(stored_directed_graph), allocatable :: original, transposed_graph, duplicate, reversed_graph
    type(stored_directed_graph) :: replaced, replacement_transpose
    type(graph_holder), allocatable :: holder, second_holder, holders(:)
    type(graph) :: vertex_identity, edge_identity, transposed_identity
    type(partition_relation) :: relation
    integer, allocatable :: outgoing(:), incoming(:)
    integer :: edge_index, vertex, i

    allocate(original, transposed_graph)
    original = stored_directed_graph(5, tails, heads, number=29, num_parts=30, vglobal=global_vertices)
    vertex_identity = original % vertex_set()
    edge_identity = original % edge_set()
    transposed_graph = original % transpose()
    deallocate(original)
    call assert_all(transposed_graph % transposed() .and. transposed_graph % id() == 29, &
         & 'the transpose keeps orientation and identity after its source is finalized')
    do edge_index = 1, 5
       call assert_all(transposed_graph % edge_tail(edge_index) == heads(edge_index) .and. &
            & transposed_graph % edge_head(edge_index) == tails(edge_index), &
            & 'the transpose reads exchanged endpoints after its source is finalized')
    end do
    transposed_identity = transposed_graph % vertex_set()
    call assert_all(transposed_identity % same_as(vertex_identity), 'the transpose keeps the vertex identity')
    transposed_identity = transposed_graph % edge_set()
    call assert_all(transposed_identity % same_as(edge_identity), 'the transpose keeps the edge identity')
    relation = transposed_graph % whole_relation()
    call assert_all(relation % has_part_relation() .and. &
         & all([(relation % global_vertex_index(vertex) == global_vertices(vertex), vertex = 1, 5)]), &
         & 'the transpose keeps the partition maps')
    active_graph = transposed_graph
    num_reads = 0
    call active_graph % read_incoming(check_incoming)
    call assert_all(num_reads == 1, 'the transpose reads its compressed incidence after its source is finalized')

    ! reverse in place equals the transpose and is an involution
    allocate(reversed_graph)
    reversed_graph = stored_directed_graph(5, tails, heads, number=29)
    call reversed_graph % reverse()
    call assert_all(reversed_graph % transposed(), 'reverse changes the orientation in place')
    do edge_index = 1, 5
       call assert_all(reversed_graph % edge_tail(edge_index) == transposed_graph % edge_tail(edge_index) .and. &
            & reversed_graph % edge_head(edge_index) == transposed_graph % edge_head(edge_index), &
            & 'the reversed value reads the transpose endpoints')
    end do
    do vertex = 1, 5
       call reversed_graph % incoming_edges(vertex, incoming)
       call transposed_graph % incoming_edges(vertex, outgoing)
       call assert_all(size(outgoing) == size(incoming) .and. all(outgoing == incoming), &
            & 'the reversed value reads the transpose incidence')
    end do
    active_graph = reversed_graph
    num_reads = 0
    call active_graph % read_incoming(check_incoming)
    call assert_all(num_reads == 1, 'the reversed value invokes one reader')
    call reversed_graph % reverse()
    call assert_all(.not. reversed_graph % transposed() .and. &
         & all([(reversed_graph % edge_tail(edge_index) == tails(edge_index) .and. &
         &        reversed_graph % edge_head(edge_index) == heads(edge_index), edge_index = 1, 5)]), &
         & 'reversing twice restores the orientation')

    ! independence of copies under every copy mechanism
    allocate(holder, second_holder, holders(3))
    holder % incidence = reversed_graph
    second_holder = holder
    call second_holder % incidence % reverse()
    call assert_all(.not. holder % incidence % transposed() .and. second_holder % incidence % transposed(), &
         & 'a container copy has an independent orientation')
    do i = 1, 3
       holders(i) % incidence = holder % incidence
    end do
    call holders(2) % incidence % reverse()
    call assert_all(.not. holders(1) % incidence % transposed() .and. holders(2) % incidence % transposed() .and. &
         & .not. holders(3) % incidence % transposed() .and. .not. holder % incidence % transposed(), &
         & 'array elements of containers have independent orientations')
    allocate(duplicate, source=reversed_graph)
    call duplicate % reverse()
    call assert_all(duplicate % transposed() .and. .not. reversed_graph % transposed() .and. &
         & duplicate % edge_tail(2) == 3 .and. reversed_graph % edge_tail(2) == 1, &
         & 'a source= copy has independent orientation and reads its own arrays')
    deallocate(reversed_graph)
    call assert_all(duplicate % edge_head(5) == 4 .and. duplicate % num_edges() == 5, &
         & 'a source= copy is valid after its source is finalized')
    block
      type(stored_directed_graph) :: local_copy
      local_copy = transpose_result(duplicate)
      call assert_all(.not. local_copy % transposed() .and. duplicate % transposed(), &
           & 'a function result is an independent value')
    end block
    call assert_all(duplicate % transposed() .and. duplicate % edge_tail(2) == 3, &
         & 'block scope end leaves other values unchanged')
    deallocate(holder, second_holder, holders)

    ! replacement: the transpose survives the replacement of its source
    replaced = stored_directed_graph(3, [1,2], [2,3], number=7)
    replacement_transpose = replaced % transpose()
    replaced = stored_directed_graph(2, [1], [2], number=8)
    call assert_all(replacement_transpose % num_vertices() == 3 .and. replacement_transpose % edge_tail(2) == 3 .and. &
         & replacement_transpose % id() == 7, 'the transpose reads its own arrays after its source is replaced')
    call assert_all(replaced % num_vertices() == 2 .and. replaced % edge_head(1) == 2 .and. replaced % id() == 8, &
         & 'the replaced value reads the new arrays')
  end subroutine check_orientation

  type(stored_directed_graph) function transpose_result(source)
    type(stored_directed_graph), intent(in) :: source
    transpose_result = source % transpose()
  end function transpose_result

  !-------------------------------------------------------------------!
  ! The counted ownership primitive: owners bound through defined
  ! assignment, released by finalization, the cell cleared with the
  ! last owner and reused; a source= copy is not an owner and observes
  ! the count. The mechanisms measured in doc/topology-ownership.md.
  !-------------------------------------------------------------------!
  subroutine check_counted_storage()
    type(counted_test_cell) :: template
    type(counted_reference), allocatable :: first, second, observer, moved
    type(counted_reference) :: replacement
    type(reference_holder), allocatable :: holder, second_holder, holders(:)
    class(counted_storage), pointer :: first_cell, second_cell
    integer :: cleared_before, i

    cleared_before = num_cleared
    allocate(first, second)
    call assert_all(.not. first % live() .and. first % num_owners() == 0 .and. .not. associated(first % storage()), &
         & 'a default reference is not live and has no storage')
    call first % acquire(template)
    first_cell => first % storage()
    select type (first_cell)
    type is (counted_test_cell)
       allocate(first_cell % payload(4))
       first_cell % payload = 7
    class default
       call assert_all(.false., 'the acquired cell has the template type')
    end select
    call assert_all(first % live() .and. first % num_owners() == 1, 'acquisition binds one owner')
    second = first
    call assert_all(second % live() .and. second % num_owners() == 2 .and. associated(second % storage(), first_cell), &
         & 'assignment binds a second owner of the same cell')
    allocate(observer, source=first)
    call assert_all(observer % num_owners() == 2, 'a source= copy observes the count and is not an owner')
    first = first
    call assert_all(first % num_owners() == 2, 'self assignment changes no owner')
    second = first
    call assert_all(second % num_owners() == 2, 'assignment between owners of one cell changes no owner')
    deallocate(first)
    call assert_all(second % num_owners() == 1 .and. num_cleared == cleared_before, &
         & 'finalizing one owner releases one binding and clears nothing')

    allocate(holder, second_holder, holders(3))
    holder % reference = second
    second_holder = holder
    call assert_all(second % num_owners() == 3, 'container assignment binds the contained reference')
    do i = 1, 3
       holders(i) % reference = second_holder % reference
    end do
    call assert_all(second % num_owners() == 6, 'array elements of containers are owners')
    deallocate(holder, second_holder)
    call assert_all(second % num_owners() == 4, 'finalizing containers releases their bindings')
    deallocate(holders)
    call assert_all(second % num_owners() == 1, 'finalizing an array of containers releases every element')
    block
      type(counted_reference) :: local_reference
      type(reference_holder) :: local_holder
      local_reference = second
      local_holder % reference = reference_result(second)
      call assert_all(second % num_owners() == 3, 'block locals and a function result are owners')
    end block
    call assert_all(second % num_owners() == 1, 'block scope end releases its owners')
    call move_alloc(second, moved)
    call assert_all(moved % num_owners() == 1, 'move_alloc transfers a binding without a new owner')
    replacement = moved
    call replacement % acquire(template)
    call assert_all(moved % num_owners() == 1 .and. replacement % num_owners() == 1 .and. &
         & .not. associated(replacement % storage(), first_cell), 'acquiring over a binding releases it')
    deallocate(moved)
    call assert_all(num_cleared == cleared_before + 1 .and. observer % num_owners() == 0 .and. .not. observer % live() &
         & .and. .not. associated(observer % storage()), &
         & 'the last owner clears the cell and the observer is not live')
    allocate(second)
    call second % acquire(template)
    second_cell => second % storage()
    call assert_all(associated(second_cell, first_cell) .and. .not. observer % live(), &
         & 'a cleared cell is reused with a new version and the stale observer stays not live')
    deallocate(observer)
    call assert_all(second % num_owners() == 1 .and. replacement % num_owners() == 1, &
         & 'finalizing a stale observer releases nothing')
  end subroutine check_counted_storage

  !-------------------------------------------------------------------!
  ! The registry: n sole owners released clear n cells; n acquisitions
  ! from that free list address n distinct cells with one owner each;
  ! n bindings on one cell made and released return its owner count
  ! to one, and n assignments over sole owners clear their n cells.
  ! The loops run on two threads when built with -fopenmp and are the
  ! serial loops otherwise; the thread count observed is printed.
  !-------------------------------------------------------------------!
  subroutine check_registry_consistency()
    integer, parameter :: num_cells = 256
    type(counted_test_cell) :: template
    type(counted_reference), allocatable :: references(:)
    type(counted_reference) :: shared
    class(counted_storage), pointer :: cell_i, cell_j
    integer :: cleared_before, i, j, num_threads
    logical :: distinct, sole_owners, bound

    cleared_before = num_cleared
    allocate(references(num_cells))
    do i = 1, num_cells
       call references(i) % acquire(template)
    end do
    deallocate(references)
    call assert_all(num_cleared == cleared_before + num_cells, 'releasing n sole owners clears n cells')

    allocate(references(num_cells))
    num_threads = 1
    !$omp parallel num_threads(2)
    !$omp single
    !$ num_threads = omp_get_num_threads()
    !$omp end single
    !$omp do
    do i = 1, num_cells
       call references(i) % acquire(template)
    end do
    !$omp end do
    !$omp end parallel
    print '(1x,a,i0,a)', 'registry operations over ', num_threads, ' thread(s)'
    distinct = .true.
    sole_owners = .true.
    do i = 1, num_cells
       sole_owners = sole_owners .and. references(i) % live() .and. references(i) % num_owners() == 1
       cell_i => references(i) % storage()
       do j = i + 1, num_cells
          cell_j => references(j) % storage()
          if (associated(cell_i, cell_j)) distinct = .false.
       end do
    end do
    call assert_all(distinct, 'n acquisitions from a free list of n cells address n distinct cells')
    call assert_all(sole_owners .and. num_cleared == cleared_before + num_cells, &
         & 'each acquired cell has one owner and none was cleared')

    call shared % acquire(template)
    bound = .true.
    !$omp parallel do num_threads(2) reduction(.and.:bound)
    do i = 1, num_cells
       bound = bound .and. bound_and_released(shared)
    end do
    !$omp end parallel do
    call assert_all(bound .and. shared % num_owners() == 1 .and. num_cleared == cleared_before + num_cells, &
         & 'n bindings made and released on one cell return its owner count to one')

    !$omp parallel do num_threads(2)
    do i = 1, num_cells
       references(i) = shared
    end do
    !$omp end parallel do
    call assert_all(shared % num_owners() == num_cells + 1 .and. num_cleared == cleared_before + 2 * num_cells, &
         & 'n assignments over sole owners bind n more owners and clear n cells')
    deallocate(references)
    call assert_all(shared % num_owners() == 1 .and. num_cleared == cleared_before + 2 * num_cells, &
         & 'finalizing the n owners leaves the original owner and clears nothing')
  end subroutine check_registry_consistency

  ! Bind a local reference to the shared cell and release it on return.
  logical function bound_and_released(shared)
    type(counted_reference), intent(in) :: shared
    type(counted_reference) :: local_reference
    local_reference = shared
    bound_and_released = local_reference % live() .and. local_reference % num_owners() >= 2 .and. &
         & associated(local_reference % storage(), shared % storage())
  end function bound_and_released

  type(counted_reference) function reference_result(source)
    type(counted_reference), intent(in) :: source
    reference_result = source
  end function reference_result

  !-------------------------------------------------------------------!
  ! A hierarchy (level_storage) is jointly owned immutable topology:
  ! every copy by defined assignment reads the same nodes with the
  ! same identities; the nodes are released with the last owner; a
  ! hierarchy with more than one owner is not extended; a source=
  ! copy observes the count and is not an owner.
  !-------------------------------------------------------------------!
  subroutine check_hierarchy_ownership()
    type(level_storage), allocatable :: original, copy, twin
    type(level_storage) :: elements(2), replacement
    type(hierarchy_holder) :: holder, holder_copy
    type(graph), pointer :: root, member, copied_root
    integer :: leaf_a, leaf_b, root_at

    allocate(original, copy)
    call assert_all(original % num_nodes() == 0 .and. original % num_owners() == 0, &
         & 'a hierarchy without nodes has no owner')
    leaf_a = original % allocate_node()
    leaf_b = original % allocate_node()
    root_at = original % assemble([leaf_a, leaf_b], 0)
    root => original % node(root_at)
    call assert_all(original % num_nodes() == 5 .and. original % num_owners() == 1 .and. &
         & level_num_members(root) == 2, 'a level of two leaves owns five nodes with one owner')
    copy = original
    copied_root => copy % node(root_at)
    member => level_member(copied_root, 2)
    call assert_all(copy % num_owners() == 2 .and. associated(copied_root, root) .and. &
         & associated(member, original % node(leaf_b)) .and. member % same_as(original % node(leaf_b)), &
         & 'a copy reads the same nodes with the same identities')
    holder % nodes = copy
    holder_copy = holder
    elements(1) = copy
    elements(2) = hierarchy_result(copy)
    call assert_all(original % num_owners() == 6, &
         & 'a container, its copy, an array element and a function result are owners')
    allocate(twin, source=copy)
    call assert_all(twin % num_owners() == 6 .and. twin % num_nodes() == 5, &
         & 'a source= copy observes the count and the nodes and is not an owner')
    deallocate(twin)
    call assert_all(original % num_owners() == 5 .and. original % num_nodes() == 5, &
         & 'finalizing the twin releases one binding, its owner twin, and the others read on')
    deallocate(original)
    root => copy % node(root_at)
    member => level_member(root, 1)
    call assert_all(copy % num_owners() == 4 .and. copy % num_nodes() == 5 .and. &
         & level_is_leaf(member) .and. member % same_as(holder % nodes % node(leaf_a)), &
         & 'the source destroyed first leaves the destination reading the nodes')
    block
      type(level_storage) :: local
      local = copy
      call assert_all(copy % num_owners() == 5, 'a block local is an owner')
    end block
    call assert_all(copy % num_owners() == 4, 'block scope end releases its owner')
    replacement = holder % nodes
    holder % nodes = elements(1)
    call assert_all(copy % num_owners() == 5 .and. replacement % num_owners() == 5, &
         & 'assignment between owners of one hierarchy keeps every owner')
    deallocate(copy)
    elements(1) = level_storage_result()
    elements(2) = elements(1)
    holder_copy = hierarchy_holder()
    call assert_all(holder % nodes % num_owners() == 2 .and. elements(1) % num_owners() == 0 .and. &
         & elements(1) % num_nodes() == 0, &
         & 'the destination destroyed first leaves the source and its remaining owners')
    root => replacement % node(root_at)
    call assert_all(level_num_members(root) == 2 .and. replacement % num_nodes() == 5, &
         & 'the last owners still read the hierarchy')
    holder % nodes = elements(1)
    call assert_all(replacement % num_owners() == 1, 'one owner remains')
    leaf_a = replacement % allocate_node()
    call assert_all(replacement % num_nodes() == 6 .and. replacement % num_owners() == 1, &
         & 'the sole owner extends its hierarchy')
    call assert_all(elements(1) % allocate_node() == 1 .and. elements(1) % num_owners() == 1, &
         & 'a released storage builds a new hierarchy')
  end subroutine check_hierarchy_ownership

  type(level_storage) function hierarchy_result(source)
    type(level_storage), intent(in) :: source
    hierarchy_result = source
  end function hierarchy_result

  type(level_storage) function level_storage_result()
    type(level_storage) :: unbound
    level_storage_result = unbound
  end function level_storage_result

  subroutine check_incoming(offsets, indices, sources)
    integer, intent(in) :: offsets(:), indices(:), sources(:)
    integer, allocatable :: incoming(:)
    integer :: vertex, position, edge_index

    num_reads = num_reads + 1
    call assert_all(size(offsets) == active_graph % num_vertices() + 1, 'incoming offsets cover all vertices')
    call assert_all(size(sources) == active_graph % num_edges(), 'source endpoints cover all edges')
    do vertex = 1, active_graph % num_vertices()
       call active_graph % incoming_edges(vertex, incoming)
       call assert_all(all(indices(offsets(vertex):offsets(vertex+1)-1) == incoming), &
            & 'read-only compressed incoming traversal equals the public incidence query')
       do position = offsets(vertex), offsets(vertex+1) - 1
          edge_index = indices(position)
          call assert_all(sources(edge_index) == active_graph % edge_tail(edge_index) .and. &
               & vertex == active_graph % edge_head(edge_index), 'compressed traversal preserves edge orientation')
       end do
    end do
  end subroutine check_incoming

end program topology_ownership
