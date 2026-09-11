module topology_statistics
  implicit none
  type :: fibre_statistics
     integer :: num_members = 0
     integer :: member_sum = 0
     integer :: num_reads = 0
  end type fibre_statistics
contains
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
  use topology_statistics, only : fibre_statistics, sum_context_members

  implicit none

  integer :: num_failures = 0
  integer :: num_reads = 0
  integer, allocatable :: expected_members(:)
  type(stored_directed_graph) :: active_graph

  call check_fibres()
  call check_incidence()
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
