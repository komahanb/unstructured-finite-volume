program assembly_refusal

  use util_precision, only : dp
  use view_directed_stored, only : stored_directed_graph
  use view_directed, only : directed_graph
  use field_calculus, only : field
  use field_stored, only : stored_field
  use graph_fractal, only : graph
  use map_set_representation, only : counted_set_representation
  use map_set_store, only : set_store
  use transform_partitioner, only : partitioner, PARTITION_LINEAR
  use transform_assembler, only : assembler
  use relation_partition, only : partition_relation

  implicit none

  type(stored_directed_graph) :: whole, destination
  type(partitioner) :: cut
  type(assembler) :: collect
  type(partition_relation) :: rel
  type(set_store) :: sets
  type(stored_field) :: data
  type(graph) :: domain
  class(directed_graph), allocatable :: part
  class(field), allocatable :: part_data, global_data
  character(len=32) :: mode

  call get_command_argument(1, mode)

  whole = stored_directed_graph(6, tails=[1, 2, 3, 4, 5], heads=[2, 3, 4, 5, 6])
  destination = whole
  domain = whole % vertex_set()
  call sets % bind(domain, counted_set_representation(6))
  data = stored_field('q', domain, 6)
  call data % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp])
  cut = partitioner(PARTITION_LINEAR, num_parts=1, part=1)
  collect = assembler()
  call cut % partition_graph(whole, part, rel)
  call sets % bind(part % vertex_set(), counted_set_representation(part % num_vertices()))
  call sets % bind(part % edge_set(), counted_set_representation(part % num_edges()))
  call cut % partition_data(rel, whole, data, part, sets, part_data)

  select case (trim(mode))
  case ('wrong-whole')
     destination = stored_directed_graph(6, tails=[1, 1, 1, 1, 1], heads=[2, 3, 4, 5, 6])
  case ('smaller-whole')
     destination = stored_directed_graph(3, tails=[1, 2], heads=[2, 3])
  case ('whole-vertex-count')
     destination % nv = 3
  case ('whole-edge-count')
     destination % ne = 4
  case ('part-count')
     select type (part)
     type is (stored_directed_graph)
        part % nv = 5
     end select
  case ('out-of-range-full', 'out-of-range-subset')
     rel = partition_relation(part % vertex_set(), 6, part % edge_set(), 5, &
          & whole % vertex_set(), 6, whole % edge_set(), 5, 1, 1, &
          & [1, 2, 3, 4, 5, 7], [1, 1, 1, 1, 1, 1], &
          & [1, 2, 3, 4, 5], [1, 1, 1, 1, 1])
     if (trim(mode) == 'out-of-range-subset') then
        call sets % declare_subobject(domain, [6], 'last member', part % vertex_set())
        data = stored_field('q subset', domain, 1)
        call data % set_real_vector([60.0_dp])
        deallocate(part_data)
        allocate(part_data, source=data)
     end if
  case ('field-count')
     data = stored_field('q short', part % vertex_set(), 5)
     call data % set_real_vector([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp])
     deallocate(part_data)
     allocate(part_data, source=data)
  case default
     error stop 'assembly_refusal: unknown case'
  end select

  call collect % assemble_data(rel, part, part_data, destination, sets, global_data)
  error stop 'assembly_refusal: invalid assembly was accepted'

end program assembly_refusal
