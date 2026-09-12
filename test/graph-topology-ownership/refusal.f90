program topology_refusal
  use graph_fractal, only : graph
  use view_directed_stored, only : stored_directed_graph
  use relation_binary, only : integer_fibre
  use view_level, only : level_storage
  implicit none
  type(stored_directed_graph) :: incidence
  type(integer_fibre) :: fibre
  type(graph) :: undeclared_domain
  type(graph), pointer :: node
  type(level_storage) :: hierarchy, shared
  type(level_storage), allocatable :: twin
  integer :: at
  integer, target :: members(2) = [3,7]
  character(len=32) :: case_name

  call get_command_argument(1, case_name)
  select case(trim(case_name))
  case('negative_count')
     incidence = stored_directed_graph(-1, [integer ::], [integer ::])
  case('endpoint_extents')
     incidence = stored_directed_graph(2, [1,2], [2])
  case('zero_tail')
     incidence = stored_directed_graph(2, [0], [2])
  case('outside_tail')
     incidence = stored_directed_graph(2, [3], [2])
  case('vertex_tags')
     incidence = stored_directed_graph(2, [1], [2], vtags=['a'])
  case('edge_tags')
     incidence = stored_directed_graph(2, [1], [2], etags=['a','b'])
  case('global_vertex_extent')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1])
  case('global_vertex_range')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,3])
  case('global_vertex_duplicates')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,1])
  case('global_edge_extent')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], eglobal=[1,2])
  case('global_edge_range')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], eglobal=[2])
  case('global_edge_duplicates')
     incidence = stored_directed_graph(2, [1,2], [2,1], vglobal=[1,2], eglobal=[1,1])
  case('vertex_owner_extent')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], vowner=[1])
  case('vertex_owner_range')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], vowner=[1,2])
  case('edge_owner_extent')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], eowner=[1,1])
  case('edge_owner_range')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], eowner=[0])
  case('missing_partition_vertices')
     incidence = stored_directed_graph(2, [1], [2], eglobal=[1])
  case('partition_number')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], number=0)
  case('undeclared_whole_vertices')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], whole_vertices=undeclared_domain)
  case('undeclared_whole_edges')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], whole_edges=undeclared_domain)
  case('inconsistent_vertex_extent')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], num_whole_vertices=3)
  case('inconsistent_edge_extent')
     incidence = stored_directed_graph(2, [1], [2], vglobal=[1,2], num_whole_edges=2)
  case('boundary_transpose')
     incidence = stored_directed_graph(2, [1], [0])
     incidence = incidence % transpose()
  case('boundary_reverse')
     incidence = stored_directed_graph(2, [1], [0])
     call incidence % reverse()
  case('zero_fibre_index')
     fibre = integer_fibre(members)
     print *, fibre % member(0)
  case('outside_fibre_index')
     fibre = integer_fibre(members)
     print *, fibre % member(3)
  case('empty_fibre_index')
     print *, fibre % member(1)
  case('hierarchy_shared_extension')
     at = hierarchy % allocate_node()
     shared = hierarchy
     at = hierarchy % allocate_node()
  case('hierarchy_released_twin')
     at = hierarchy % allocate_node()
     allocate(twin, source=hierarchy)
     deallocate(twin)
     node => hierarchy % node(at)
  case('hierarchy_empty_index')
     node => hierarchy % node(1)
  case default
     error stop 'unknown refusal case'
  end select
  print *, 'Invalid topology input was admitted.'
end program topology_refusal
