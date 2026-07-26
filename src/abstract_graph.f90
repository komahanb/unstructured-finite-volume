!=====================================================================!
! Abstract graph architecture interfaces.
!
! This module contains only abstract types and deferred type-bound
! procedures. Concrete graph storage, partitioning algorithms, field
! layouts, reductions, transforms, and operations must extend these
! contracts in separate class modules.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module abstract_graph_types

  use iso_fortran_env, only : dp => REAL64

  implicit none

  private

  public :: graph
  public :: graph_entity
  public :: graph_vertex
  public :: graph_edge
  public :: graph_support
  public :: graph_vertex_support
  public :: graph_edge_support

  public :: graph_data
  public :: graph_field
  public :: graph_vertex_field
  public :: graph_edge_field
  public :: graph_functional
  public :: graph_reduction

  public :: graph_transform
  public :: graph_partition_transform
  public :: graph_partitioner
  public :: graph_assembler
  public :: graph_resolution_transform
  public :: graph_coarsener
  public :: graph_refiner

  public :: graph_operation
  public :: graph_vertex_operation
  public :: graph_edge_operation
  public :: graph_vertex_field_operation
  public :: graph_vertex_functional_operation
  public :: graph_edge_field_operation
  public :: graph_edge_functional_operation

  !===================================================================!
  ! Structure interfaces
  !===================================================================!

  type, abstract :: graph_entity
   contains
     procedure(entity_id_interface)   , deferred :: id
     procedure(entity_label_interface), deferred :: label
  end type graph_entity

  type, abstract, extends(graph_entity) :: graph_vertex
  end type graph_vertex

  type, abstract, extends(graph_entity) :: graph_edge
   contains
     procedure(edge_tail_interface), deferred :: tail
     procedure(edge_head_interface), deferred :: head
     procedure(edge_has_head_interface), deferred :: has_head
  end type graph_edge

  type, abstract :: graph_support
   contains
     procedure(support_kind_interface)    , deferred :: kind
     procedure(support_size_interface)    , deferred :: size
     procedure(support_ids_interface)     , deferred :: ids
     procedure(support_contains_interface), deferred :: contains_id
  end type graph_support

  type, abstract, extends(graph_support) :: graph_vertex_support
  end type graph_vertex_support

  type, abstract, extends(graph_support) :: graph_edge_support
  end type graph_edge_support

  type, abstract :: graph
   contains
     procedure(graph_id_interface)               , deferred :: id
     procedure(graph_num_vertices_interface)     , deferred :: num_vertices
     procedure(graph_num_edges_interface)        , deferred :: num_edges
     procedure(graph_vertex_interface)           , deferred :: vertex
     procedure(graph_edge_interface)             , deferred :: edge
     procedure(graph_vertex_support_interface)   , deferred :: vertex_support
     procedure(graph_edge_support_interface)     , deferred :: edge_support
     procedure(graph_incident_edges_interface)   , deferred :: incident_edges
     procedure(graph_adjacent_vertices_interface), deferred :: adjacent_vertices
  end type graph

  !===================================================================!
  ! Data interfaces
  !===================================================================!

  type, abstract :: graph_data
   contains
     procedure(data_name_interface) , deferred :: name
     procedure(data_units_interface), deferred :: units
  end type graph_data

  type, abstract, extends(graph_data) :: graph_field
   contains
     procedure(field_support_interface)       , deferred :: support
     procedure(field_num_components_interface), deferred :: num_components
     procedure(field_num_entries_interface)   , deferred :: num_entries
     procedure(field_get_vector_interface)    , deferred :: get_vector
     procedure(field_set_vector_interface)    , deferred :: set_vector
  end type graph_field

  type, abstract, extends(graph_field) :: graph_vertex_field
  end type graph_vertex_field

  type, abstract, extends(graph_field) :: graph_edge_field
  end type graph_edge_field

  type, abstract, extends(graph_data) :: graph_functional
   contains
     procedure(functional_value_interface)    , deferred :: value
     procedure(functional_set_value_interface), deferred :: set_value
  end type graph_functional

  type, abstract :: graph_reduction
   contains
     procedure(reduction_identity_interface)  , deferred :: identity
     procedure(reduction_accumulate_interface), deferred :: accumulate
     procedure(reduction_combine_interface)   , deferred :: combine
     procedure(reduction_finalize_interface)  , deferred :: finalize
     procedure(reduction_reduce_interface)    , deferred :: reduce
  end type graph_reduction

  !===================================================================!
  ! Transform interfaces
  !===================================================================!

  type, abstract :: graph_transform
   contains
     procedure(transform_check_graph_interface), deferred :: check_graph
     procedure(transform_check_data_interface) , deferred :: check_data
  end type graph_transform

  type, abstract, extends(graph_transform) :: graph_partition_transform
  end type graph_partition_transform

  type, abstract, extends(graph_partition_transform) :: graph_partitioner
   contains
     procedure(partitioner_partition_graph_interface)   , deferred :: partition_graph
     procedure(partitioner_partition_data_interface)    , deferred :: partition_data
     procedure(partitioner_num_parts_interface)         , deferred :: num_parts
     procedure(partitioner_vertex_support_interface)    , deferred :: owned_vertices
     procedure(partitioner_vertex_support_interface)    , deferred :: borrowed_vertices
     procedure(partitioner_vertex_support_interface)    , deferred :: overlap_vertices
     procedure(partitioner_edge_support_interface)      , deferred :: owned_edges
     procedure(partitioner_edge_support_interface)      , deferred :: borrowed_edges
     procedure(partitioner_edge_support_interface)      , deferred :: overlap_edges
     procedure(partitioner_owning_part_interface)       , deferred :: owning_part
     procedure(partitioner_part_to_full_interface)      , deferred :: part_to_full_vertices
     procedure(partitioner_part_to_full_interface)      , deferred :: part_to_full_edges
     procedure(partitioner_full_to_part_interface)      , deferred :: full_to_part_vertices
     procedure(partitioner_full_to_part_interface)      , deferred :: full_to_part_edges
     procedure(partitioner_borrowed_owners_interface)   , deferred :: borrowed_value_owners
     procedure(partitioner_owned_borrowers_interface)   , deferred :: owned_value_borrowers
  end type graph_partitioner

  type, abstract, extends(graph_partition_transform) :: graph_assembler
   contains
     procedure(assembler_assemble_graph_interface)     , deferred :: assemble_graph
     procedure(assembler_assemble_data_interface)      , deferred :: assemble_data
     procedure(assembler_assemble_functional_interface), deferred :: assemble_functional
     procedure(assembler_check_ownership_interface)    , deferred :: check_ownership
     procedure(assembler_check_reduction_interface)    , deferred :: check_reduction
  end type graph_assembler

  type, abstract, extends(graph_transform) :: graph_resolution_transform
  end type graph_resolution_transform

  type, abstract, extends(graph_resolution_transform) :: graph_coarsener
   contains
     procedure(coarsener_graph_interface), deferred :: coarsen_graph
     procedure(coarsener_data_interface) , deferred :: coarsen_data
  end type graph_coarsener

  type, abstract, extends(graph_resolution_transform) :: graph_refiner
   contains
     procedure(refiner_graph_interface), deferred :: refine_graph
     procedure(refiner_data_interface) , deferred :: refine_data
  end type graph_refiner

  !===================================================================!
  ! Operation interfaces
  !===================================================================!

  type, abstract :: graph_operation
   contains
     procedure(operation_name_interface), deferred :: name
  end type graph_operation

  type, abstract, extends(graph_operation) :: graph_vertex_operation
   contains
     procedure(vertex_operation_support_interface), deferred :: support
  end type graph_vertex_operation

  type, abstract, extends(graph_operation) :: graph_edge_operation
   contains
     procedure(edge_operation_support_interface), deferred :: support
  end type graph_edge_operation

  type, abstract, extends(graph_vertex_operation) :: graph_vertex_field_operation
   contains
     procedure(vertex_field_operation_apply_interface), deferred :: apply
  end type graph_vertex_field_operation

  type, abstract, extends(graph_vertex_operation) :: graph_vertex_functional_operation
   contains
     procedure(vertex_functional_operation_apply_interface), deferred :: apply
  end type graph_vertex_functional_operation

  type, abstract, extends(graph_edge_operation) :: graph_edge_field_operation
   contains
     procedure(edge_field_operation_apply_interface), deferred :: apply
  end type graph_edge_field_operation

  type, abstract, extends(graph_edge_operation) :: graph_edge_functional_operation
   contains
     procedure(edge_functional_operation_apply_interface), deferred :: apply
  end type graph_edge_functional_operation

  !===================================================================!
  ! Deferred type-bound procedure interfaces
  !===================================================================!

  abstract interface

     pure integer function entity_id_interface(this)
       import :: graph_entity
       class(graph_entity), intent(in) :: this
     end function entity_id_interface

     pure function entity_label_interface(this) result(label)
       import :: graph_entity
       class(graph_entity), intent(in) :: this
       character(len=:), allocatable :: label
     end function entity_label_interface

     pure integer function edge_tail_interface(this)
       import :: graph_edge
       class(graph_edge), intent(in) :: this
     end function edge_tail_interface

     pure integer function edge_head_interface(this)
       import :: graph_edge
       class(graph_edge), intent(in) :: this
     end function edge_head_interface

     pure logical function edge_has_head_interface(this)
       import :: graph_edge
       class(graph_edge), intent(in) :: this
     end function edge_has_head_interface

     pure integer function support_kind_interface(this)
       import :: graph_support
       class(graph_support), intent(in) :: this
     end function support_kind_interface

     pure integer function support_size_interface(this)
       import :: graph_support
       class(graph_support), intent(in) :: this
     end function support_size_interface

     pure subroutine support_ids_interface(this, ids)
       import :: graph_support
       class(graph_support), intent(in) :: this
       integer, allocatable, intent(out) :: ids(:)
     end subroutine support_ids_interface

     pure logical function support_contains_interface(this, id)
       import :: graph_support
       class(graph_support), intent(in) :: this
       integer, intent(in) :: id
     end function support_contains_interface

     pure integer function graph_id_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_id_interface

     pure integer function graph_num_vertices_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_num_vertices_interface

     pure integer function graph_num_edges_interface(this)
       import :: graph
       class(graph), intent(in) :: this
     end function graph_num_edges_interface

     subroutine graph_vertex_interface(this, id, vertex)
       import :: graph, graph_vertex
       class(graph), intent(in) :: this
       integer, intent(in) :: id
       class(graph_vertex), allocatable, intent(out) :: vertex
     end subroutine graph_vertex_interface

     subroutine graph_edge_interface(this, id, edge)
       import :: graph, graph_edge
       class(graph), intent(in) :: this
       integer, intent(in) :: id
       class(graph_edge), allocatable, intent(out) :: edge
     end subroutine graph_edge_interface

     subroutine graph_vertex_support_interface(this, selector, support, tag)
       import :: graph, graph_vertex_support
       class(graph), intent(in) :: this
       integer, intent(in) :: selector
       class(graph_vertex_support), allocatable, intent(out) :: support
       character(len=*), intent(in), optional :: tag
     end subroutine graph_vertex_support_interface

     subroutine graph_edge_support_interface(this, selector, support, tag)
       import :: graph, graph_edge_support
       class(graph), intent(in) :: this
       integer, intent(in) :: selector
       class(graph_edge_support), allocatable, intent(out) :: support
       character(len=*), intent(in), optional :: tag
     end subroutine graph_edge_support_interface

     subroutine graph_incident_edges_interface(this, vertex_id, support)
       import :: graph, graph_edge_support
       class(graph), intent(in) :: this
       integer, intent(in) :: vertex_id
       class(graph_edge_support), allocatable, intent(out) :: support
     end subroutine graph_incident_edges_interface

     subroutine graph_adjacent_vertices_interface(this, vertex_id, support)
       import :: graph, graph_vertex_support
       class(graph), intent(in) :: this
       integer, intent(in) :: vertex_id
       class(graph_vertex_support), allocatable, intent(out) :: support
     end subroutine graph_adjacent_vertices_interface

     pure function data_name_interface(this) result(name)
       import :: graph_data
       class(graph_data), intent(in) :: this
       character(len=:), allocatable :: name
     end function data_name_interface

     pure function data_units_interface(this) result(units)
       import :: graph_data
       class(graph_data), intent(in) :: this
       character(len=:), allocatable :: units
     end function data_units_interface

     subroutine field_support_interface(this, support)
       import :: graph_field, graph_support
       class(graph_field), intent(in) :: this
       class(graph_support), allocatable, intent(out) :: support
     end subroutine field_support_interface

     pure integer function field_num_components_interface(this)
       import :: graph_field
       class(graph_field), intent(in) :: this
     end function field_num_components_interface

     pure integer function field_num_entries_interface(this)
       import :: graph_field
       class(graph_field), intent(in) :: this
     end function field_num_entries_interface

     pure subroutine field_get_vector_interface(this, values)
       import :: graph_field, dp
       class(graph_field), intent(in) :: this
       real(dp), allocatable, intent(out) :: values(:)
     end subroutine field_get_vector_interface

     pure subroutine field_set_vector_interface(this, values)
       import :: graph_field, dp
       class(graph_field), intent(inout) :: this
       real(dp), intent(in) :: values(:)
     end subroutine field_set_vector_interface

     pure real(dp) function functional_value_interface(this)
       import :: graph_functional, dp
       class(graph_functional), intent(in) :: this
     end function functional_value_interface

     pure subroutine functional_set_value_interface(this, value)
       import :: graph_functional, dp
       class(graph_functional), intent(inout) :: this
       real(dp), intent(in) :: value
     end subroutine functional_set_value_interface

     pure subroutine reduction_identity_interface(this, state)
       import :: graph_reduction, dp
       class(graph_reduction), intent(in) :: this
       real(dp), allocatable, intent(out) :: state(:)
     end subroutine reduction_identity_interface

     pure subroutine reduction_accumulate_interface(this, input_graph, field, support, state, measure)
       import :: graph_reduction, graph, graph_field
       import :: graph_support, dp
       class(graph_reduction), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_field), intent(in) :: field
       class(graph_support), intent(in) :: support
       real(dp), intent(inout) :: state(:)
       class(graph_field), intent(in), optional :: measure
     end subroutine reduction_accumulate_interface

     pure subroutine reduction_combine_interface(this, left, right, combined)
       import :: graph_reduction, dp
       class(graph_reduction), intent(in) :: this
       real(dp), intent(in) :: left(:), right(:)
       real(dp), allocatable, intent(out) :: combined(:)
     end subroutine reduction_combine_interface

     pure subroutine reduction_finalize_interface(this, state, value)
       import :: graph_reduction, dp
       class(graph_reduction), intent(in) :: this
       real(dp), intent(in) :: state(:)
       real(dp), intent(out) :: value
     end subroutine reduction_finalize_interface

     subroutine reduction_reduce_interface(this, input_graph, field, support, functional, measure)
       import :: graph_reduction, graph, graph_field
       import :: graph_support, graph_functional
       class(graph_reduction), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_field), intent(in) :: field
       class(graph_support), intent(in) :: support
       class(graph_functional), allocatable, intent(out) :: functional
       class(graph_field), intent(in), optional :: measure
     end subroutine reduction_reduce_interface

     pure logical function transform_check_graph_interface(this, input_graph)
       import :: graph_transform, graph
       class(graph_transform), intent(in) :: this
       class(graph), intent(in) :: input_graph
     end function transform_check_graph_interface

     pure logical function transform_check_data_interface(this, input_graph, input_data)
       import :: graph_transform, graph, graph_data
       class(graph_transform), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_data), intent(in) :: input_data
     end function transform_check_data_interface

     subroutine partitioner_partition_graph_interface(this, full_graph)
       import :: graph_partitioner, graph
       class(graph_partitioner), intent(inout) :: this
       class(graph), intent(in) :: full_graph
     end subroutine partitioner_partition_graph_interface

     subroutine partitioner_partition_data_interface(this, full_graph, full_data, part_data)
       import :: graph_partitioner, graph, graph_data
       class(graph_partitioner), intent(inout) :: this
       class(graph), intent(in) :: full_graph
       class(graph_data), intent(in) :: full_data
       class(graph_data), allocatable, intent(out) :: part_data
     end subroutine partitioner_partition_data_interface

     pure integer function partitioner_num_parts_interface(this)
       import :: graph_partitioner
       class(graph_partitioner), intent(in) :: this
     end function partitioner_num_parts_interface

     subroutine partitioner_vertex_support_interface(this, part, support)
       import :: graph_partitioner, graph_vertex_support
       class(graph_partitioner), intent(in) :: this
       integer, intent(in) :: part
       class(graph_vertex_support), allocatable, intent(out) :: support
     end subroutine partitioner_vertex_support_interface

     subroutine partitioner_edge_support_interface(this, part, support)
       import :: graph_partitioner, graph_edge_support
       class(graph_partitioner), intent(in) :: this
       integer, intent(in) :: part
       class(graph_edge_support), allocatable, intent(out) :: support
     end subroutine partitioner_edge_support_interface

     pure integer function partitioner_owning_part_interface(this, entity_kind, entity_id)
       import :: graph_partitioner
       class(graph_partitioner), intent(in) :: this
       integer, intent(in) :: entity_kind, entity_id
     end function partitioner_owning_part_interface

     pure subroutine partitioner_part_to_full_interface(this, part, ids)
       import :: graph_partitioner
       class(graph_partitioner), intent(in) :: this
       integer, intent(in) :: part
       integer, allocatable, intent(out) :: ids(:)
     end subroutine partitioner_part_to_full_interface

     pure subroutine partitioner_full_to_part_interface(this, part, ids)
       import :: graph_partitioner
       class(graph_partitioner), intent(in) :: this
       integer, intent(in) :: part
       integer, allocatable, intent(out) :: ids(:)
     end subroutine partitioner_full_to_part_interface

     pure subroutine partitioner_borrowed_owners_interface(this, part, owner_part, owner_index)
       import :: graph_partitioner
       class(graph_partitioner), intent(in) :: this
       integer, intent(in) :: part
       integer, allocatable, intent(out) :: owner_part(:), owner_index(:)
     end subroutine partitioner_borrowed_owners_interface

     pure subroutine partitioner_owned_borrowers_interface(this, part, ptr, borrower_part, borrower_index)
       import :: graph_partitioner
       class(graph_partitioner), intent(in) :: this
       integer, intent(in) :: part
       integer, allocatable, intent(out) :: ptr(:), borrower_part(:), borrower_index(:)
     end subroutine partitioner_owned_borrowers_interface

     subroutine assembler_assemble_graph_interface(this, partitioner, part_graph, full_graph)
       import :: graph_assembler, graph_partitioner, graph
       class(graph_assembler), intent(in) :: this
       class(graph_partitioner), intent(in) :: partitioner
       class(graph), intent(in) :: part_graph
       class(graph), allocatable, intent(out) :: full_graph
     end subroutine assembler_assemble_graph_interface

     subroutine assembler_assemble_data_interface(this, partitioner, part_data, full_data)
       import :: graph_assembler, graph_partitioner, graph_data
       class(graph_assembler), intent(in) :: this
       class(graph_partitioner), intent(in) :: partitioner
       class(graph_data), intent(in) :: part_data
       class(graph_data), allocatable, intent(out) :: full_data
     end subroutine assembler_assemble_data_interface

     subroutine assembler_assemble_functional_interface(this, partitioner, part_functional, full_functional)
       import :: graph_assembler, graph_partitioner
       import :: graph_functional
       class(graph_assembler), intent(in) :: this
       class(graph_partitioner), intent(in) :: partitioner
       class(graph_functional), intent(in) :: part_functional
       class(graph_functional), allocatable, intent(out) :: full_functional
     end subroutine assembler_assemble_functional_interface

     pure logical function assembler_check_ownership_interface(this, partitioner)
       import :: graph_assembler, graph_partitioner
       class(graph_assembler), intent(in) :: this
       class(graph_partitioner), intent(in) :: partitioner
     end function assembler_check_ownership_interface

     pure logical function assembler_check_reduction_interface(this, reduction)
       import :: graph_assembler, graph_reduction
       class(graph_assembler), intent(in) :: this
       class(graph_reduction), intent(in) :: reduction
     end function assembler_check_reduction_interface

     subroutine coarsener_graph_interface(this, fine_graph, coarse_graph)
       import :: graph_coarsener, graph
       class(graph_coarsener), intent(in) :: this
       class(graph), intent(in) :: fine_graph
       class(graph), allocatable, intent(out) :: coarse_graph
     end subroutine coarsener_graph_interface

     subroutine coarsener_data_interface(this, fine_graph, fine_data, coarse_graph, coarse_data)
       import :: graph_coarsener, graph, graph_data
       class(graph_coarsener), intent(in) :: this
       class(graph), intent(in) :: fine_graph
       class(graph_data), intent(in) :: fine_data
       class(graph), intent(in) :: coarse_graph
       class(graph_data), allocatable, intent(out) :: coarse_data
     end subroutine coarsener_data_interface

     subroutine refiner_graph_interface(this, coarse_graph, fine_graph)
       import :: graph_refiner, graph
       class(graph_refiner), intent(in) :: this
       class(graph), intent(in) :: coarse_graph
       class(graph), allocatable, intent(out) :: fine_graph
     end subroutine refiner_graph_interface

     subroutine refiner_data_interface(this, coarse_graph, coarse_data, fine_graph, fine_data)
       import :: graph_refiner, graph, graph_data
       class(graph_refiner), intent(in) :: this
       class(graph), intent(in) :: coarse_graph
       class(graph_data), intent(in) :: coarse_data
       class(graph), intent(in) :: fine_graph
       class(graph_data), allocatable, intent(out) :: fine_data
     end subroutine refiner_data_interface

     pure function operation_name_interface(this) result(name)
       import :: graph_operation
       class(graph_operation), intent(in) :: this
       character(len=:), allocatable :: name
     end function operation_name_interface

     subroutine vertex_operation_support_interface(this, input_graph, support)
       import :: graph_vertex_operation, graph, graph_vertex_support
       class(graph_vertex_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_vertex_support), allocatable, intent(out) :: support
     end subroutine vertex_operation_support_interface

     subroutine edge_operation_support_interface(this, input_graph, support)
       import :: graph_edge_operation, graph, graph_edge_support
       class(graph_edge_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_edge_support), allocatable, intent(out) :: support
     end subroutine edge_operation_support_interface

     subroutine vertex_field_operation_apply_interface(this, input_graph, input_data, output)
       import :: graph_vertex_field_operation, graph
       import :: graph_data, graph_vertex_field
       class(graph_vertex_field_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_data), intent(in), optional :: input_data(:)
       class(graph_vertex_field), allocatable, intent(out) :: output
     end subroutine vertex_field_operation_apply_interface

     subroutine vertex_functional_operation_apply_interface(this, input_graph, input_data, output)
       import :: graph_vertex_functional_operation, graph
       import :: graph_data, graph_functional
       class(graph_vertex_functional_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_data), intent(in), optional :: input_data(:)
       class(graph_functional), allocatable, intent(out) :: output
     end subroutine vertex_functional_operation_apply_interface

     subroutine edge_field_operation_apply_interface(this, input_graph, input_data, output)
       import :: graph_edge_field_operation, graph
       import :: graph_data, graph_edge_field
       class(graph_edge_field_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_data), intent(in), optional :: input_data(:)
       class(graph_edge_field), allocatable, intent(out) :: output
     end subroutine edge_field_operation_apply_interface

     subroutine edge_functional_operation_apply_interface(this, input_graph, input_data, output)
       import :: graph_edge_functional_operation, graph
       import :: graph_data, graph_functional
       class(graph_edge_functional_operation), intent(in) :: this
       class(graph), intent(in) :: input_graph
       class(graph_data), intent(in), optional :: input_data(:)
       class(graph_functional), allocatable, intent(out) :: output
     end subroutine edge_functional_operation_apply_interface

  end interface

end module abstract_graph_types
