module interface_mesh_loader

  use iso_fortran_env , only : dp => real64
  use class_string    , only : string

  implicit none

  private
  public :: mesh_loader

  !-------------------------------------------------------------------!
  ! Abstract type for mesh loading
  !-------------------------------------------------------------------!
  
  type, abstract :: mesh_loader

   contains

     ! Type bound procedure that returns all information needed for
     ! mesh creation
     procedure(get_mesh_data_interface), deferred :: get_mesh_data

  end type mesh_loader

  interface

     subroutine get_mesh_data_interface(this, &
          & num_vertices, vertex_numbers, vertex_tags , vertices ,  & 
          & num_edges   , edge_numbers  , edge_tags   , edge_vertices , num_edge_vertices , &
          & num_faces   , face_numbers  , face_tags   , face_vertices , num_face_vertices , &
          & num_cells   , cell_numbers  , cell_tags   , cell_vertices , num_cell_vertices , &
          & num_tags    , tag_numbers , tag_info )

       import mesh_loader
       import dp
       import string
       
       ! Arguments
       class(mesh_loader)  , intent(in)   :: this

       ! Vertices
       integer , intent(out)              :: num_vertices
       integer , intent(out), allocatable :: vertex_numbers(:)
       integer , intent(out), allocatable :: vertex_tags(:)
       real(dp), intent(out), allocatable :: vertices(:,:)

       ! Edges
       integer, intent(out)              :: num_edges
       integer, intent(out), allocatable :: edge_numbers(:)
       integer, intent(out), allocatable :: edge_tags(:)
       integer, intent(out), allocatable :: edge_vertices(:,:)
       integer, intent(out), allocatable :: num_edge_vertices(:)

       ! Faces    
       integer, intent(out)              :: num_faces
       integer, intent(out), allocatable :: face_numbers(:)
       integer, intent(out), allocatable :: face_tags(:)
       integer, intent(out), allocatable :: face_vertices(:,:)
       integer, intent(out), allocatable :: num_face_vertices(:)

       ! Cells
       integer, intent(out)              :: num_cells
       integer, intent(out), allocatable :: cell_numbers(:)
       integer, intent(out), allocatable :: cell_tags(:)
       integer, intent(out), allocatable :: cell_vertices(:,:)
       integer, intent(out), allocatable :: num_cell_vertices(:)
       
       ! Tagging boundaries and domain with integers/strings
       integer     , intent(out)              :: num_tags              
       integer     , allocatable, intent(out) :: tag_numbers(:)       
       type(string), allocatable, intent(out) :: tag_info(:)

     end subroutine get_mesh_data_interface

  end interface

contains

end module interface_mesh_loader
