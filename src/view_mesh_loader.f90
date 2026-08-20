!=====================================================================!
! The loader contract: whoever claims to be a mesh source must hand
! over the raw incidence lists - vertices with coordinates, then the
! edge -> vertex, face -> vertex and cell -> vertex tables, plus the
! physical tags. That is a graph described by its edges, delivered
! before any geometry is measured; class_mesh does the wiring and
! the measuring. One deferred procedure, mesh_data, carries it
! all - a gmsh file on disk and an in-memory array satisfy the same
! contract.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module view_mesh_loader

  use iso_fortran_env , only : dp => real64
  use util_string    , only : string

  implicit none

  private
  public :: mesh_loader

  !-------------------------------------------------------------------!
  ! Every mesh loader extends this abstract type.
  !-------------------------------------------------------------------!

  type, abstract :: mesh_loader

   contains

     ! This deferred procedure returns all the information needed for
     ! mesh creation.
     procedure(mesh_data_interface), deferred :: mesh_data

  end type mesh_loader

  interface

     !================================================================!
     ! This is the one deferred procedure of the contract: hand over
     ! the raw incidence lists, the type codes, and the tag table that
     ! describe the mesh graph.
     !================================================================!

     subroutine mesh_data_interface(this, &
          & num_vertices, vertex_numbers, vertex_tags , vertices ,  & 
          & num_edges   , edge_numbers  , edge_tags   , edge_vertices , num_edge_vertices , &
          & num_faces   , face_numbers  , face_tags   , face_vertices , num_face_vertices , &
          & num_cells   , cell_numbers  , cell_tags   , cell_vertices , num_cell_vertices , &
          & cell_types  , face_types    , edge_types  , &
          & num_tags    , tag_numbers   , tag_physical_dimensions, tag_info )

       import mesh_loader
       import dp
       import string
       
       ! This is the loader being asked.
       class(mesh_loader)  , intent(in)   :: this

       ! These arguments carry the vertices.
       integer , intent(out)              :: num_vertices
       integer , intent(out), allocatable :: vertex_numbers(:)
       integer , intent(out), allocatable :: vertex_tags(:)
       real(dp), intent(out), allocatable :: vertices(:,:)

       ! These arguments carry the edges.
       integer, intent(out)              :: num_edges
       integer, intent(out), allocatable :: edge_numbers(:)
       integer, intent(out), allocatable :: edge_tags(:)
       integer, intent(out), allocatable :: edge_vertices(:,:)
       integer, intent(out), allocatable :: num_edge_vertices(:)
       integer, intent(out), allocatable :: edge_types(:)

       ! These arguments carry the faces.
       integer, intent(out)              :: num_faces
       integer, intent(out), allocatable :: face_numbers(:)
       integer, intent(out), allocatable :: face_tags(:)
       integer, intent(out), allocatable :: face_vertices(:,:)
       integer, intent(out), allocatable :: num_face_vertices(:)
       integer, intent(out), allocatable :: face_types(:)

       ! These arguments carry the cells.
       integer, intent(out)              :: num_cells
       integer, intent(out), allocatable :: cell_numbers(:)
       integer, intent(out), allocatable :: cell_tags(:)
       integer, intent(out), allocatable :: cell_vertices(:,:)
       integer, intent(out), allocatable :: num_cell_vertices(:)
       integer, intent(out), allocatable :: cell_types(:)

       ! These arguments tag the boundaries and the domain with
       ! integers and strings.
       integer     , intent(out)              :: num_tags
       integer     , allocatable, intent(out) :: tag_numbers(:)
       integer     , allocatable, intent(out) :: tag_physical_dimensions(:)
       type(string), allocatable, intent(out) :: tag_info(:)

     end subroutine mesh_data_interface

  end interface

contains

end module view_mesh_loader
