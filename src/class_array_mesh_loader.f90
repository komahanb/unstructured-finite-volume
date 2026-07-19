!=====================================================================!
! A mesh loader whose data is arrays, filled in memory by whoever
! builds them - no file behind it. This is the front door for meshes
! the code derives itself (a refinement of a loaded mesh, say): fill
! the components, hand the loader to the mesh constructor, and the
! derived mesh arrives with everything a loaded mesh has - faces,
! centres, volumes, neighbours - built by the same machinery.
!
! The components mirror the loader contract's outputs one for one;
! get_mesh_data does nothing but copy them out.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module class_array_mesh_loader

  use iso_fortran_env      , only : dp => real64
  use class_string         , only : string
  use interface_mesh_loader, only : mesh_loader

  implicit none

  private
  public :: array_mesh_loader

  type, extends(mesh_loader) :: array_mesh_loader

     ! vertices
     integer               :: num_vertices = 0
     integer , allocatable :: vertex_numbers(:)
     integer , allocatable :: vertex_tags(:)
     real(dp), allocatable :: vertices(:,:)

     ! edges
     integer              :: num_edges = 0
     integer, allocatable :: edge_numbers(:)
     integer, allocatable :: edge_tags(:)
     integer, allocatable :: edge_vertices(:,:)
     integer, allocatable :: num_edge_vertices(:)
     integer, allocatable :: edge_types(:)

     ! faces (the boundary faces; interior faces are derived)
     integer              :: num_faces = 0
     integer, allocatable :: face_numbers(:)
     integer, allocatable :: face_tags(:)
     integer, allocatable :: face_vertices(:,:)
     integer, allocatable :: num_face_vertices(:)
     integer, allocatable :: face_types(:)

     ! cells
     integer              :: num_cells = 0
     integer, allocatable :: cell_numbers(:)
     integer, allocatable :: cell_tags(:)
     integer, allocatable :: cell_vertices(:,:)
     integer, allocatable :: num_cell_vertices(:)
     integer, allocatable :: cell_types(:)

     ! the tag table
     integer                   :: num_tags = 0
     integer     , allocatable :: tag_numbers(:)
     integer     , allocatable :: tag_physical_dimensions(:)
     type(string), allocatable :: tag_info(:)

   contains

     procedure :: get_mesh_data

  end type array_mesh_loader

contains

  !===================================================================!
  ! The contract, honoured by copying the stored arrays out
  !===================================================================!

  subroutine get_mesh_data(this, &
       & num_vertices, vertex_numbers, vertex_tags , vertices ,  &
       & num_edges   , edge_numbers  , edge_tags   , edge_vertices , num_edge_vertices , &
       & num_faces   , face_numbers  , face_tags   , face_vertices , num_face_vertices , &
       & num_cells   , cell_numbers  , cell_tags   , cell_vertices , num_cell_vertices , &
       & cell_types  , face_types    , edge_types  , &
       & num_tags    , tag_numbers   , tag_physical_dimensions, tag_info )

    class(array_mesh_loader), intent(in) :: this

    integer , intent(out)              :: num_vertices
    integer , intent(out), allocatable :: vertex_numbers(:)
    integer , intent(out), allocatable :: vertex_tags(:)
    real(dp), intent(out), allocatable :: vertices(:,:)

    integer, intent(out)              :: num_edges
    integer, intent(out), allocatable :: edge_numbers(:)
    integer, intent(out), allocatable :: edge_tags(:)
    integer, intent(out), allocatable :: edge_vertices(:,:)
    integer, intent(out), allocatable :: num_edge_vertices(:)
    integer, intent(out), allocatable :: edge_types(:)

    integer, intent(out)              :: num_faces
    integer, intent(out), allocatable :: face_numbers(:)
    integer, intent(out), allocatable :: face_tags(:)
    integer, intent(out), allocatable :: face_vertices(:,:)
    integer, intent(out), allocatable :: num_face_vertices(:)
    integer, intent(out), allocatable :: face_types(:)

    integer, intent(out)              :: num_cells
    integer, intent(out), allocatable :: cell_numbers(:)
    integer, intent(out), allocatable :: cell_tags(:)
    integer, intent(out), allocatable :: cell_vertices(:,:)
    integer, intent(out), allocatable :: num_cell_vertices(:)
    integer, intent(out), allocatable :: cell_types(:)

    integer     , intent(out)              :: num_tags
    integer     , allocatable, intent(out) :: tag_numbers(:)
    integer     , allocatable, intent(out) :: tag_physical_dimensions(:)
    type(string), allocatable, intent(out) :: tag_info(:)

    num_vertices   = this % num_vertices
    vertex_numbers = this % vertex_numbers
    vertex_tags    = this % vertex_tags
    vertices       = this % vertices

    num_edges         = this % num_edges
    edge_numbers      = this % edge_numbers
    edge_tags         = this % edge_tags
    edge_vertices     = this % edge_vertices
    num_edge_vertices = this % num_edge_vertices
    edge_types        = this % edge_types

    num_faces         = this % num_faces
    face_numbers      = this % face_numbers
    face_tags         = this % face_tags
    face_vertices     = this % face_vertices
    num_face_vertices = this % num_face_vertices
    face_types        = this % face_types

    num_cells         = this % num_cells
    cell_numbers      = this % cell_numbers
    cell_tags         = this % cell_tags
    cell_vertices     = this % cell_vertices
    num_cell_vertices = this % num_cell_vertices
    cell_types        = this % cell_types

    num_tags                = this % num_tags
    tag_numbers             = this % tag_numbers
    tag_physical_dimensions = this % tag_physical_dimensions
    tag_info                = this % tag_info

  end subroutine get_mesh_data

end module class_array_mesh_loader
