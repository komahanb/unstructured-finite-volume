module class_mesh_loader

  ! import dependencies
  use iso_fortran_env, only : dp => real64
  use file_class     , only : file
  use string_class, only : string

  implicit none

  private
  public :: mesh_loader

  !-------------------------------------------------------------------!
  ! Derived type for mesh load
  !-------------------------------------------------------------------!
  
  type :: mesh_loader

     type(file) :: file ! mesh file

   contains

     ! Type bound procedure that returns all information needed for
     ! mesh creation
     procedure :: get_mesh_data

  end type mesh_loader

  !-------------------------------------------------------------------!
  ! interface to construct a mesh_loader
  !-------------------------------------------------------------------!

  interface mesh_loader
     module procedure create
  end interface mesh_loader

contains

  type(mesh_loader) function create(filename) result (this)
    
    type(character(*)), intent(in) :: filename

    this % file = file(filename)

  end function create
  
  !====================================================================!
  ! Supply all information needed to create a mesh object
  !====================================================================!

  subroutine get_mesh_data(this, &
       & num_vertices, vertices, vertex_numbers, vertex_tags, & 
       & num_edges, edge_vertices, num_edge_vertices, edge_tags, &
       & num_faces, face_vertices, num_face_vertices, face_tags, &
       & num_cells, cell_vertices, num_cell_vertices, cell_tags)

    ! Arguments
    class(mesh_loader)  , intent(in)   :: this

    integer , intent(out)              :: num_vertices
    real(dp), intent(out), allocatable :: vertices(:,:)
    integer , intent(out), allocatable :: vertex_numbers(:)
    integer , intent(out), allocatable :: vertex_tags(:)

    integer, intent(out)              :: num_faces
    integer, intent(out), allocatable :: face_vertices(:,:)
    integer, intent(out), allocatable :: num_face_vertices(:)
    integer, intent(out), allocatable :: face_tags(:)

    integer, intent(out)              :: num_edges
    integer, intent(out), allocatable :: edge_vertices(:,:)
    integer, intent(out), allocatable :: num_edge_vertices(:)
    integer, intent(out), allocatable :: edge_tags(:)

    integer, intent(out)              :: num_cells
    integer, intent(out), allocatable :: cell_vertices(:,:)
    integer, intent(out), allocatable :: num_cell_vertices(:)
    integer, intent(out), allocatable :: cell_tags(:)

    ! Local
    type(string), allocatable, dimension(:) :: lines

    ! Load the mesh into memory
    call this % file % read_lines(lines)
    
    ! Process nodes

    ! process elements

    deallocate(lines)

  end subroutine get_mesh_data

end module class_mesh_loader
