module class_mesh_loader

  ! import dependencies
  use iso_fortran_env, only : dp => real64
  use class_file     , only : file
  use class_string, only : string

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
    
    ! Mesh tag
    character(len=*), parameter :: BEGIN_MESH = "$MeshFormat"
    character(len=*), parameter :: END_MESH   = "$EndMeshFormat"  
    integer :: idx_start_mesh
    integer :: idx_end_mesh

    ! Physical_Names tag
    character(len=*), parameter :: BEGIN_PHYSICAL_NAMES = "$PhysicalNames"
    character(len=*), parameter :: END_PHYSICAL_NAMES   = "$EndPhysicalNames"  
    integer :: idx_start_physical_names
    integer :: idx_end_physical_names

    ! Nodes tag
    character(len=*), parameter   :: BEGIN_NODES = "$Nodes"
    character(len=*), parameter   :: END_NODES   = "$EndNodes"  
    integer :: idx_start_nodes
    integer :: idx_end_nodes

    ! Elements tag
    character(len=*), parameter   :: BEGIN_ELEMENTS = "$Elements"
    character(len=*), parameter   :: END_ELEMENTS   = "$EndElements"  
    integer :: idx_start_elements
    integer :: idx_end_elements

    integer :: iline, num_lines

    write(*,'(a,a)') "Loading mesh file : ", this % file % filename
    
    ! Load the mesh into memory
    call this % file % read_lines(lines)
    ! call lines % print()  

    ! Extract start and end indices of different mesh tags used by
    ! GMSH
    num_lines = size(lines)
    do concurrent(iline = 1:num_lines)

       ! Find mesh start and end
       if (index(lines(iline) % str, BEGIN_MESH) .eq. 1) then
          idx_start_mesh = iline
       end if
       if (index(lines(iline) % str, END_MESH) .eq. 1) then
          idx_end_mesh = iline
       end if

       ! Find physical_names start and end
       if (index(lines(iline) % str, BEGIN_PHYSICAL_NAMES) .eq. 1) then
          idx_start_physical_names = iline
       end if
       if (index(lines(iline) % str, END_PHYSICAL_NAMES) .eq. 1) then
          idx_end_physical_names = iline
       end if

       ! Find nodes start and end
       if (index(lines(iline) % str, BEGIN_NODES) .eq. 1) then
          idx_start_nodes = iline
       end if
       if (index(lines(iline) % str, END_NODES) .eq. 1) then
          idx_end_nodes = iline
       end if

       ! Find elements start and end
       if (index(lines(iline) % str, BEGIN_ELEMENTS) .eq. 1) then
          idx_start_elements = iline
       end if
       if (index(lines(iline) % str, END_ELEMENTS) .eq. 1) then
          idx_end_elements = iline
       end if

    end do

    write(*,'(a,i8,i8)') "mesh           : " , idx_start_mesh           , idx_end_mesh
    write(*,'(a,i8,i8)') "physical names : " , idx_start_physical_names , idx_end_physical_names
    write(*,'(a,i8,i8)') "nodes          : " , idx_start_nodes          , idx_end_nodes
    write(*,'(a,i8,i8)') "elements       : " , idx_start_elements       , idx_end_elements   

    ! Process nodes
!!$    num_vertices = idx_end_nodes - idx_start_nodes
!!$    allocate(vertices(num_vertices))
!!$    vertex_numbers = lines(iline) % tokenize(" ", tokens, num_tokens)

    ! process elements

    deallocate(lines)

  end subroutine get_mesh_data

end module class_mesh_loader
