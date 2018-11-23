!=====================================================================!
! Test loading of mesh and mesh pre-processing
!=====================================================================!

program test_mesh
  
  use class_mesh             , only : mesh
  use class_gmsh_loader      , only : gmsh_loader
  use class_test_mesh_loader , only : test_mesh_loader
  use class_file             , only : file
  use class_string           , only : string

  implicit none
  
  !===================================================================!
  ! Test the functionalities of Class GMSH_LOADER
  !===================================================================!
  
  test_gmsh_delaunay: block

    character(len=*)  , parameter   :: filename = 'rectangle.msh'
    class(gmsh_loader), allocatable :: gmsh_loader_obj
    class(mesh)       , allocatable :: mesh_obj

    ! Create a mesh loader for mesh file
    allocate(gmsh_loader_obj, source =  gmsh_loader(filename))
    allocate(mesh_obj, source = mesh(gmsh_loader_obj))
    call mesh_obj % to_string()
    deallocate(mesh_obj)
    deallocate(gmsh_loader_obj)

  end block test_gmsh_delaunay

  !===================================================================!
  ! Test the functionalities of Class GMSH_LOADER
  !===================================================================!
  
  test_gmsh_triangle: block
    
    character(len=*)  , parameter   :: filename = 'frontal.msh'
    class(gmsh_loader), allocatable :: gmsh_loader_obj
    class(mesh)       , allocatable :: mesh_obj

    ! Create a mesh loader for mesh file
    allocate(gmsh_loader_obj, source =  gmsh_loader(filename))
    allocate(mesh_obj, source = mesh(gmsh_loader_obj))
    call mesh_obj % to_string()
    deallocate(mesh_obj)
    deallocate(gmsh_loader_obj)

  end block test_gmsh_triangle
  
  test_mesh1 : block

    character(len=*), parameter :: coord_file = 'coordinates_10.input'
    character(len=*), parameter :: elem_file  = 'elements_10.input'

    class(test_mesh_loader), allocatable :: test_mesh_loader_obj
    class(mesh), allocatable             :: mesh_obj

    allocate(test_mesh_loader_obj, source = test_mesh_loader(coord_file, elem_file))
    allocate(mesh_obj, source = mesh(test_mesh_loader_obj))
    !call mesh_obj % to_string()
    deallocate(mesh_obj)
    deallocate(test_mesh_loader_obj)

  end block test_mesh1

  test_mesh2 : block

    character(len=*), parameter :: coord_file = 'coordinates_20.input'
    character(len=*), parameter :: elem_file  = 'elements_20.input'

    class(test_mesh_loader), allocatable :: test_mesh_loader_obj
    class(mesh), allocatable             :: mesh_obj

    allocate(test_mesh_loader_obj, source = test_mesh_loader(coord_file, elem_file))
    allocate(mesh_obj, source = mesh(test_mesh_loader_obj))
    !call mesh_obj % to_string()
    deallocate(mesh_obj)
    deallocate(test_mesh_loader_obj)

  end block test_mesh2

  stop

  !===================================================================!
  ! Test the functionalities of Class String and Class File
  !===================================================================!

  test_file_string : block

    character(len=*), parameter :: filename = 'rectangle.msh'
    type(string), allocatable :: lines(:)
    type(file)                :: file_obj
    integer                   :: iline, num_lines

    ! Create a file object
    file_obj = file(filename)
    
    ! Do a line by line read and print contents
    num_lines = file_obj % get_num_lines()
    allocate(lines(num_lines))
    call file_obj % open()
    do iline = 1, num_lines
       call file_obj % read_line(lines(iline))
    end do
    call file_obj % close()
    call lines % print()
    deallocate(lines)

    ! Read everything and print content
    call file_obj % read_lines(lines)
    call lines % print()  

  end block test_file_string

end program test_mesh
