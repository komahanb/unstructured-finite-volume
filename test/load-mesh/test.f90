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
   
  character(len=*), parameter :: coord_file = 'coordinates_10.input'
  character(len=*), parameter :: elem_file  = 'elements_10.input'

  type(test_mesh_loader)      :: test_mesh_loader_obj
  type(mesh)                  :: mesh_obj

  test_mesh_loader_obj = test_mesh_loader(coord_file, elem_file)
  mesh_obj             = mesh(test_mesh_loader_obj)
  call mesh_obj % to_string()

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

  !===================================================================!
  ! Test the functionalities of Class GMSH_LOADER
  !===================================================================!

  test_gmsh: block

    character(len=*), parameter :: filename = 'rectangle.msh'

    call test_gmsh_loader(filename)

  end block test_gmsh
  
contains
  
  subroutine test_gmsh_loader(filename)
    
    character(len=*), intent(in) :: filename
    type(gmsh_loader) :: gmsh_loader_obj
    type(mesh)        :: mesh_obj
    
    ! Create a mesh loader for mesh file
    gmsh_loader_obj =  gmsh_loader(filename)

    ! Get the mesh using loader object
    mesh_obj = mesh(gmsh_loader_obj)
    call mesh_obj % to_string()

  end subroutine test_gmsh_loader

end program test_mesh
