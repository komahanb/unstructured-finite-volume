!=====================================================================!
! Test loading of mesh and mesh pre-processing
!=====================================================================!

program test_mesh

  use class_mesh_loader , only : mesh_loader
  use class_mesh        , only : mesh
!!$  use class_file        , only : file
!!$  use class_string      , only : string

  implicit none

  character(len=*), parameter :: filename = 'rectangle.msh'
  type(mesh_loader) :: gmsh_loader
  type(mesh)        :: grid

!!$  type(file)        :: ofile
!!$
!!$  type(string), allocatable :: lines(:)
!!$  integer :: iline
!!$
!!$  ofile = file(filename)
!!$  allocate(lines(ofile % get_num_lines()))
!!$  call ofile % open()
!!$  do iline = 1, 42
!!$     call ofile % read_line(lines(iline))
!!$  end do
!!$  call ofile % close()
!!$  call lines % print()
!!$  deallocate(lines)
!!$
!!$  call ofile % read_lines(lines)
!!$  call lines % print()

  ! Create a mesh loader for mesh file
  gmsh_loader = mesh_loader(filename)

  ! Get the mesh using loader object
  grid = mesh(gmsh_loader)
  !call grid % to_string()

end program test_mesh
