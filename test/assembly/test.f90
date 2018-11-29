!=====================================================================!
! Test loading of mesh and mesh pre-processing
!=====================================================================!

program test_mesh

  use iso_fortran_env   , only : dp => REAL64
  use class_gmsh_loader , only : gmsh_loader
  use class_mesh        , only : mesh
  use class_assembler   , only : assembler, solve_conjugate_gradient

  implicit none
  
  character(len=*)  , parameter   :: filename = "delaunay.msh"
  class(gmsh_loader), allocatable :: gmsh_loader_obj
  class(mesh)       , allocatable :: mesh_obj

  ! Solution parameters
  real(dp) , parameter   :: max_tol = 1.0d-4
  integer  , parameter   :: max_it = 10000
  real(dp) , allocatable :: x(:)
  class(assembler), allocatable :: FVMAssembler
  integer :: npts
  
  ! Create a mesh object
  allocate(gmsh_loader_obj, source =  gmsh_loader(filename))
  allocate(mesh_obj, source = mesh(gmsh_loader_obj))
  deallocate(gmsh_loader_obj)

  ! Create an assembler object for assembling the linear system
  ! Geometry and meshing
  allocate(FVMAssembler, source = assembler(mesh_obj))

  npts = FVMAssembler % num_state_vars

  block

    real(8), allocatable :: phinoskew(:), phic(:), ss(:)

    ! Initial solution guess
    allocate(phic(npts))
    phic = 0;
    allocate(phinoskew(npts))    
    call random_number(phinoskew)    
    call FVMassembler % write_solution("mesh-in.dat", phinoskew)

    allocate(ss(npts)); ss = 0.0d0;
    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phinoskew)
    print *, "delphi", norm2(phic-phinoskew)
    call FVMassembler % write_solution("mesh-noskew.dat", phinoskew)
    phic = phinoskew
    call FVMAssembler % get_skew_source(ss, phic)
    print *, "norm of skew source before solution", norm2(ss)
    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phic)
    print *, "delphi", norm2(phic-phinoskew)
    call FVMassembler % write_solution("mesh-skew.dat", phic)

    call FVMAssembler % get_skew_source(ss, phic)
    print *, "norm of skew source before solution", norm2(ss)
    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phic)
    print *, "delphi", norm2(phic-phinoskew)


    call FVMAssembler % get_skew_source(ss, phic)
    print *, "norm of skew source before solution", norm2(ss)
    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phic)
    print *, "delphi", norm2(phic-phinoskew)


    call FVMAssembler % get_skew_source(ss, phic)
    print *, "norm of skew source before solution", norm2(ss)
    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phic)
    print *, "delphi", norm2(phic-phinoskew)


    call FVMAssembler % get_skew_source(ss, phic)
    print *, "norm of skew source before solution", norm2(ss)
    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phic)
    print *, "delphi", norm2(phic-phinoskew)



    call FVMassembler % write_solution("mesh-diff.dat", phic-phinoskew)
    
    ! Writes the mesh for tecplot
    call FVMassembler % write_solution("mesh-out.dat", phic)  

  end block

  stop
  
!!$  ! Iterate until skew source is zero
!!$
!!$  ! Create a solver object to solve the linear system  
!!$  allocate(x(npts))
!!$  call random_number(x)
!!$  !print *, 'xinit', x
!!$  call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ssx)
!!$  !print *, 'solution', x
!!$
!!$  print *, max_it, max_tol
!!$
!!$  ! Writes the mesh for tecplot
!!$  call FVMassembler % write_solution("mesh.dat", x)

  deallocate(mesh_obj)
  deallocate(FVMAssembler)
  deallocate(x)

  contains

end program test_mesh
