!=====================================================================!
! Test loading of mesh and mesh pre-processing
!=====================================================================!

program test_mesh

  use iso_fortran_env          , only : dp => REAL64
  use class_gmsh_loader        , only : gmsh_loader
  use class_mesh               , only : mesh
  use class_assembler          , only : assembler
  use interface_linear_solver  , only : linear_solver
  use class_conjugate_gradient , only : conjugate_gradient

  implicit none
  
  character(len=*)     , parameter   :: filename = "rectangle.msh"
  class(gmsh_loader)   , allocatable :: gmsh
  class(mesh)          , allocatable :: grid
  class(linear_solver) , allocatable :: CG
  class(assembler)     , allocatable :: FVMAssembler

  meshing : block

    ! Create a mesh object
    allocate(gmsh, source =  gmsh_loader(filename))
    allocate(grid, source = mesh(gmsh))
    deallocate(gmsh)
    
  end block meshing

  assembly : block
    ! Create an assembler object for assembling the linear system
    ! Geometry and meshing
    allocate(FVMAssembler, source = assembler(grid))
    
  end block assembly

  solver : block

    real(dp) , parameter   :: max_tol = 1.0d-4
    integer  , parameter   :: max_it = 10000

    real(dp) , allocatable :: x(:)
    integer :: npts
    npts = FVMAssembler % num_state_vars
    
    allocate(CG, source = conjugate_gradient(&
         & FVAssembler=FVMassembler, &
         & max_tol=max_tol, &
         & max_it=max_it))
    
    deallocate(CG)

  end block solver

  deallocate(grid)
  deallocate(FVMAssembler)

!!$
!!$
!!$  block
!!$
!!$    real(8), allocatable :: phinoskew(:), phic(:), ss(:)
!!$
!!$    ! Initial solution guess
!!$    allocate(phic(npts))
!!$    phic = 0;
!!$    allocate(phinoskew(npts))    
!!$    call random_number(phinoskew)    
!!$    call FVMassembler % write_solution("mesh-in.dat", phinoskew)
!!$
!!$    allocate(ss(npts)); ss = 0.0d0;
!!$    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phinoskew)
!!$    print *, "delphi", norm2(phic-phinoskew)
!!$    call FVMassembler % write_solution("mesh-noskew.dat", phinoskew)
!!$    phic = phinoskew
!!$    call FVMAssembler % get_skew_source(ss, phic)
!!$    print *, "norm of skew source before solution", norm2(ss)
!!$    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phic)
!!$    print *, "delphi", norm2(phic-phinoskew)
!!$    call FVMassembler % write_solution("mesh-skew.dat", phic)
!!$
!!$    call FVMAssembler % get_skew_source(ss, phic)
!!$    print *, "norm of skew source before solution", norm2(ss)
!!$    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phic)
!!$    print *, "delphi", norm2(phic-phinoskew)
!!$
!!$
!!$    call FVMAssembler % get_skew_source(ss, phic)
!!$    print *, "norm of skew source before solution", norm2(ss)
!!$    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phic)
!!$    print *, "delphi", norm2(phic-phinoskew)
!!$
!!$
!!$    call FVMAssembler % get_skew_source(ss, phic)
!!$    print *, "norm of skew source before solution", norm2(ss)
!!$    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phic)
!!$    print *, "delphi", norm2(phic-phinoskew)
!!$
!!$
!!$    call FVMAssembler % get_skew_source(ss, phic)
!!$    print *, "norm of skew source before solution", norm2(ss)
!!$    call solve_conjugate_gradient(FVMAssembler, max_it, max_tol, ss, phic)
!!$    print *, "delphi", norm2(phic-phinoskew)
!!$
!!$
!!$
!!$    call FVMassembler % write_solution("mesh-diff.dat", phic-phinoskew)
!!$    
!!$    ! Writes the mesh for tecplot
!!$    call FVMassembler % write_solution("mesh-out.dat", phic)  
!!$
!!$  end block
!!$
!!$  stop
  
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

  contains

end program test_mesh
