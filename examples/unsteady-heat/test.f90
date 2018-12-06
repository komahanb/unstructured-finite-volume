!=====================================================================!
! Test loading of mesh and mesh pre-processing
!=====================================================================!

program test_mesh

  use iso_fortran_env             , only : dp => REAL64
  use class_gmsh_loader           , only : gmsh_loader
  use class_mesh                  , only : mesh
  use interface_physics           , only : physics
  use interface_linear_solver     , only : linear_solver
  use class_assembler             , only : assembler
  use class_physics_unsteady_heat , only : unsteady_heat
  use class_conjugate_gradient    , only : conjugate_gradient

  implicit none

  character(len=*)     , parameter   :: filename = "square-10msh"
  class(gmsh_loader)   , allocatable :: gmsh
  class(mesh)          , allocatable :: grid
  class(linear_solver) , allocatable :: solver
  class(assembler)     , allocatable :: FVMAssembler
  class(physics)       , allocatable :: heat

  meshing : block

    ! Create a mesh object
    allocate(gmsh, source =  gmsh_loader(filename))
    allocate(grid, source = mesh(gmsh))
    deallocate(gmsh)
    
  end block meshing

  physics: block

    allocate(heat, source = unsteady_heat(diff_coeff=1.0d0))   
    
  end block physics

  assembly : block

    ! Create an assembler object for assembling the linear system
    ! Geometry and meshing
    allocate(FVMAssembler, source = assembler(grid, heat))

  end block assembly

  cg_solver : block

    real(dp) , parameter   :: max_tol     = 100.0d0*epsilon(1.0d0)
    integer  , parameter   :: max_it      = 100
    integer  , parameter   :: print_level = -1
    real(dp) , allocatable :: fhatc(:), fhatf(:), fhatv(:)
    real(dp) , allocatable :: fc(:), ff(:), fv(:)
    real(dp) , allocatable :: ec(:), ef(:), ev(:)
    real(dp) :: rmse
    integer  :: i, npts

    allocate(solver, &
         & source      = conjugate_gradient( &
         & FVAssembler = FVMassembler, &
         & max_tol     = max_tol, &
         & max_it      = max_it, &
         & print_level = print_level))

    ! Solve using solver method
    call solver % solve(fhatc)
    print *, 'cg solution = '
    do i = 1, min(10, size(fhatc))
       print *, i,  fhatc(i)
    end do

    ! Writes the mesh for tecplot
    call FVMassembler % write_solution("poission-cg-40.dat", fhatc)
    call FVMAssembler % evaluate_vertex_flux(fhatv, fhatc)
    call FVMAssembler % evaluate_face_flux(fhatf, fhatc)
    
    ! Get exact cell center values of solution
    call get_exact_solution(ff, FVMAssembler % grid % face_centers(1:2,:))
    call get_exact_solution(fc, FVMAssembler % grid % cell_centers(1:2,:))
    call get_exact_solution(fv, FVMAssembler % grid % vertices(1:2,:))

    ! Write tecplot output of error
    call FVMassembler % create_vector(ec)
    ec = abs(fc-fhatc)
    call FVMassembler % write_solution("poission-ec-40.dat", ec)

    print *, "num_faces      ", FVMAssembler % grid % num_faces
    print *, "num_cells      ", FVMAssembler % grid % num_cells
    print *, "num_vertices   ", FVMAssembler % grid % num_vertices
    print *, "num_state_vars ", FVMAssembler % num_state_vars

    npts = FVMAssembler % grid % num_cells &
         & + FVMAssembler % grid % num_vertices  &
         & + FVMAssembler % grid % num_faces

    rmse = sum(abs(fv-fhatv)**2)
    rmse = sqrt(rmse/dble(FVMAssembler % grid % num_vertices))

    print *, "rmse = ", rmse
    print *, "cell volume = ", &
         & sum(FVMAssembler % grid % cell_volumes)/dble(FVMAssembler % grid % num_cells), &
         & minval(FVMAssembler % grid % cell_volumes), &
         & maxval(FVMAssembler % grid % cell_volumes)
    
    deallocate(fhatc, fhatv, fhatf)
    deallocate(fc, fv, ff)
    deallocate(solver)
    
  end block cg_solver

  deallocate(heat)  
  deallocate(grid)
  deallocate(FVMAssembler)
  
contains

end program
