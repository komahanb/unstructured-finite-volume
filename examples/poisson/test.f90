!=====================================================================!
! Test loading of mesh and mesh pre-processing
!=====================================================================!

program test_mesh

  use iso_fortran_env          , only : dp => REAL64
  use class_string             , only : string
  use class_gmsh_loader        , only : gmsh_loader
  use class_mesh               , only : mesh
  use class_assembler          , only : assembler
  use interface_linear_solver  , only : linear_solver
  use class_conjugate_gradient , only : conjugate_gradient
  use class_gauss_jacobi       , only : gauss_jacobi
  use class_gauss_seidel       , only : gauss_seidel
  use class_sor                , only : sor

  implicit none

!!$  
!!$  type(string)         , parameter   :: fname = string("rectangle-40.msh")

  character(len=*)     , parameter   :: filename = "square-40.msh"
  class(gmsh_loader)   , allocatable :: gmsh
  class(mesh)          , allocatable :: grid
  class(linear_solver) , allocatable :: solver
  class(assembler)     , allocatable :: FVMAssembler

  ! Arguments for "solver-type", "input-mesh",

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

    ! Also supply
    ! allocate(FVMAssembler, source = assembler(grid,physics_list)) 
    ! physics with tags Assembler combines Geometry and Physics ( EQNS
    ! + BC) to provide linear/nonlinear systems
    
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

    rmse =  sum(abs(fv-fhatv)**2)
    !rmse = rmse + sum(abs(ff-fhatf)**2)
    !rmse = rmse + sum(abs(fc-fhatc)**2)
    rmse = sqrt(rmse/dble(FVMAssembler % grid % num_vertices))

    print *, "rmse = ", rmse
    print *, "cell volume = ", &
         & sum(FVMAssembler % grid % cell_volumes)/dble(FVMAssembler % grid % num_cells), &
         & minval(FVMAssembler % grid % cell_volumes), &
         & maxval(FVMAssembler % grid % cell_volumes)
    
!!$
!!$    rmse = sqrt(&
!!$         & sum(abs(fc-fhatc)**2) + sum(abs(ff-fhatf)**2) + sum(abs(fv-fhatv)**2)&
!!$         & /dble(npts) &
!!$         & )
    

!!$
!!$    rmse = sqrt(sum(error**2.0d0)/dble(FVMassembler % num_state_vars))
!!$    print *, "rmse cell center", rmse, rmse_vertices

    deallocate(fhatc, fhatv, fhatf)
    deallocate(fc, fv, ff)
    deallocate(solver)
    
  end block cg_solver

  stop

  sor_solver : block

    real(dp) , parameter   :: max_tol     = 100.0d0*epsilon(1.0d0)
    integer  , parameter   :: max_it      = 100
    integer  , parameter   :: print_level = -1
    real(dp) , allocatable :: x(:)
    integer :: i

    allocate(solver, &
         & source      = sor( &
         & FVAssembler = FVMassembler, &
         & omega       = 1.8545d0, &
         & max_tol     = max_tol, &
         & max_it      = max_it, &
         & print_level = print_level))

    ! Solve using Seidel method
    call solver % solve(x)
    print *, 'sor solution = '
    do i = 1, min(10, size(x))
       print *, i,  x(i)
    end do

    ! Writes the mesh for tecplot
    call FVMassembler % write_solution("poission-sor-40.dat", x)

    deallocate(x)   
    deallocate(solver)
    
  end block sor_solver
  
  seidel_solver : block

    real(dp) , parameter   :: max_tol     = 100.0d0*epsilon(1.0d0)
    integer  , parameter   :: max_it      = 100
    integer  , parameter   :: print_level = -1
    real(dp) , allocatable :: x(:)
    integer :: i

    allocate(solver, &
         & source      = gauss_seidel( &
         & FVAssembler = FVMassembler, &
         & max_tol     = max_tol, &
         & max_it      = max_it, &
         & print_level = print_level))

    ! Solve using Seidel method
    call solver % solve(x)
    print *, 'seidel solution = '
    do i = 1, min(10, size(x))
       print *, i,  x(i)
    end do

    ! Writes the mesh for tecplot
    call FVMassembler % write_solution("poission-seidel-40.dat", x)

    deallocate(x)   
    deallocate(solver)
    
  end block seidel_solver

  jacobi_solver : block

    real(dp) , parameter   :: max_tol     = 100.0d0*epsilon(1.0d0)
    integer  , parameter   :: max_it      = 100
    integer  , parameter   :: print_level = -1
    real(dp) , allocatable :: x(:)
    integer :: i

    allocate(solver, &
         & source      = gauss_jacobi( &
         & FVAssembler = FVMassembler, &
         & max_tol     = max_tol, &
         & max_it      = max_it, &
         & print_level = print_level))

    ! Solve using Jacobi method
    call solver % solve(x)
    print *, 'jacobi solution = '
    do i = 1, min(10, size(x))
       print *, i,  x(i)
    end do

    ! Writes the mesh for tecplot
    call FVMassembler % write_solution("poission-jacobi-40.dat", x)

    deallocate(x)   
    deallocate(solver)
    
  end block jacobi_solver
  
  deallocate(grid)
  deallocate(FVMAssembler)
  
contains

  subroutine get_exact_solution(fexact, x)

    ! Arguments
    real(dp)              , intent(in)  :: x(:,:)
    real(dp), allocatable , intent(out) :: fexact(:)

    ! Locals
    real(dp), parameter :: PI = 4.0d0*atan(1.0d0)
    real(dp), parameter :: alpha = 16.0d0/(PI**4.0d0)
    integer  :: i, ii, jj, mm, nn, npts
            
    npts = size(x, dim=2)
    allocate(fexact(npts))
    fexact = 0.0d0

    ! Exact solution as a function of coordinate
    do i = 1, npts

       ! Evalute f(x) as a summation
       do ii = 0, 99
          do jj = 0, 99
             mm = 2*ii+1
             nn = 2*jj+1
             fexact(i) = fexact(i) + &
                  & alpha*sin(dble(mm)*x(1,i)*pi) *sin(dble(nn)*x(2,i)*pi) &
                  &/(dble(mm*nn)*dble(mm**2+nn**2))
          end do
       end do

    end do

  end subroutine get_exact_solution

end program
