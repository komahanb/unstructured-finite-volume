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
!!$  type(string)         , parameter   :: fname = string("square-40.msh")

  character(len=*)     , parameter   :: filename = "square-10.msh"
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
    real(dp) , allocatable :: x(:), exact(:), error(:)
    real(dp) :: rmse
    integer :: i

    allocate(solver, &
         & source      = conjugate_gradient( &
         & FVAssembler = FVMassembler, &
         & max_tol     = max_tol, &
         & max_it      = max_it, &
         & print_level = print_level))

    ! Solve using solver method
    call solver % solve(x)
    print *, 'cg solution = '
    do i = 1, min(10, size(x))
       print *, i,  x(i)
    end do

    ! Writes the mesh for tecplot
    call FVMassembler % write_solution("poission-cg-10.dat", x)

    call get_exact_solution(FVMassembler, exact)
    call FVMassembler % write_solution("poission-exact-10.dat", exact)

    call FVMassembler % create_vector(error)
    error = exact-x
    call FVMassembler % write_solution("poission-error-10.dat", error)

    rmse = sqrt(sum(error**2.0d0)/dble(FVMassembler % num_state_vars))
    print *, "rmse", rmse

    deallocate(x)   
    deallocate(exact)
    deallocate(error)
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

!!$
!!$  ! Exact solution for homogenous poisson equation
!!$  exact: block
!!$
!!$  end block exact
!!$
!!$  ! Compute the root mean square error
!!$  rmse: block
!!$
!!$    rmse = sqrt(sum(error**2.0d0)/dble(ncells))
!!$  end block rmse
  
!!$  write(*,*) ncells, error
  
  deallocate(grid)
  deallocate(FVMAssembler)
  
contains

  subroutine get_exact_solution(FVAssembler, exact)

    ! Arguments
    class(assembler)      , intent(in)  :: FVAssembler
    real(dp), allocatable , intent(out) :: exact(:)

    ! Locals
    real(dp), parameter :: PI = 4.0d0*atan(1.0d0)
    real(dp), parameter :: alpha = 16.0d0/(PI**4.0d0)
    integer     :: icell, ii, jj, mm, nn
    real(dp)    :: x(2)

    ! Allocate space for exact solution
    call FVAssembler % create_vector(exact, 0.0d0)

    loop_cells: do icell = 1, FVAssembler % num_state_vars

       ! Get the cell center vertex coordinates
       x = FVAssembler % grid % cell_centers(1:2,icell)

       ! Exact solution as a function of coordinate
       do ii = 0, 9
          do jj = 0, 9
             mm = 2*ii+1
             nn = 2*jj+1
             exact(icell) = exact(icell) + &
                  & alpha*sin(dble(mm)*x(1)*pi) *sin(dble(nn)*x(2)*pi) &
                  &/(dble(mm*nn)*dble(mm**2+nn**2))
          end do
       end do

    end do loop_cells

  end subroutine get_exact_solution

end program
