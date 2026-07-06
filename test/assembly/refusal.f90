!=====================================================================!
! Expected-failure driver: a REVERSE product on a system that neither
! declares its operator symmetric nor overrides transpose_product MUST
! refuse (error stop, nonzero exit). The harness asserts this program
! dies; if it returns, the refusal law is broken.
!=====================================================================!

program transpose_refusal

  use iso_fortran_env  , only : dp => REAL64
  use class_gmsh_loader, only : gmsh_loader
  use class_mesh       , only : mesh
  use class_assembler  , only : assembler
  use module_solve_mode, only : REVERSE

  implicit none

  class(gmsh_loader), allocatable :: gmsh
  class(mesh)       , allocatable :: grid
  class(assembler)  , allocatable :: fvm
  real(dp)          , allocatable :: v(:), w(:)

  allocate(gmsh, source = gmsh_loader("../box-36.msh"))
  allocate(grid, source = mesh(gmsh))
  allocate(fvm , source = assembler(grid))

  allocate(v(fvm % num_state_vars), w(fvm % num_state_vars))
  v = 1.0_dp

  ! no symmetry declared, no override: this call must not return
  call fvm % get_jacobian_residual_product(w, v, mode = REVERSE)

  write(*,'(a)') "refusal did not fire"

end program transpose_refusal
