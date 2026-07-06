!=====================================================================!
! Expected-failure driver: an illegal tag must die at the door with its
! name, never be silently reinterpreted. DIAGONAL passed positionally
! where mode is expected is exactly the collision B4 closes - the part
! range is disjoint from the mode range, so validation refuses it.
!=====================================================================!

program tag_refusal

  use iso_fortran_env  , only : dp => REAL64
  use class_gmsh_loader, only : gmsh_loader
  use class_mesh       , only : mesh
  use class_assembler  , only : assembler, DIAGONAL

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

  ! a part tag in the mode position: this call must not return
  call fvm % get_jacobian_residual_product(w, v, DIAGONAL)

  write(*,'(a)') "tag refusal did not fire"

end program tag_refusal
