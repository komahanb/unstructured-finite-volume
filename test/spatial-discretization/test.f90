!=====================================================================!
! Spatial order-of-accuracy study by grid refinement.
!
! Solve the 2d poisson problem  -div(grad u) = 1  on the unit square
! with u = 0 on all four sides, on a sequence of structured quad meshes
! (square-10/20/40/80, h = 1/n), and measure the discrete L2 error at
! the cell centres against the exact fourier-series solution. the
! observed order  p = log(e_coarse/e_fine)/log(h_coarse/h_fine) must
! approach 2 - the cell-centred laplacian's design order.
!
! A low-degree polynomial would be reproduced exactly by the 2nd-order
! scheme (zero error, nothing to measure); the fourier solution has all
! derivatives nonzero, so it exposes the true order.
!
! A nonzero exit (error stop) means the observed order is off.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program spatial_discretization

  use iso_fortran_env          , only : dp => REAL64
  use class_gmsh_loader        , only : gmsh_loader
  use class_mesh               , only : mesh
  use class_assembler          , only : assembler
  use class_diffusion          , only : diffusion
  use class_conjugate_gradient , only : conjugate_gradient

  implicit none

  integer , parameter :: nlevel         = 4
  integer , parameter :: levels(nlevel) = [10, 20, 40, 80]
  real(dp), parameter :: tol            = 100.0_dp*epsilon(1.0_dp)

  real(dp) :: h(nlevel), err(nlevel)
  integer  :: ncell(nlevel)
  real(dp) :: order
  integer  :: k
  logical  :: ok

  ! Refinement sequence
  do k = 1, nlevel
     call solve_poisson(levels(k), h(k), ncell(k), err(k))
  end do

  ! Convergence table
  write(*,'(a)') " spatial order of accuracy - 2d poisson on the unit square"
  write(*,'(2x,a6,2x,a8,2x,a14,2x,a8)') "n", "ncells", "L2 error", "order"
  write(*,'(2x,a6,2x,a8,2x,a14,2x,a8)') "------", "--------", "--------------", "--------"

  do k = 1, nlevel
     if (k .eq. 1) then
        write(*,'(2x,i6,2x,i8,2x,es14.6,2x,a8)') levels(k), ncell(k), err(k), "      -"
     else
        order = log(err(k-1)/err(k)) / log(h(k-1)/h(k))
        write(*,'(2x,i6,2x,i8,2x,es14.6,2x,f8.3)') levels(k), ncell(k), err(k), order
     end if
  end do

  ! The finest observed order must approach 2
  order = log(err(nlevel-1)/err(nlevel)) / log(h(nlevel-1)/h(nlevel))
  ok    = (order .gt. 1.7_dp) .and. (order .lt. 2.3_dp)

  write(*,*) "============================================="
  if (ok) then
     write(*,'(1x,a,f6.3,a)') "PASS : observed order ", order, " ~ 2 (2nd-order accurate)"
  else
     write(*,'(1x,a,f6.3,a)') "FAIL : observed order ", order, " not ~ 2"
     error stop
  end if

contains

  !===================================================================!
  ! Solve poisson on square-<n>; return h, cell count and the discrete
  ! L2 error at the cell centres against the exact solution.
  !===================================================================!

  subroutine solve_poisson(n, h, ncell, l2err)

    integer , intent(in)  :: n
    real(dp), intent(out) :: h
    integer , intent(out) :: ncell
    real(dp), intent(out) :: l2err

    class(gmsh_loader)       , allocatable :: gl
    class(mesh)              , allocatable :: grid
    class(assembler)         , allocatable :: fvm
    class(conjugate_gradient), allocatable :: cg
    real(dp)                 , allocatable :: uh(:), uex(:)
    character(len=64)                      :: meshfile, nstr

    write(nstr, '(i0)') n
    meshfile = "../square-"//trim(nstr)//".msh"

    allocate(gl  , source = gmsh_loader(trim(meshfile)))
    allocate(grid, source = mesh(gl))
    allocate(fvm , source = assembler(grid))

    ! Poisson: unit source (source=-1 is the assembler's sign convention
    ! matching the +1 fourier rhs), homogeneous dirichlet on all sides
    call fvm % set_equation(diffusion(1.0_dp, source = -1.0_dp))
    call fvm % set_dirichlet("BoundaryLeft"  , 0.0_dp)
    call fvm % set_dirichlet("BoundaryRight" , 0.0_dp)
    call fvm % set_dirichlet("BoundaryTop"   , 0.0_dp)
    call fvm % set_dirichlet("BoundaryBottom", 0.0_dp)

    allocate(cg, source = conjugate_gradient(fvm, 500, tol, 0))
    call cg % solve(uh)

    call get_exact_solution(uex, fvm % grid % cell_centers(1:2,:))

    ncell = fvm % grid % num_cells
    h     = 1.0_dp/real(n, dp)
    l2err = sqrt(sum(fvm % grid % cell_volumes*(uh - uex)**2) &
         &       / sum(fvm % grid % cell_volumes))

  end subroutine solve_poisson

  !===================================================================!
  ! Exact fourier-series solution of -div(grad u) = 1, u = 0 on the
  ! unit square (as in examples/poisson)
  !===================================================================!

  subroutine get_exact_solution(fexact, x)

    real(dp)             , intent(in)  :: x(:,:)
    real(dp), allocatable, intent(out) :: fexact(:)

    real(dp), parameter :: PI    = 4.0_dp*atan(1.0_dp)
    real(dp), parameter :: alpha = 16.0_dp/(PI**4)
    integer             :: i, ii, jj, mm, nn, npts

    npts = size(x, dim=2)
    allocate(fexact(npts))
    fexact = 0.0_dp

    do i = 1, npts
       do ii = 0, 99
          do jj = 0, 99
             mm = 2*ii + 1
             nn = 2*jj + 1
             fexact(i) = fexact(i) + &
                  & alpha*sin(real(mm,dp)*x(1,i)*PI)*sin(real(nn,dp)*x(2,i)*PI) &
                  & /(real(mm*nn,dp)*real(mm**2+nn**2,dp))
          end do
       end do
    end do

  end subroutine get_exact_solution

end program spatial_discretization
