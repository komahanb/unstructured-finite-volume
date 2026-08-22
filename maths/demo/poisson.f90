! Poisson on the unit square through the constitution: a material,
! four walls, a mesh; the compiled stencil closes under GMRES.
! Manufactured solution u = sin(pi x) sin(pi y), k = 1, f = 2 pi^2 u,
! u = 0 on every wall. Reports the volume-weighted L2 and the max
! error per mesh, the observed order, and wall-clock seconds for
! mesh derivation, operator assembly and the solve.
program poisson
  use iso_fortran_env, only : dp => REAL64, int64
  use view_mesh, only : mesh
  use view_mesh_builder, only : mesh_from_gmsh
  use operation_conduction, only : conduction
  use operation_robin_condition, only : robin_condition, dirichlet
  use operation_diffusion, only : diffusion_stencil
  use operation_stencil, only : stencil
  use operation_gmres, only : gmres
  use field_stored, only : stored_field
  implicit none

  character(len=*), parameter :: dir = 'mesh/'
  real(dp), parameter :: pi = 3.141592653589793_dp
  integer :: sizes(4) = [10, 20, 40, 80]
  type(mesh) :: m
  type(stencil) :: op
  type(gmres) :: gm
  type(stored_field) :: centres, volumes
  real(dp), allocatable :: xc(:), vol(:), u(:), f(:), g(:), r(:), rhs(:), x(:)
  real(dp) :: l2(4), linf(4), t_mesh, t_asm, t_solve, achieved, sgn
  integer(int64) :: c0, c1, c2, c3, rate
  integer :: k, n, i, ne
  character(len=8) :: tag

  print '(a)', 'TABLE poisson'
  print '(a)', 'N & cells & faces & L2 error & order & max error & order & t_mesh & t_assemble & t_solve \\'
  do k = 1, 4
     write(tag, '(i0)') sizes(k)
     call system_clock(c0, rate)
     m = mesh_from_gmsh(dir // 'square-' // trim(tag) // '.msh')
     call system_clock(c1)
     n  = m % num_vertices()
     ne = m % num_edges()

     block
       type(robin_condition) :: walls(4)
       walls = [dirichlet('BoundaryLeft', 0.0_dp), dirichlet('BoundaryRight', 0.0_dp), &
            &   dirichlet('BoundaryTop', 0.0_dp), dirichlet('BoundaryBottom', 0.0_dp)]
       op = diffusion_stencil(m, conduction(1.0_dp), walls)
     end block
     call system_clock(c2)

     centres = m % cell_centre(); call centres % real_vector(xc)
     volumes = m % cell_volume(); call volumes % real_vector(vol)
     allocate(u(n), f(n))
     do i = 1, n
        u(i) = sin(pi * xc(3 * (i - 1) + 1)) * sin(pi * xc(3 * (i - 1) + 2))
        f(i) = 2.0_dp * pi * pi * u(i) * vol(i)
     end do

     call gm % attach(op, m, m % vertex_set(), n)
     gm % tolerance = 1.0d-12
     gm % max_iterations = 5000
     call gm % constant(g)

     ! the sign the compiled balance gives the laplacian, read off
     ! the exact solution rather than assumed: A u + g ~ sgn * (k lap u) V
     call gm % matvec(u, r)
     r = r + g
     sgn = sign(1.0_dp, sum(r * (-f)))
     rhs = sgn * (-f) - g

     allocate(x(n)); x = 0.0_dp
     call gm % solve(rhs, x, achieved)
     call system_clock(c3)

     l2(k)   = sqrt(sum(vol * (x - u) ** 2))
     linf(k) = maxval(abs(x - u))
     t_mesh  = real(c1 - c0, dp) / real(rate, dp)
     t_asm   = real(c2 - c1, dp) / real(rate, dp)
     t_solve = real(c3 - c2, dp) / real(rate, dp)

     if (k == 1) then
        print '(i3,a,i6,a,i6,a,es9.2,a,a,a,es9.2,a,a,a,f7.3,a,f8.3,a,f8.3,a)', sizes(k), ' & ', n, ' & ', ne, &
             & ' & ', l2(k), ' & ', '--', ' & ', linf(k), ' & ', '--', ' & ', t_mesh, ' & ', t_asm, ' & ', t_solve, ' \\'
     else
        print '(i3,a,i6,a,i6,a,es9.2,a,f5.2,a,es9.2,a,f5.2,a,f7.3,a,f8.3,a,f8.3,a)', sizes(k), ' & ', n, ' & ', ne, &
             & ' & ', l2(k), ' & ', log(l2(k-1) / l2(k)) / log(2.0_dp), ' & ', linf(k), ' & ', &
             & log(linf(k-1) / linf(k)) / log(2.0_dp), ' & ', t_mesh, ' & ', t_asm, ' & ', t_solve, ' \\'
     end if
     print '(a,es9.2,a,f5.1,a,i0,a)', '   achieved residual ', achieved, '  sign ', sgn, '  entries ', size(op % rows)
     deallocate(u, f, x)
  end do
end program poisson
