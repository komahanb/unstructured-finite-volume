!=====================================================================!
! THE HEAT EQUATION ON THE UNIT SQUARE AND THE UNIT CUBE BY LAGRANGE
! ELEMENTS: the order of the error of the elements of degree one.
!
!      u_t - div(k grad u) - f = 0      on [0, T] x B,   k = 1 + x / 2,
!
! B the unit square or the unit cube, u = 0 on every face of B, u at
! t = 0 the exact solution's, T = 1/2, and f the law applied to the
! exact solution
!
!      u = sin(pi x) sin(pi y) cos(t),   times sin(pi z) on the cube.
!
! The law is written in flux form, each flux k u_x differentiated
! along its own coordinate, and the region is discretized at the
! corners of a mesh, whose cells are the elements:
!
!      omega_h = omega % discretize(time=points(n), space=corners(mesh(file)))
!      L_h     = L % discretize(omega_h, time=..., space=finite_element(degree=1))
!
! The arguments: the mesh, default square_q8.msh, and the instant
! count, default 21. The meshes come from square.geo and cube.geo:
! quadrangles (square_q*) or triangles (square_t*), hexahedra (cube_h*)
! or tetrahedra (cube_t*), n cells on each side. Measured, the error in
! L2 over [0, T] x B:
!
!      quadrangles, 81 instants, n = 4 to 32:  1.98e-2, 4.94e-3, 1.23e-3, 3.08e-4
!      triangles, 81 instants, n = 4 to 16:    5.16e-2, 1.37e-2, 3.47e-3
!      hexahedra, 41 instants, n = 2 to 8:     6.41e-2, 1.54e-2, 3.81e-3
!      tetrahedra, 41 instants, n = 2 to 8:    1.67e-1, 6.58e-2, 1.92e-2
!
! falling by 4 per halving, order 2, the order of the elements of
! degree one; on the tetrahedra by 2.5 and 3.4, still approaching it.
! The digits are the same on 1, 2 and 4 images.
!=====================================================================!

program heat_conduction_study

  use util_precision          , only : dp
  use operation_manifold      , only : continuous_manifold, discrete_manifold, interval, region, points, mesh, corners
  use operation_field         , only : continuous_field, discrete_field, integral
  use operation_field         , only : sin, cos, operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual      , only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family        , only : dirk, bdf, adams, chain
  use operation_finite_element, only : finite_element

  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp), T_final = 0.5_dp
  integer  :: num_instants = 21, k, d
  character(len=64) :: file = 'square_q8.msh', item

  type(continuous_manifold) :: omega, initial
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, x, y, z, u, u_exact, exact, conductivity, flux_law, flux_exact, square
  type(continuous_field)    :: q(3), q_exact(3), coordinate(3)
  type(continuous_residual) :: r, ic, Lagrangian
  type(continuous_residual), allocatable :: face(:)
  type(discrete_field)      :: estimate, solution, error, functional
  type(discrete_residual)   :: L_h
  real(dp) :: J(1)

  if (command_argument_count() >= 1) call get_command_argument(1, file)
  if (command_argument_count() >= 2) then
     call get_command_argument(2, item)
     read(item, *) num_instants
  end if

  ! the region of the geometry the mesh was made from
  omega   = continuous_manifold(time=interval(0.0_dp, T_final), &
       &                        space=region(trim(merge('cube.geo  ', 'square.geo', file(1:4) == 'cube'))))
  omega_h = omega % discretize(time=points(num_instants), space=corners(mesh(trim(file))))
  d       = omega_h % dimension
  t       = omega % coordinate('t')
  x       = omega % coordinate('x')
  y       = omega % coordinate('y')
  z       = omega % coordinate('z')
  initial = omega % boundary(t, 0.0_dp)

  conductivity = 1.0_dp + x / 2.0_dp
  u_exact = sin(pi * x) * sin(pi * y) * cos(t)
  if (d == 3) u_exact = u_exact * sin(pi * z)
  exact = continuous_field(omega, [u_exact])
  u     = omega % unknown('u')

  ! the law in flux form: u_t - div(k grad u) - f, f the law applied to
  ! the exact solution
  coordinate = [x, y, z]
  do k = 1, d
     q(k)       = conductivity * u % derivative([coordinate(k)])
     q_exact(k) = conductivity * u_exact % derivative([coordinate(k)])
  end do
  flux_law   = q(1) % derivative([x]) + q(2) % derivative([y])
  flux_exact = q_exact(1) % derivative([x]) + q_exact(2) % derivative([y])
  if (d == 3) then
     flux_law   = flux_law   + q(3) % derivative([z])
     flux_exact = flux_exact + q_exact(3) % derivative([z])
  end if
  r  = continuous_residual(omega, [u % derivative([t]) - flux_law - (u_exact % derivative([t]) - flux_exact)])
  ic = continuous_residual(initial, [u - u_exact])

  ! u = 0 on every face of the square or the cube
  allocate(face(2 * d))
  face(1) = continuous_residual(omega % boundary(x, 0.0_dp), [u])
  face(2) = continuous_residual(omega % boundary(x, 1.0_dp), [u])
  face(3) = continuous_residual(omega % boundary(y, 0.0_dp), [u])
  face(4) = continuous_residual(omega % boundary(y, 1.0_dp), [u])
  if (d == 3) then
     face(5) = continuous_residual(omega % boundary(z, 0.0_dp), [u])
     face(6) = continuous_residual(omega % boundary(z, 1.0_dp), [u])
  end if
  Lagrangian = r + ic
  do k = 1, 2 * d
     Lagrangian = Lagrangian + face(k)
  end do

  L_h = Lagrangian % discretize(omega_h, &
       & time=chain([(dirk(2), k = 1, 4), (bdf(2), k = 5, 7), (adams(2), k = 8, num_instants)], &
       &            from=[(k, k = 1, num_instants)]), space=finite_element(degree=1))
  estimate = exact % discretize(omega_h)
  call L_h % minimize(estimate, solution)

  error      = solution % fields(['u']) - estimate
  square     = integral((u - u_exact)**2 / 2.0_dp, over=omega)
  functional = square % at(solution)
  call functional % values(J)
  if (this_image() == 1) print '(a, a, i0, a, es12.4, a, es12.4)', trim(file), ', ', num_instants, &
       & ' instants: error at the nodes ', error % norm(), ', in L2 over the manifold ', sqrt(2.0_dp * J(1))

end program heat_conduction_study
