!=====================================================================!
! THE CONTINUOUS ADJOINT OF THE VAN DER POL OSCILLATOR, against the
! discrete one. The Lagrangian of the energy J = integral of (u^2 +
! u_t^2)/2 under the oscillator r = u_tt - nu (1 - u^2) u_t + u = 0,
!
!      L = J + integral of lambda r  +  rho_0 . g_0,
!
! has two stationarities in u. Taken before the discretization, with
! the integral by parts, it is the continuous adjoint equation and
! its natural conditions at T,
!
!      lambda_tt + nu (1 - u^2) lambda_t + lambda + u - u_tt = 0,
!      lambda(T) = 0,   lambda_t(T) + nu (1 - u^2) lambda(T) - u_t(T) = 0,
!
! formed here by the library from the density (u^2 + u_t^2)/2 +
! lambda r: its stationarity in u along t, and the boundary terms of
! the integration by parts. The oscillator with its initial data and
! the adjoint with its terminal conditions are one two-point problem
! on [0, T], solved as one block over every instant by dirk(2), the
! same family at both ends. The sensitivity by the continuous adjoint
! is then dJ/dnu = integral of lambda (dr/dnu) with dr/dnu the partial
! of r in nu, formed by the library too.
!
! Taken after the discretization, the stationarity of L_h in u_h is
! the discrete adjoint: the multipliers of the discrete rows, which
! minimize forms by the reverse sweep, and kappa = -dJ_h/dnu exactly.
! The two agree as the instants are refined, and the expansion along
! nu gives dJ_h/dnu a third way. The instant count is read from the
! command line (default 21).
!
! The terminal conditions occupy the rows of lambda's components at
! the first instant, which no family determines, their columns at T:
! the two-point structure, which one block over every instant admits
! and a chain of blocks does not. The condition on lambda_t reads
! lambda through the coefficient nu (1 - u^2), evaluated at zero
! state in the prescribed row; with lambda(T) = 0 by the other
! condition the row is exact.
!
!      instants   continuous adjoint   discrete adjoint  -kappa
!         21        -0.6495986           -0.6568391
!         81        -0.6435874           -0.6437486
!        321        -0.6432127           -0.6432196
!
! and lambda(t_k) against mu_k / h at the instants of the multistep
! families agree to 2 per cent at 21 instants.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program van_der_pol_continuous_adjoint

  use util_precision    , only : dp
  use util_verbosity    , only : set_verbosity
  use operation_manifold, only : continuous_manifold, discrete_manifold, interval, instants, parameter, expansion
  use operation_field   , only : continuous_field, discrete_field, integral, &
       &                         operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual, only : continuous_residual, discrete_residual, operator(+), operator(*)
  use operation_family  , only : dirk, bdf, adams, chain
  use view_expression   , only : expression_view

  implicit none

  real(dp), parameter :: T_final = 2.0_dp, nu_design = 1.0_dp, u_initial = 2.0_dp, v_initial = 0.0_dp
  integer  :: num_instants = 21

  ! the continuous adjoint: one manifold, two unknowns, one block
  type(continuous_manifold) :: omega, domega_0, domega_T
  type(discrete_manifold)   :: omega_h
  type(continuous_field)    :: t, nu, u, lambda, r, density, adjoint, u0, v0, rest, sensitivity
  type(continuous_field), allocatable :: terminal(:)
  type(discrete_field)      :: estimate, solution, lambda_h, functional
  type(continuous_residual) :: system, g_0, g_T, L
  type(discrete_residual)   :: L_h
  type(continuous_field)    :: rho_0, rho_T

  ! the discrete adjoint: the manifold of the sensitivity contract
  type(continuous_manifold) :: omega_d, domega_d
  type(discrete_manifold)   :: omega_dh
  type(continuous_field)    :: t_d, nu_d, u_d, lambda_r, lambda_g, kappa, r_d, u0_d, v0_d, energy, rest_d
  type(discrete_field)      :: estimate_d, solution_d, kappa_h, adjoint_h, energy_h
  type(continuous_residual) :: r_term, g_term, d_term, L_d
  type(discrete_residual)   :: L_dh

  type(expression_view) :: view
  real(dp) :: J(1), sensitivities(2)
  character(len=32) :: item
  integer :: k

  if (command_argument_count() >= 1) then
     call get_command_argument(1, item)
     read(item, *) num_instants
  end if

  ! THE CONTINUOUS ADJOINT. Omega = [0, T] x {nu}; its faces at 0 and
  ! at T; Omega_h the instants, the design expanded to order zero
  omega    = continuous_manifold(time=interval(0.0_dp, T_final), design=parameter('nu', nu_design))
  domega_0 = omega % boundary(time=0.0_dp)
  domega_T = omega % boundary(time=T_final)
  omega_h  = omega % discretize(time=instants(num_instants), design=expansion(order=0))

  t      = omega % coordinate('t')
  nu     = omega % coordinate('nu')
  u      = omega % unknown('u', [t, nu])
  lambda = omega % unknown('lambda', [t, nu])

  ! the oscillator, the density of the Lagrangian, and its
  ! stationarity in u: the adjoint equation, with the boundary terms
  ! of the integration by parts as the conditions at T
  r        = u % derivative([t, t]) - nu*(1.0_dp - u*u)*u % derivative([t]) + u
  density  = (u*u + u % derivative([t])**2)/2.0_dp + lambda*r
  adjoint  = density % stationarity(u, t)
  terminal = density % boundary_terms(u, t)
  view = expression_view(adjoint % graph(1))
  print '(a)', 'the adjoint equation, the stationarity of (u^2 + u_t^2)/2 + lambda r in u:'
  print '(a)', '   ' // view % formula()
  do k = 1, size(terminal)
     view = expression_view(terminal(k) % graph(1))
     print '(a, i0, a)', '   condition at T, term ', k, ': ' // view % formula()
  end do

  ! the initial data of u at 0, the natural conditions of lambda at
  ! T, the two-point system on Omega, and its Lagrangian
  u0  = u - u_initial
  v0  = u % derivative([t]) - v_initial
  system = continuous_residual(omega,    [r, adjoint])
  g_0    = continuous_residual(domega_0, [u0, v0])
  g_T    = continuous_residual(domega_T, terminal)
  rho_0  = g_0 % multiplier('rho_0', [nu])
  rho_T  = g_T % multiplier('rho_T', [nu])
  L      = system + rho_0*g_0 + rho_T*g_T

  ! one block over every instant by dirk(2): the state marches from
  ! its data at 0, the adjoint from its conditions at T, together
  L_h = L % discretize(omega_h, time=chain([dirk(2)], from=[1]))

  rest     = continuous_field(omega, [u_initial, 0.0_dp])
  estimate = rest % discretize(omega_h)
  if (this_image() == 1 .and. num_instants <= 21) call set_verbosity(1)
  call L_h % minimize(estimate, solution)

  ! dJ/dnu by the continuous adjoint: the integral of lambda dr/dnu
  sensitivity = integral(lambda * r % partial(nu), over=omega)
  functional  = sensitivity % at(solution)
  call functional % values(J)
  sensitivities(1) = J(1)
  lambda_h = solution % fields(['lambda'])

  ! THE DISCRETE ADJOINT: the statement of the sensitivity contract,
  ! one design coordinate, the chain dirk(2), bdf(2), adams(2) over
  ! the same instants, the expansion of order 1
  omega_d  = continuous_manifold(time=interval(0.0_dp, T_final), design=parameter('nu', nu_design))
  domega_d = omega_d % boundary(time=0.0_dp)
  omega_dh = omega_d % discretize(time=instants(num_instants), design=expansion(order=1))
  t_d  = omega_d % coordinate('t')
  nu_d = omega_d % coordinate('nu')
  u_d  = omega_d % unknown('u', [t_d, nu_d])
  r_d  = u_d % derivative([t_d, t_d]) - nu_d*(1.0_dp - u_d*u_d)*u_d % derivative([t_d]) + u_d
  u0_d = u_d - u_initial
  v0_d = u_d % derivative([t_d]) - v_initial
  energy   = integral((u_d*u_d + u_d % derivative([t_d])**2)/2.0_dp, over=omega_d % time())
  r_term   = continuous_residual(omega_d,  [r_d])
  g_term   = continuous_residual(domega_d, [u0_d, v0_d])
  d_term   = continuous_residual(omega_d % design(), [nu_d - nu_design])
  lambda_r = r_term % multiplier('adjoint', [t_d, nu_d])
  lambda_g = g_term % multiplier('lambda', [nu_d])
  kappa    = d_term % multiplier('kappa')
  L_d      = energy + lambda_r*r_term + lambda_g*g_term + kappa*d_term
  L_dh     = L_d % discretize(omega_dh, time=chain([(dirk(2), k = 1, 3), (bdf(2), k = 4, 7), &
       &                                            (adams(2), k = 8, num_instants)], from=[(k, k = 1, num_instants)]))
  rest_d     = continuous_field(omega_d, [u_initial])
  estimate_d = rest_d % discretize(omega_dh)
  call set_verbosity(0)
  call L_dh % minimize(estimate_d, solution_d)
  kappa_h   = solution_d % fields(['kappa'])
  adjoint_h = solution_d % fields(['adjoint'])
  energy_h  = energy % at(solution_d)
  call energy_h % values(sensitivities)
  sensitivities(1) = J(1)

  ! the adjoint at the instants: the continuous lambda(t_k), and the
  ! discrete multiplier of the equation's row at t_k over the step,
  ! at the instants of the multistep families
  print '(a)', ' t        lambda continuous     mu discrete / h'
  do k = 8, num_instants, max(1, (num_instants - 8) / 6)
     print '(f6.2, 2es22.10)', omega_h % instant(k), lambda_h % value(k, 1), &
          & adjoint_h % value(k, 1) / omega_dh % step(k)
  end do
  print '(a, i0, a)', 'dJ/dnu at ', num_instants, ' instants:'
  print '(a, es24.16)', '   continuous adjoint, integral of lambda dr/dnu = ', sensitivities(1)
  print '(a, es24.16)', '   discrete adjoint, -kappa                      = ', -kappa_h % value(1, 1)
  print '(a, es24.16)', '   expansion along nu                            = ', sensitivities(2)

end program van_der_pol_continuous_adjoint
