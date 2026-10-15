!=====================================================================!
! THE VAN DER POL OSCILLATOR AS THE TEMPLATE OF A PROBLEM: the
! second-order equation
!
!      u_tt - nu (1 - u^2) u_t + u = 0,   u = u_0,  u_t = u_t_0  at t = 0,
!
! on the interval [0, T], with the energy
!
!      J = integral over [0, T] of (u^2 + u_t^2) / 2
!
! as the objective. The problem is stated once, in the order in which
! each statement reads the ones before it:
!
!   THE FIELD        the unknown u, a function on the manifold
!                    [0, T] x {nu} of the coordinates t and nu
!   THE LAWS         every law of the problem as a residual r_l that
!                    vanishes on a domain of its own: r_omega the
!                    oscillator upon the manifold, the initial face
!                    included, r_ic the pin of (u, u_t) to its datum
!                    upon the initial face t = 0, r_s = nu - nu_0 the
!                    pin of the design upon the factor {nu}
!   THE OBJECTIVE    a functional of the field, J
!   THE MULTIPLIERS  every law paired with a multiplier lambda_l of the
!                    same index
!   THE DISCRETE IMAGE chosen coordinate by coordinate, after the
!                    statement and without a change to it: along t the
!                    instants t_bar and a chain of families, dirk(2)
!                    over the first three, bdf(2) over the next four,
!                    adams(2) to the end; along nu the expansion s_bar
!                    of an order at nu_0. The instants are the
!                    manifold's image, equally spaced, or, in (e), the
!                    equidistribution, a placement the discrete image
!                    of L states, since it reads the error of J
!
! SIX PROBLEMS are then six readings of one Lagrangian, each
! differing from the one before it by one term or by the placement of
! its instants:
!
!   (a) physics          L = lambda_omega . r_omega + lambda_ic . r_ic + lambda_s . r_s
!   (b) adjoint, and
!   (c) sensitivity      L = J + lambda_omega . r_omega + lambda_ic . r_ic + lambda_s . r_s
!   (d) design           L = J + lambda_omega . r_omega + lambda_ic . r_ic
!   (e) design on the equidistribution, L of (d)
!   (f) design to a tolerance, L of (d), the count of instants free
!
! (a) The stationarity of L in the multipliers is the oscillator with
! its data and its design, solved forward for u_bar and, the
! expansion being shared, for its derivatives along nu to the order.
! Without J the stationarity in u has no source: the reverse sweep,
! which runs because r_omega is paired, returns every multiplier
! zero. At eleven instants max |lambda_omega| = 0, and u(T) =
! 0.29118558490570512, the state of (b) to every digit.
!
! (b) With J, the stationarity in u is the adjoint equation, solved in
! reverse for lambda_omega(t), the sensitivity of J to a source added
! to the oscillator at the instant t; the stationarity in the pinned
! values of u is read for lambda_ic = -(dJ/du_0, dJ/du_t_0). The
! energy is integrated by the schemes' own quadrature: the tableau's
! weights at the stages of a DIRK step, the interpolatory rule over
! the instants of a multistep step. The expansion of order zero is
! the oscillator at nu_0 alone, and the state, the multipliers and
! the energy are the same at every order: at eleven instants J =
! 3.1356779007580480 and lambda_ic = (-1.3772000188879279,
! -0.10639333936774159). The solution converges with the count: u(2)
! = 0.32331 at 401 and at 1601 instants. THE ESTIMATE of J - J_bar is
! the disagreement along time of the solution with the chain one
! order higher, dirk(3), bdf(3), adams(3), on the same instants,
! J+(P Q) - J(Q) + mu+ . R+(P Q), the stages of each step of dirk(3)
! solved from the solution's instant behind it: one contribution per
! instant. Against J by the chain of order three at 2561 instants,
! 3.0969362409469, which 5121 instants meet to 2.8e-12, its ratio to
! the error is 1.259, 1.102, 1.032, 1.0099, 1.0029 and 1.0008 at 11,
! 21, 41, 81, 161 and 321 instants, and to the difference from the
! solve of order three on the same instants 1.0215, 1.0039, 1.00068,
! 1.00013, 1.00003 and 1.00001; on instants graded toward either end
! the ratio to the error is 0.93 to 1.30 at 21 instants and 0.988 to
! 1.025 at 161.
!
! (c) In the same solve, the stationarity in nu is read for lambda_s
! = -dJ/dnu. THE AUDIT is the agreement of the two halves of L_bar:
! d^k J / dnu^k by the expansion of the state, and by the expansion of
! the adjoint, the derivative of order k - 1 of -lambda_s. At eleven
! instants and order three they are apart by 5.6e-16, 4.4e-16 and
! 1.1e-15. At 21 instants and order three the jets agree with central
! differences of the solution at nu_0 = 1 +- 0.01, +- 0.02, to the
! differences' own truncation error: du(T)/dnu = 1.2502856 against
! 1.2502527, d^2u/dnu^2 = -0.598952 against -0.598900, d^3u/dnu^3 =
! -1.9712 against -1.9710; dJ/dnu = -0.6565620 against -0.6565450,
! d^2J/dnu^2 = 1.548406 against 1.548272, d^3J/dnu^3 = 1.0211 against
! 1.0215. Two design coordinates, the mixed derivative of J and the
! symmetry of its hessian are audited in
! test/graph-scheme/van_der_pol_design_audit.f90.
!
! (d) Without r_s the design is free, and the same stationarity in nu
! is solved for nu, where dJ/dnu = 0: by Newton's iteration from nu_0,
! one forward and one reverse sweep at each iterate, the gradient
! from lambda_s and the hessian from the expansion of order one. At
! eleven instants |dJ/dnu| falls as 7.07e-1, 3.84e-2, 1.82e-3,
! 4.90e-6, 3.58e-11, quadratically, to nu = 1.5269121988 with J =
! 2.9656219 against 3.1356779 at nu_0; at order zero, by the BFGS
! update, in ten iterations to the same nu within 4.1e-11. The design
! found converges with the instants at second order, the order of
! the chain: nu = 1.49414469, 1.48607327, 1.48427255, 1.48387461 and
! 1.48378216 at 21, 41, 81, 161 and 321 instants, apart by 8.1e-3,
! 1.8e-3, 4.0e-4 and 9.2e-5. The estimate of J - J_bar at the design
! found, against J at that design by the chain of order three at 2561
! instants, has the ratio 1.259, 1.136, 1.066 and 1.028 to the error
! at 21, 41, 81 and 161 instants.
!
! (e) The Lagrangian of (d), its discrete image with the placement of
! the equidistribution: the instants of manifold_bar moved, their
! count kept, minimize solving (d), estimating and placing the
! instants again round after round. With the design free the error
! placed is that of the design found, nu* - nu_h ~ -H^(-1) E'(nu_h):
! each step's part the derivative along nu of its eta, by central
! differences of the solutions at nu_h +- delta, delta = eps^(1/3)
! max(1, |nu_h|), and H by the same differences of dJ/dnu; the step
! density that part to the power 1/3 over the step, the instants its
! quantiles, a step at most 1 + sqrt(2) times the step before it,
! Grigorieff's bound for bdf(2). A round that does not lower the sum
! of the parts is rejected and the round before returned; one that
! lowers it by less than the uncertainty of the estimate, its distance
! from the estimate by the chain of order three on the same instants,
! is returned. The instants are not design coordinates of L
! (doc/unified-residuals.md, section 16). Against 1.4837523856, nu
! by the chain of order three at 2561 instants, the error of nu found
! falls from 4.316e-2, 1.039e-2, 2.321e-3 and 5.202e-4 on equally
! spaced instants to 3.398e-2, 7.346e-3, 1.606e-3 and 3.935e-4 on the
! equidistribution, at 11, 21, 41 and 81 instants, by 21, 29, 31 and
! 24 per cent, and the estimate of nu* - nu_h reads 0.942, 0.948,
! 0.969 and 0.982 of it. With the design pinned at nu_0 the
! equidistribution of eta lowers the error of J by 33, 30, 25 and 20
! per cent at 21, 41, 81 and 161 instants. The equidistribution of eta
! with the design free lowers the error of J at the design found
! three- to tenfold, and raises that of nu found by 13 to 34 per cent
! at 41 to 161 instants: the design found reads the gradient of J,
! not its value.
!
! (f) The Lagrangian of (d) on the equidistribution with a tolerance,
! 1e-3 of the design, its count of instants free. At each count the
! rounds of (e) settle the placement; where the estimate of the error
! of the design found is not within the tolerance times the norm of
! the design, the count of steps is raised to the one at which the
! estimate, falling as the steps to the power -r, would meet it - r
! two, the order of the chain, until two counts observe the order -
! the chain's last family continued to the new count, the instants
! placed there and the next solve begun from the solution interpolated
! along time. A count whose estimate does not fall below that of the
! count before, or falls within its own uncertainty, ends the rounds
! unmet. From eleven instants the tolerance is met at 47, nu found
! 1.4849560023, its error 1.204e-3, where equally spaced instants
! meet that error at 57; at 1e-4 it is met at 147 instants, the error
! 1.171e-4, where equally spaced instants need 165. The same rule on
! the error of J at nu_0 meets 1e-3 of sum |w f| at 33 instants, the
! error 2.98e-3, against 34 equally spaced, and 1e-4 at 100, the
! error 2.15e-4, against 116.
!
! THE LAGRANGIAN of (b) and (c) is printed before its solve: J and
! every law with the multiplier it is paired with, each equation as a
! formula, then the whole as a tree.
!
! The instant count, the data u_0, u_t_0 and the order of the
! expansion are read from the command line (default 11, 2, 0 and 1).
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program van_der_pol

  use util_precision     , only : dp
  use util_verbosity     , only : set_verbosity
  use operation_manifold , only : continuous_manifold, discrete_manifold, interval, points, parameter
  use operation_expansion, only : expansion
  use operation_field    , only : continuous_field, discrete_field, integral, &
       &                          operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_residual , only : continuous_residual, discrete_residual, operator(+), operator(*), &
       &                          equidistribution
  use operation_family   , only : dirk, bdf, adams, chain
  use view_residual      , only : formula, tree

  implicit none

  real(dp), parameter :: T_final = 2.0_dp, nu_0 = 1.0_dp
  real(dp) :: u_0, u_t_0
  integer  :: num_instants, order

  type(continuous_manifold) :: manifold
  type(discrete_manifold)   :: manifold_bar
  type(continuous_field)    :: t, nu, u, J, lambda_omega, lambda_ic, lambda_s, rest
  type(continuous_residual) :: r_omega, r_ic, r_s, L
  type(discrete_residual)   :: L_bar
  type(discrete_field)      :: estimate, solution, u_bar, lambda_omega_bar, lambda_ic_bar, lambda_s_bar, J_bar, nu_bar
  real(dp), allocatable     :: eta(:)
  integer                   :: k

  num_instants = nint(argument(1, 11.0_dp))
  u_0          = argument(2, 2.0_dp)
  u_t_0        = argument(3, 0.0_dp)
  order        = nint(argument(4, 1.0_dp))

  ! THE FIELD. The manifold [0, T] x {nu}, the design nu at its datum
  ! nu_0; the coordinate functions t, nu; the unknown function u(t, nu)
  manifold = continuous_manifold(time=interval(0.0_dp, T_final), design=parameter('nu', nu_0))
  t        = manifold % coordinate('t')
  nu       = manifold % coordinate('nu')
  u        = manifold % unknown('u', [t, nu])

  ! THE LAWS, each upon its own domain: r_omega the oscillator, a
  ! function of the jet of u to order two along t, reading nu as a
  ! coordinate; r_ic upon the face t = 0, the pin of the value and the
  ! first derivative of u to their data; r_s upon the factor {nu}, the
  ! pin of the design to nu_0
  r_omega = continuous_residual(manifold, [u % derivative([t, t]) - nu*(1.0_dp - u*u)*u % derivative([t]) + u])
  r_ic    = continuous_residual(manifold % boundary(t, 0.0_dp), [u - u_0, u % derivative([t]) - u_t_0])
  r_s     = continuous_residual(manifold % design(), [nu - nu_0])

  ! THE OBJECTIVE. The energy J, the integral over the time factor
  J = integral((u*u + u % derivative([t])**2)/2.0_dp, over=manifold % time())

  ! THE MULTIPLIERS. Every law declares its multiplier, one component
  ! per equation, a function of the coordinates of its domain:
  ! lambda_omega(t, nu) the adjoint, lambda_ic(nu) one number per
  ! pinned field, lambda_s one number
  lambda_omega = r_omega % multiplier('lambda_omega', [t, nu])
  lambda_ic    = r_ic % multiplier('lambda_ic', [nu])
  lambda_s     = r_s % multiplier('lambda_s')

  ! THE DISCRETE IMAGE of the manifold: the instants t_bar = {t(1)..t(N)}
  ! along t, the expansion s_bar of the order at nu_0 along nu; and
  ! the estimate u = u_0 at every instant, the same for every problem
  manifold_bar = manifold % discretize(time=points(num_instants), design=expansion(order=order))
  rest         = continuous_field(manifold, [u_0])
  estimate     = rest % discretize(manifold_bar)

  ! the convergence of every solve printed by image 1
  if (this_image() == 1) call set_verbosity(1)

  !===================================================================!
  ! (a) PHYSICS. L without J: the stationarity in the multipliers
  ! solved for u_bar; the multipliers zero
  !===================================================================!

  L     = lambda_omega*r_omega + lambda_ic*r_ic + lambda_s*r_s
  L_bar = discrete_image(L)
  call L_bar % minimize(estimate, solution)

  u_bar            = solution % fields(['u'])
  lambda_omega_bar = solution % fields(['lambda_omega'])
  lambda_ic_bar    = solution % fields(['lambda_ic'])
  lambda_s_bar     = solution % fields(['lambda_s'])
  print '(a)', '(a) physics: L = lambda_omega . r_omega + lambda_ic . r_ic + lambda_s . r_s'
  print '(a, es24.16)', 'u(T)                    = ', u_bar % value(num_instants, 1)
  print '(a, es24.16)', 'max |lambda_omega|      = ', &
       & maxval([(abs(lambda_omega_bar % value(k, 1)), k = 1, num_instants)])
  print '(a, 3es24.16)', 'lambda_ic, lambda_s     = ', lambda_ic_bar % value(1, 1), lambda_ic_bar % value(1, 2), &
       & lambda_s_bar % value(1, 1)

  !===================================================================!
  ! (b) ADJOINT and (c) SENSITIVITY. J added: the stationarity in u
  ! solved in reverse for lambda_omega, and that in the pinned fields
  ! read for lambda_ic and lambda_s
  !===================================================================!

  L     = J + L
  L_bar = discrete_image(L)

  ! THE LAGRANGIAN: J and every law with its multiplier, as a formula
  ! and as a tree
  print '(a)', 'the Lagrangian of (b) and (c):'
  print '(a)', formula(L)
  print '(a)', tree(L)

  call L_bar % minimize(estimate, solution)

  ! (b) u_bar and lambda_omega_bar at every tenth of the interval, and
  ! the multipliers of the initial pin
  u_bar            = solution % fields(['u'])
  lambda_omega_bar = solution % fields(['lambda_omega'])
  lambda_ic_bar    = solution % fields(['lambda_ic'])
  print '(a)', '(b) adjoint: L = J + lambda_omega . r_omega + lambda_ic . r_ic + lambda_s . r_s'
  do k = 1, num_instants
     if (mod(k - 1, max(1, (num_instants - 1) / 10)) == 0) then
        print '(a, f8.4, a, es16.8, a, es16.8)', 't = ', manifold_bar % instant(k), '   u = ', u_bar % value(k, 1), &
             & '   lambda_omega = ', lambda_omega_bar % value(k, 1)
     end if
  end do
  J_bar = J % at(solution)
  print '(a, es24.16)', 'u(T)                    = ', u_bar % value(num_instants, 1)
  print '(a, es24.16)', 'J                       = ', J_bar % value(1, 1)
  print '(a, es24.16)', 'lambda_ic = -dJ / du_0   = ', lambda_ic_bar % value(1, 1)
  print '(a, es24.16)', 'lambda_ic = -dJ / du_t_0 = ', lambda_ic_bar % value(1, 2)

  ! (b) THE ESTIMATE of J - J_bar, one contribution per instant
  eta = L_bar % disagreement(solution)
  print '(a, es24.16)', 'estimate of J - J_bar    = ', sum(eta)

  ! (c) the multiplier of the design pin, and THE AUDIT: d^k J / dnu^k
  ! by the expansion of the state, against the derivative of order
  ! k - 1 of -lambda_s, the expansion of the adjoint
  lambda_s_bar = solution % fields(['lambda_s'])
  print '(a)', '(c) sensitivity: the same solve'
  print '(a, es24.16)', 'lambda_s = -dJ / dnu     = ', lambda_s_bar % value(1, 1)
  if (order >= 1) print '(a)', 'k, d^k J / dnu^k by the expansion of the state, by that of the adjoint, their difference'
  do k = 1, order
     if (k > 1) lambda_s_bar = lambda_s_bar % derivative([nu])
     print '(i3, 2es24.16, es12.2)', k, J_bar % value(1, 1 + k), -lambda_s_bar % value(1, 1), &
          & J_bar % value(1, 1 + k) + lambda_s_bar % value(1, 1)
  end do

  ! (c) d^k u(T) / dnu^k by the expansion of the state
  do k = 1, order
     u_bar = u_bar % derivative([nu])
     print '(a, i0, a, i0, a, es24.16)', 'd^', k, ' u(T) / dnu^', k, '        = ', u_bar % value(num_instants, 1)
  end do

  !===================================================================!
  ! (d) DESIGN. r_s removed: nu free, solved for where dJ/dnu = 0
  !===================================================================!

  L     = J + lambda_omega*r_omega + lambda_ic*r_ic
  L_bar = discrete_image(L)
  call L_bar % minimize(estimate, solution)

  nu_bar = solution % fields(['nu'])
  J_bar  = J % at(solution)
  print '(a)', '(d) design: L = J + lambda_omega . r_omega + lambda_ic . r_ic'
  print '(a, es24.16)', 'nu found                = ', nu_bar % value(1, 1)
  print '(a, es24.16)', 'J                       = ', J_bar % value(1, 1)
  if (order >= 1) print '(a, es24.16)', 'dJ / dnu                = ', J_bar % value(1, 2)
  eta = L_bar % disagreement(solution)
  print '(a, es24.16)', 'estimate of J - J_bar   = ', sum(eta)

  !===================================================================!
  ! (e) DESIGN ON THE EQUIDISTRIBUTION. L of (d), the instants of
  ! manifold_bar moved, their count kept, until the error of the
  ! design found is the same over every step
  !===================================================================!

  L_bar = discrete_image(L, placement=equidistribution())
  call L_bar % minimize(estimate, solution)

  ! the design found, and the instants placed with the contribution of
  ! each to the estimate of J - J_bar, every tenth
  u_bar  = solution % fields(['u'])
  nu_bar = solution % fields(['nu'])
  J_bar  = J % at(solution)
  eta    = L_bar % disagreement(solution)
  print '(a)', '(e) design on the equidistribution: L of (d)'
  do k = 1, num_instants
     if (mod(k - 1, max(1, (num_instants - 1) / 10)) == 0) then
        print '(a, f8.4, a, es12.4)', 't = ', u_bar % instant(k), '   eta = ', eta(k)
     end if
  end do
  print '(a, es24.16)', 'nu found                = ', nu_bar % value(1, 1)
  print '(a, es24.16)', 'J                       = ', J_bar % value(1, 1)
  print '(a, es24.16)', 'estimate of J - J_bar   = ', sum(eta)

  !===================================================================!
  ! (f) DESIGN TO A TOLERANCE. L of (d), the count of the instants free:
  ! as many, equidistributed, as place the estimate of the error of the
  ! design found within 1e-3 of the design
  !===================================================================!

  L_bar = discrete_image(L, placement=equidistribution(tolerance=1.0e-3_dp))
  call L_bar % minimize(estimate, solution)

  nu_bar = solution % fields(['nu'])
  J_bar  = J % at(solution)
  print '(a)', '(f) design to a tolerance: L of (d)'
  print '(a, i24)', 'instants                = ', solution % num_points()
  print '(a, es24.16)', 'nu found                = ', nu_bar % value(1, 1)
  print '(a, es24.16)', 'J                       = ', J_bar % value(1, 1)

contains

  ! the discrete image of a Lagrangian upon manifold_bar: along t the
  ! chain dirk(2), bdf(2), adams(2), one block per instant, each block
  ! reading the instants before it; along nu the derivatives of the
  ! laws exact; with a placement, the instants placed by it
  function discrete_image(L, placement) result(L_bar)
    type(continuous_residual), intent(in) :: L
    type(equidistribution)   , intent(in), optional :: placement
    type(discrete_residual) :: L_bar
    integer :: k
    L_bar = L % discretize(manifold_bar, time=chain([(dirk(2), k = 1, 3), (bdf(2), k = 4, 7), &
         &                                           (adams(2), k = 8, num_instants)], from=[(k, k = 1, num_instants)]), &
         &                 placement=placement)
  end function discrete_image

  ! the k-th number of the command line, or the default without one
  real(dp) function argument(k, default)
    integer , intent(in) :: k
    real(dp), intent(in) :: default
    character(len=32) :: item
    argument = default
    if (command_argument_count() >= k) then
       call get_command_argument(k, item)
       read(item, *) argument
    end if
  end function argument

end program van_der_pol
