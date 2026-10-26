!=====================================================================!
! THE VAN DER POL OSCILLATOR: MOTION, ENERGY AND ITS DESIGN
!
! Consider changing the parameter nu of the Van der Pol oscillator.
! It changes the nonlinear damping and hence the motion. To choose nu,
! we need to compare motions over the interval of interest, [0, T].
! We shall account for both the displacement u and its rate of change
! u_t throughout that interval by integrating their quadratic energy:
!
!      J = integral over [0, T] of (u^2 + u_t^2) / 2.
!
! We can then ask which nu reduces J. We shall first determine the
! motion, then how J changes with nu, and then the value of nu that
! makes J stationary. A computed choice also raises two questions:
! how accurately have we determined it, and how much work is needed?
! Problems (a)-(h) develop these questions through one construction.
!
! The motion depends on time and on the parameter we wish to vary.
! Its domain is [0, T] x {nu}, represented by manifold, with coordinate
! functions t and nu. The unknown field u(t, nu) must satisfy
!
!      u_tt - nu (1 - u^2) u_t + u = 0,
!      u = u_0,  u_t = u_t_0  at t = 0.
!
! We state these required relations as laws whose residuals vanish.
! The motion equation gives r_omega on manifold, including its initial
! face. The two conditions fixing the initial data give r_ic on t = 0.
! To examine a prescribed nu before choosing it, we also require
! r_s = nu - nu_0 = 0 on the design factor {nu}. To hold the positions
! of the instants fixed, we require r_t = t_map - t = 0 on the time
! factor. Here t_map gives those positions. We call these conditions
! fixing the data, the design and the time map pins; removing one
! allows us to ask what determines the quantity it held fixed.
!
! To obtain the required relations by stationarity, we pair each
! residual with a multiplier of the same index. Once J is included,
! these multipliers will also let us determine how the energy responds
! to changes in the relations and their data. Each multiplier has one
! component per equation and depends on its domain's coordinates:
! lambda_omega(t, nu), lambda_ic(nu), lambda_s and lambda_t(t).
! Thus lambda_ic has one component per pinned field, lambda_s is one
! number, and lambda_t has one value per instant after discretization.
!
! (a) We begin by asking for the motion at the prescribed data and nu.
! For this question the pairings alone form the Lagrangian:
!
!      L = lambda_omega . r_omega + lambda_ic . r_ic
!        + lambda_s . r_s + lambda_t . r_t.
!
! Stationarity in the multipliers recovers the oscillator and its
! conditions. We can now choose how to represent this same statement
! numerically. Along t we use equally spaced instants, fixed by r_t,
! and a chain: dirk(2) over the first three instants, bdf(2) over the
! next four, and adams(2) to the end. Along nu we use an expansion at
! nu_0 of the requested order. These choices form manifold_bar and
! the discrete image L_bar without changing the continuous statement.
! The initial estimate is u = u_0 at every instant for every problem.
!
! The forward solve determines u_bar and, through the same expansion,
! its derivatives along nu up to that order. With no J in L, the
! stationarity equation in u has no source. The reverse sweep still
! runs because r_omega is paired, but all its multipliers are zero.
! The recorded result at eleven instants is max |lambda_omega| = 0
! and u(T) = 0.29118558490570512, identical to (b) in every digit.
!
! To solve L_bar, we pair it with a discrete minimizer M_bar from
! operation_minimization. It applies Newton's method to each block of
! instants with tolerance 1e-13. A Newton step first eliminates rows
! that fix an unknown from earlier unknowns through the history, the
! time scheme or the data; these are the tying rows. A dense LU solve
! determines the remaining unknowns. We retain M_bar as the questions
! grow: it solves every problem below, the placement rounds and the
! enrichment used to estimate error.
!
! (b) The motion now lets us evaluate J. To determine its response to
! changes in the motion equation and initial data, add J to L of (a).
! Stationarity in u then supplies an equation for that response, the
! adjoint equation. A reverse sweep determines lambda_omega(t), the
! sensitivity of J to a source added to the oscillator at t.
! Stationarity in the pinned initial values determines
! lambda_ic = -(dJ/du_0, dJ/du_t_0).
!
! The energy must be integrated consistently with the computed motion.
! A DIRK step uses its tableau's weights at its stages; a multistep
! step uses the interpolatory rule over its instants. At expansion
! order zero we solve only at nu_0. Increasing the order retains the
! same state, multipliers and energy whilst adding their derivatives.
! At eleven instants the recorded J = 3.1356779007580480 and
! lambda_ic = (-1.3772000188879279, -0.10639333936774159).
! As the count increases, u(2) = 0.32331 at both 401 and 1601 instants.
!
! Before using J to choose nu, we need to assess its numerical error.
! We compare the solution with a chain one order higher on the same
! instants: dirk(3), bdf(3), adams(3). Write Q for the solution and +
! for this enriched chain. To form P Q, retain Q at the instants and
! solve each dirk(3) step's stages from the instant behind that step.
! The change in quadrature and the residual of the enriched equations,
! weighted by their adjoint mu+, give the disagreement
!
!      J+(P Q) - J(Q) + mu+ . R+(P Q).
!
! Its contributions, one per instant, estimate J - J_bar. This gives
! both an estimate of the total error and its distribution along time,
! which we shall use when the positions of the instants become free.
!
! The recorded reference is J = 3.0969362409469 from the chain of
! order three at 2561 instants; 5121 instants agree to 2.8e-12.
! At 11, 21, 41, 81, 161 and 321 instants, the disagreement divided
! by the error relative to that reference is 1.259, 1.102, 1.032,
! 1.0099, 1.0029 and 1.0008. Divided instead by the difference from
! the order-three solve on the same instants, it is 1.0215, 1.0039,
! 1.00068, 1.00013, 1.00003 and 1.00001. Grading the instants toward
! either end gives ratios to the reference error of 0.93 to 1.30 at
! 21 instants and 0.988 to 1.025 at 161.
!
! (c) Suppose we now ask how the energy changes with nu. No term needs
! to change in L of (b): stationarity in the pinned design already
! gives lambda_s = -dJ/dnu. The expansion of the state gives the same
! energy derivatives by a second route. We can therefore audit
! d^k J / dnu^k against the derivative of order k - 1 of -lambda_s,
! obtained from the adjoint expansion in the same solve.
!
! At eleven instants and order three, the recorded differences are
! 5.6e-16, 4.4e-16 and 1.1e-15. At 21 instants and order three, the
! state and energy derivatives also agree with central differences
! at nu_0 = 1 +- 0.01, +- 0.02 to their truncation error. The
! expansion and central differences give, respectively,
!
!      du(T)/dnu = 1.2502856 against 1.2502527,
!      d^2u/dnu^2 = -0.598952 against -0.598900,
!      d^3u/dnu^3 = -1.9712 against -1.9710;
!      dJ/dnu = -0.6565620 against -0.6565450,
!      d^2J/dnu^2 = 1.548406 against 1.548272,
!      d^3J/dnu^3 = 1.0211 against 1.0215.
!
! The audit in test/graph-scheme/van_der_pol_design_audit.f90 extends
! this comparison to two design coordinates, the mixed derivative of
! J and the symmetry of its Hessian.
!
! The same argument applies to the positions at which we compute the
! motion. Stationarity in the pinned time map gives lambda_t =
! -dJ_bar/dt_k. Moving an interior instant changes J_bar whilst the
! continuous J remains fixed, so dJ_bar/dt_k is minus the derivative
! of the error J - J_bar. Moving an endpoint changes the interval and
! hence J as well. These multipliers will let us examine the error
! through changes of the instants in (g).
!
! We can check this interpretation because the oscillator has no
! explicit dependence on t. Shifting all instants together leaves
! J_bar unchanged, so the sum of lambda_t should vanish. At order
! three its recorded sum is -1.1e-15 at eleven instants and -4.9e-15
! at 21. At 21 instants the interior values agree with -dJ_bar/dt_k
! from the state expansion to 4e-12 on dirk(2) and 1.7e-11 on the
! program's chain. At the two endpoints of the chain, with the
! interval moving, central differences give agreement to 6e-12
! and 6e-11.
!
! (d) Knowing the energy derivative allows us to choose nu. Remove
! lambda_s . r_s from L of (b)-(c), keeping the initial data and the
! instants fixed. The stationarity condition in nu now determines nu
! where dJ/dnu = 0, instead of determining the reaction to its pin.
! Starting at nu_0, Newton's iteration uses one forward and one
! reverse sweep per iterate: lambda_s supplies the gradient and its
! order-one expansion supplies the Hessian. At order zero, a BFGS
! update supplies the approximation needed for the design iteration.
!
! At eleven instants the recorded |dJ/dnu| falls quadratically as
! 7.07e-1, 3.84e-2, 1.82e-3, 4.90e-6, 3.58e-11. The result is
! nu = 1.5269121988 and J = 2.9656219, compared with 3.1356779 at
! nu_0. BFGS at order zero reaches the same nu within 4.1e-11 in ten
! iterations. Thus the derivative supplies a design choice, but that
! choice still depends on the instants used to compute it.
!
! At 21, 41, 81, 161 and 321 instants, the recorded designs are
! nu = 1.49414469, 1.48607327, 1.48427255, 1.48387461 and 1.48378216.
! Their successive differences are 8.1e-3, 1.8e-3, 4.0e-4 and 9.2e-5,
! consistent with the chain's second order. At each design found, we
! can also assess J using the disagreement from (b). Relative to J
! at that same design from the order-three chain at 2561 instants,
! the estimate-to-error ratios are 1.259, 1.136, 1.066 and 1.028 at
! 21, 41, 81 and 161 instants.
!
! (e) Let us now improve the accuracy of the chosen nu by moving the
! instants whilst keeping their count. We remove lambda_t . r_t from
! L of (d), leaving J + lambda_omega . r_omega + lambda_ic . r_ic.
! The required motion and initial data remain the same. The new
! freedom belongs to the numerical representation: the instants of
! manifold_bar are placed by equidistribution in the discrete image
! of L. They are not coordinates of L's physical design factor; see
! doc/unified-residuals.md, section 16. Their own design problem will
! become explicit in (g).
!
! Which error should determine their placement? For a fixed nu_0,
! distributing the contributions eta to the error of J lowers that
! error by 33, 30, 25 and 20 per cent at 21, 41, 81 and 161 instants.
! With nu free, distributing eta lowers the error of J at the design
! found three- to tenfold, yet raises the error of nu by 13 to 34 per
! cent at 41 to 161 instants. The design is determined by the gradient
! of J, so it is the error in that gradient that we need here.
!
! Write nu_h for the computed design and nu* for its reference value.
! With E = J - J_bar and H the derivative of the design gradient,
!
!      nu* - nu_h ~ -H^(-1) E'(nu_h).
!
! The expansion of the disagreement along nu supplies each step's
! derivative of eta. The reverse sweep's expansion supplies H, using
! one solve at nu_h. Thus each step contributes to an estimate of the
! design error. Its magnitude to the power 1/3, divided by the step
! length, gives a density whose quantiles place the instants, seeking
! equal contributions across the steps. We restrict each step to at
! most 1 + sqrt(2) times its predecessor, using Grigorieff's bound
! for bdf(2).
!
! Each round of minimize solves (d), estimates the design error and
! places the instants again. A round that fails to lower the sum of
! the contribution magnitudes is rejected and the preceding round
! returned. A round that lowers it by less than the uncertainty of
! the estimate is accepted and returned. That uncertainty is the
! distance from the design-error estimate obtained using the
! order-three chain on the same instants.
!
! Against nu = 1.4837523856 from the order-three chain at 2561
! instants, the recorded errors at 11, 21, 41 and 81 equally spaced
! instants are 4.316e-2, 1.039e-2, 2.321e-3 and 5.202e-4. With
! equidistribution they become 3.398e-2, 7.346e-3, 1.606e-3 and
! 3.935e-4: reductions of 21, 29, 31 and 24 per cent. The estimates
! of nu* - nu_h are 0.942, 0.948, 0.969 and 0.982 times those errors.
!
! (f) A fixed count permits a better placement, but does not specify
! the accuracy of nu. Suppose instead that we require its estimated
! error to be within 1e-3 times the norm of the design. We retain L
! of (e) and let equidistribution choose the count as well as the
! positions. At each count, the rounds of (e) first settle placement.
! If the tolerance is unmet, the estimated fall with step count to
! the power -r predicts the larger count needed. We use r equal to
! two, the chain's order, until two counts provide an observed order.
! The chain's last family continues to the new count, and the previous
! solution interpolated along time supplies the next initial estimate.
! If the estimate fails to decrease from the preceding count, or
! falls within its own uncertainty, the rounds end with the tolerance
! unmet.
!
! Starting from eleven instants, the recorded 1e-3 tolerance is met
! at 47, with nu = 1.4849560023 and error 1.204e-3; equally spaced
! instants attain that error at 57. At 1e-4 the count is 147 with
! error 1.171e-4, compared with 165 equally spaced instants. The same
! rule can answer our original energy question at fixed nu_0: for J,
! the tolerance is relative to sum |w f|, with quadrature weights w
! and integrand values f. It meets 1e-3 at 33 instants with error
! 2.98e-3, compared with 34 equally spaced, and 1e-4 at 100 with
! error 2.15e-4, compared with 116.
!
! (g) The energy tolerance can also guide the choice of time scheme.
! Return to nu_0 by adding lambda_s . r_s to L of (f); this is L of
! (b) without lambda_t . r_t. Now choose each step's family together
! with the instants to obtain J within 1e-6 of sum |w f|. The choices
! form pool: dirk(2), dirk(3), dirk(4), adams(2), adams(3), adams(4),
! bdf(2) and bdf(3).
!
! We have another design problem of the same kind. Its variables s
! are the instants and each step's weights on the families in pool.
! Its objective is the work, measured by the solves of the steps.
! Its constraint requires sum |e_k| to meet the tolerance, where e_k
! is step k's contribution to the error of J. Its Lagrangian is the
! work plus lambda times the error sum less the tolerance. A family
! is excluded wherever its history would precede the first instant;
! bdf(2) is excluded where a step exceeds 1 + sqrt(2) times the
! preceding step. These conditions specify which choices are allowed.
!
! To compare those choices, we need the error of each family at each
! step. Recall from (c) that J is independent of interior instant
! positions. More generally J does not depend on s, so a derivative
! of J_h along s is minus that of J - J_h; here J_h denotes the
! discrete energy. We first solve dirk(4), the one-step family of
! highest order in pool, with r_t restored to pin the instants.
! Its multipliers give -dJ_h/dt_k. Differentiating them along three
! directions, each moving every third interior instant, recovers the
! band of the Hessian of J_h. At leading order the error of an
! interior step of length h_k and family order p is
!
!      h_k^2 (d^2 J_h / dt_(k-1) dt_k) / (p (p + 1)).
!
! At the two steps touching the fixed endpoints, the diagonal and
! neighbouring mixed curvature give the corresponding contributions.
! For another family we subtract from the base family's error the
! change in J_h if that step used the other family. The disagreement
! with that family's chain supplies this change. Thus the multipliers
! and disagreement already introduced supply the information needed
! to choose both the positions and the families.
!
! The leading-order stationarity conditions choose, at each instant,
! the family with the least (1 + 1/p) w n, its work per unit time,
! where w is the solves per step, p the order, a the error density,
! and n = (lambda p a / w)^(1/(p + 1)) the steps per unit time.
! The multiplier lambda is chosen to meet the error tolerance. Each
! round solves the base family on the current instants, then places
! them again and chooses their families using these conditions.
!
! From eleven instants, the recorded 1e-6 tolerance gives 52 instants
! and 55 solves. The first two steps use dirk(4), since no multistep
! family in pool yet has its history. The rest use adams(4), except
! for one step of adams(3). Against J from the order-three chain at
! 2561 instants, the error is -2.430e-6 and its estimate -2.295e-6.
! The program's chain requires 887 solves on equally spaced instants
! to meet that tolerance. The best single family in pool, adams(3)
! after dirk(4), requires 53. At tolerances 1e-7, 1e-5 and 1e-4 the
! rounds return 97, 31 and 17 solves; the best family alone takes
! 96, 30 and 17, and the program's chain takes 2787, 286 and 95.
!
! Choosing the representation also costs work. At 1e-4, 1e-5, 1e-6
! and 1e-7 the recorded times for the rounds are 0.34, 0.74, 1.29
! and 2.25 s. One solve of the program's chain at its required count
! takes 0.03, 0.10, 0.30 and 0.96 s; adams(3) takes at most 0.03 s.
! Each round's cost grows with the count because the curvature uses
! only three directions of the pins' multipliers. Computing it with
! 2N - 5 order-two state expansions took 0.81, 2.70, 8.0 and 24.8 s
! at those same tolerances.
!
! (h) Suppose the available work is prescribed instead of the error.
! Keep the physical Lagrangian of (g) and its pool, and allow 55
! solves, the work (g) used. In the Lagrangian for s, make the error
! the objective and the work its constraint:
!
!      L_s = sum |e_k| + mu (work - 55).
!
! With mu = 1/lambda its leading-order stationarity conditions are
! those of (g). We can therefore use the same rounds to return the
! design of least estimated error within 55 solves. Our original
! question about energy has become a question about how accurately
! we can determine it with the work available.
!
! From eleven instants, the recorded result has 51 instants and error
! -2.45e-6 relative to the order-three chain at 2561 instants. In
! comparison, (g) at 1e-6 uses 55 solves with error -2.43e-6. At 31
! and 97 solves, (h) gives errors -2.33e-5 and -2.24e-7, whilst (g)
! at 1e-5 and 1e-7 uses 31 and 97 solves with errors -2.52e-5 and
! -2.27e-7. Least work within a tolerance and least error at that
! work are thus two readings of the same design construction.
!
! The program prints the Lagrangian of (b) and (c) before solving it:
! J and each law paired with its multiplier, each equation as a
! formula, then the whole as a tree. It then prints tree(M_bar),
! showing the march over the blocks and the minimizer of a block.
! The command line sets the instant count, u_0, u_t_0 and expansion
! order, with defaults 11, 2, 0 and 1.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

program van_der_pol

  use util_precision        , only : dp
  use util_verbosity        , only : set_verbosity
  use operation_manifold    , only : continuous_manifold, discrete_manifold, interval, points, parameter
  use operation_expansion   , only : expansion
  use operation_field       , only : continuous_field, discrete_field, integral, &
       &                             operator(+), operator(-), operator(*), operator(/), operator(**)
  use operation_lagrangian  , only : continuous_lagrangian, discrete_lagrangian, operator(+), operator(*), &
       &                             equidistribution
  use operation_family      , only : family, dirk, bdf, adams, chain
  use view_residual         , only : formula, tree
  use operation_newton      , only : newton
  use operation_elimination , only : elimination
  use operation_dense_direct, only : dense_direct
  use operation_minimization, only : discrete_minimizer

  implicit none

  real(dp), parameter :: T_final = 2.0_dp, nu_0 = 1.0_dp
  real(dp) :: u_0, u_t_0
  integer  :: num_instants, order

  type(continuous_manifold)   :: manifold
  type(discrete_manifold)     :: manifold_bar
  type(continuous_field)      :: t, t_map, nu, u, J, lambda_omega, lambda_ic, lambda_s, lambda_t, rest
  type(continuous_lagrangian) :: r_omega, r_ic, r_s, r_t, L
  type(discrete_lagrangian)   :: L_bar
  type(family), allocatable   :: pool(:)
  type(discrete_field)        :: estimate, solution, u_bar, lambda_omega_bar, lambda_ic_bar, lambda_s_bar, lambda_t_bar, &
       &                         J_bar, nu_bar
  real(dp), allocatable       :: eta(:)
  type(discrete_minimizer)    :: M_bar
  integer                     :: k

  num_instants = nint(argument(1, 11.0_dp))
  u_0          = argument(2, 2.0_dp)
  u_t_0        = argument(3, 0.0_dp)
  order        = nint(argument(4, 1.0_dp))

  ! We need the motion over time and its dependence on nu. Define
  ! that domain first, then its coordinates and the unknown u.
  manifold = continuous_manifold(time=interval(0.0_dp, T_final), design=parameter('nu', nu_0))
  t        = manifold % coordinate('t')
  nu       = manifold % coordinate('nu')
  u        = manifold % unknown('u', [t, nu])

  ! Require the motion equation and the prescribed initial data.
  ! The conditions on nu and t_map hold the design and instants fixed
  ! until we choose to release them.
  t_map   = manifold % time_map()
  r_omega = continuous_lagrangian(manifold, [u % derivative([t, t]) - nu*(1.0_dp - u*u)*u % derivative([t]) + u])
  r_ic    = continuous_lagrangian(manifold % boundary(t, 0.0_dp), [u - u_0, u % derivative([t]) - u_t_0])
  r_s     = continuous_lagrangian(manifold % design(), [nu - nu_0])
  r_t     = continuous_lagrangian(manifold % time(), [t_map - t])

  ! Integrate the energy over time to compare the motions through J.
  J = integral((u*u + u % derivative([t])**2)/2.0_dp, over=manifold % time())

  ! Pair each relation with a multiplier so that stationarity can
  ! recover the motion and, with J included, its energy sensitivities.
  lambda_omega = r_omega % multiplier('lambda_omega', [t, nu])
  lambda_ic    = r_ic % multiplier('lambda_ic', [nu])
  lambda_s     = r_s % multiplier('lambda_s')
  lambda_t     = r_t % multiplier('lambda_t', [t])

  ! Represent time by instants and dependence on nu by an expansion.
  ! Start every problem from the constant estimate u = u_0.
  manifold_bar = manifold % discretize(time=points(num_instants), design=expansion(order=order))
  rest         = continuous_field(manifold, [u_0])
  estimate     = rest % discretize(manifold_bar)

  ! Let image 1 report convergence so we can follow each solve.
  if (this_image() == 1) call set_verbosity(1)

  !===================================================================!
  ! (a) PHYSICS. Determine the motion with data, nu and instants fixed.
  ! Without J, the multiplier equations have no source.
  !===================================================================!

  L     = lambda_omega*r_omega + lambda_ic*r_ic + lambda_s*r_s + lambda_t*r_t
  L_bar = discrete_image(L)

  ! Use Newton with elimination and a dense solve for every block.
  ! Retain this M_bar for the questions that follow.
  M_bar = discrete_minimizer(L_bar, newton(elimination(dense_direct()), tolerance=1.0e-13_dp))
  call L_bar % minimize(estimate, solution, by=M_bar)

  u_bar            = solution % fields(['u'])
  lambda_omega_bar = solution % fields(['lambda_omega'])
  lambda_ic_bar    = solution % fields(['lambda_ic'])
  lambda_s_bar     = solution % fields(['lambda_s'])
  print '(a)', '(a) physics: L = lambda_omega . r_omega + lambda_ic . r_ic + lambda_s . r_s + lambda_t . r_t'
  print '(a, es24.16)', 'u(T)                    = ', u_bar % value(num_instants, 1)
  print '(a, es24.16)', 'max |lambda_omega|      = ', &
       & maxval([(abs(lambda_omega_bar % value(k, 1)), k = 1, num_instants)])
  print '(a, 3es24.16)', 'lambda_ic, lambda_s     = ', lambda_ic_bar % value(1, 1), lambda_ic_bar % value(1, 2), &
       & lambda_s_bar % value(1, 1)

  !===================================================================!
  ! (b) ADJOINT and (c) SENSITIVITY. Add J to ask how energy responds
  ! to the motion equation, initial data, nu and instant positions.
  !===================================================================!

  L     = J + L
  L_bar = discrete_image(L)

  ! Show the objective and its paired relations before solving them,
  ! first as formulas and then as a tree.
  print '(a)', 'the Lagrangian of (b) and (c):'
  print '(a)', formula(L)
  print '(a)', tree(L)

  ! Show how M_bar marches over blocks and solves each block.
  print '(a)', 'the minimizer of (b) and (c):'
  print '(a)', tree(M_bar)

  call L_bar % minimize(estimate, solution, by=M_bar)

  ! (b) Sample the motion and its adjoint across tenths of the interval.
  ! Report the initial-data multipliers as energy sensitivities.
  u_bar            = solution % fields(['u'])
  lambda_omega_bar = solution % fields(['lambda_omega'])
  lambda_ic_bar    = solution % fields(['lambda_ic'])
  print '(a)', '(b) adjoint: L = J + lambda_omega . r_omega + lambda_ic . r_ic + lambda_s . r_s + lambda_t . r_t'
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

  ! (b) Estimate the error of J from the disagreement at each instant.
  eta = L_bar % disagreement(solution, by=M_bar)
  print '(a, es24.16)', 'estimate of J - J_bar    = ', sum(eta)

  ! (c) Read the energy sensitivity from the design pin. Compare each
  ! state-expansion derivative of J with the derivative of order
  ! k - 1 of -lambda_s to audit the same quantity by the adjoint.
  lambda_s_bar = solution % fields(['lambda_s'])
  print '(a)', '(c) sensitivity: the same solve'
  print '(a, es24.16)', 'lambda_s = -dJ / dnu     = ', lambda_s_bar % value(1, 1)
  if (order >= 1) print '(a)', 'k, d^k J / dnu^k by the expansion of the state, by that of the adjoint, their difference'
  do k = 1, order
     if (k > 1) lambda_s_bar = lambda_s_bar % derivative([nu])
     print '(i3, 2es24.16, es12.2)', k, J_bar % value(1, 1 + k), -lambda_s_bar % value(1, 1), &
          & J_bar % value(1, 1 + k) + lambda_s_bar % value(1, 1)
  end do

  ! (c) Differentiate u(T) to show how the final motion depends on nu.
  do k = 1, order
     u_bar = u_bar % derivative([nu])
     print '(a, i0, a, i0, a, es24.16)', 'd^', k, ' u(T) / dnu^', k, '        = ', u_bar % value(num_instants, 1)
  end do

  ! (c) Sample -dJ_bar/dt_k across tenths of the interval. Check its
  ! sum: shifting every instant together leaves J_bar unchanged, so
  ! the time-map multipliers should sum to zero.
  lambda_t_bar = solution % fields(['lambda_t'])
  do k = 1, num_instants
     if (mod(k - 1, max(1, (num_instants - 1) / 10)) == 0) then
        print '(a, f8.4, a, es16.8)', 't = ', manifold_bar % instant(k), '   lambda_t = -dJ_bar / dt_k = ', &
             & lambda_t_bar % value(k, 1)
     end if
  end do
  print '(a, es24.16)', 'sum of lambda_t         = ', sum([(lambda_t_bar % value(k, 1), k = 1, num_instants)])

  !===================================================================!
  ! (d) DESIGN. Remove r_s so the energy derivative determines nu
  ! through dJ/dnu = 0, with the data and instants still fixed.
  !===================================================================!

  L     = J + lambda_omega*r_omega + lambda_ic*r_ic + lambda_t*r_t
  L_bar = discrete_image(L)
  call L_bar % minimize(estimate, solution, by=M_bar)

  nu_bar = solution % fields(['nu'])
  J_bar  = J % at(solution)
  print '(a)', '(d) design: L = J + lambda_omega . r_omega + lambda_ic . r_ic + lambda_t . r_t'
  print '(a, es24.16)', 'nu found                = ', nu_bar % value(1, 1)
  print '(a, es24.16)', 'J                       = ', J_bar % value(1, 1)
  if (order >= 1) print '(a, es24.16)', 'dJ / dnu                = ', J_bar % value(1, 2)
  eta = L_bar % disagreement(solution, by=M_bar)
  print '(a, es24.16)', 'estimate of J - J_bar   = ', sum(eta)

  !===================================================================!
  ! (e) DESIGN ON THE EQUIDISTRIBUTION. Remove r_t to place the same
  ! number of instants by the estimated design error, seeking equal
  ! contributions across the steps.
  !===================================================================!

  L     = J + lambda_omega*r_omega + lambda_ic*r_ic
  L_bar = discrete_image(L, placement=equidistribution())
  call L_bar % minimize(estimate, solution, by=M_bar)

  ! Report nu and J after placement. Sample the placed instants and
  ! their contributions to the estimate of J - J_bar across tenths.
  u_bar  = solution % fields(['u'])
  nu_bar = solution % fields(['nu'])
  J_bar  = J % at(solution)
  eta    = L_bar % disagreement(solution, by=M_bar)
  print '(a)', '(e) design on the equidistribution: L = J + lambda_omega . r_omega + lambda_ic . r_ic'
  do k = 1, num_instants
     if (mod(k - 1, max(1, (num_instants - 1) / 10)) == 0) then
        print '(a, f8.4, a, es12.4)', 't = ', u_bar % instant(k), '   eta = ', eta(k)
     end if
  end do
  print '(a, es24.16)', 'nu found                = ', nu_bar % value(1, 1)
  print '(a, es24.16)', 'J                       = ', J_bar % value(1, 1)
  print '(a, es24.16)', 'estimate of J - J_bar   = ', sum(eta)

  !===================================================================!
  ! (f) DESIGN TO A TOLERANCE. Keep L of (e), and let placement choose
  ! the count needed for an estimated error within 1e-3 of the design.
  !===================================================================!

  L_bar = discrete_image(L, placement=equidistribution(tolerance=1.0e-3_dp))
  call L_bar % minimize(estimate, solution, by=M_bar)

  nu_bar = solution % fields(['nu'])
  J_bar  = J % at(solution)
  print '(a)', '(f) design to a tolerance: L of (e)'
  print '(a, i24)', 'instants                = ', solution % num_points()
  print '(a, es24.16)', 'nu found                = ', nu_bar % value(1, 1)
  print '(a, es24.16)', 'J                       = ', J_bar % value(1, 1)

  !===================================================================!
  ! (g) J TO A TOLERANCE, EACH STEP'S FAMILY FROM A POOL.
  ! Fix nu again. Choose instants and families together to seek the
  ! least work for an estimated error within 1e-6 of sum |w f|.
  !===================================================================!

  L     = J + lambda_omega*r_omega + lambda_ic*r_ic + lambda_s*r_s
  pool  = [dirk(2), dirk(3), dirk(4), adams(2), adams(3), adams(4), bdf(2), bdf(3)]
  L_bar = discrete_image(L, placement=equidistribution(tolerance=1.0e-6_dp, pool=pool))
  call L_bar % minimize(estimate, solution, by=M_bar)

  J_bar = J % at(solution)
  print '(a)', '(g) J to a tolerance, each step''s family from a pool: &
       &L = J + lambda_omega . r_omega + lambda_ic . r_ic + lambda_s . r_s'
  print '(a, i24)', 'instants                = ', solution % num_points()
  print '(a, es24.16)', 'J                       = ', J_bar % value(1, 1)

  !===================================================================!
  ! (h) J AT A WORK. Keep L and the pool of (g). Prescribe its work,
  ! 55 solves, and seek the least estimated energy error within it.
  !===================================================================!

  L_bar = discrete_image(L, placement=equidistribution(work=55, pool=pool))
  call L_bar % minimize(estimate, solution, by=M_bar)

  J_bar = J % at(solution)
  print '(a)', '(h) J at a work of 55 solves, each step''s family from the pool of (g): L of (g)'
  print '(a, i24)', 'instants                = ', solution % num_points()
  print '(a, es24.16)', 'J                       = ', J_bar % value(1, 1)

contains

  ! Represent L on manifold_bar without changing its statement.
  ! The chain dirk(2), bdf(2), adams(2) has one block per instant and
  ! uses earlier instants as history. Along nu the law derivatives
  ! are exact. If supplied, placement determines the instants.
  function discrete_image(L, placement) result(L_bar)
    type(continuous_lagrangian), intent(in) :: L
    type(equidistribution)   , intent(in), optional :: placement
    type(discrete_lagrangian) :: L_bar
    integer :: k
    L_bar = L % discretize(manifold_bar, time=chain([(dirk(2), k = 1, 3), (bdf(2), k = 4, 7), &
         &                                           (adams(2), k = 8, num_instants)], from=[(k, k = 1, num_instants)]), &
         &                 placement=placement)
  end function discrete_image

  ! Read the k-th command-line value, retaining its default if absent.
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
