!=====================================================================!
! THE TEMPLATE of a problem of a flexible mechanism. A case supplies
! the assembly and its data under the names imported below - HART II
! in hart2_case.f90 - and this text is the same for every case.
!
! THE MECHANISM, a graph: members [ ] and flexible members ( ), joints
! on the edges, # the ground.
!
!                        controls u(t)
!                             |
!                             v
!   # --support-- [ member ] <==driver== #
!                     |
!                   joint
!                     |
!                 [ member ] --joint-- [ member ] --support-- #
!                     |
!                 attachment: the incidence, and the map of a face of D_m
!                     |
!              ( member, on its domain D_m )
!
! THE EQUATIONS. A rigid member's unknowns are its pose, fields of
! time; a flexible member's, its displacements on [0, T] x D_m in the
! axes its attachment gives it. With q the unknowns and u(t) the
! controls:
!
!   density  sum over the members of (T - V) + lambda . g(q) + mu . c(q, u(t)),
!            the flexible members' terms integrals over their domains
!   joints   g(q) = 0, their reactions lambda; a support is a joint
!   drivers  c(q, u(t)) = 0 on the actuated coordinates alone - a(q) - A u
!            when linear - their efforts e = -(dc/du)^T mu, A^T mu when linear
!   energy   dE/dt = e^T u', E = T + V; over a window the work of the
!            efforts less the energy gained is the energy balance's defect
!   data     at t = 0, q' = 0 and q'' = 0, and q an equilibrium fixed by
!            the reference phase condition K^T M (q - q_0) = 0: K the
!            kernel of the equilibrium's jacobian at the reference q_0 -
!            the coordinates neither the joints nor the loads determine -
!            M the mass at q_0, the metric of the projection. The
!            equilibrium's jacobian augmented by the condition is to be
!            nonsingular at q, which the audit checks. It is not the
!            equilibrium nearest q_0, which a minimization would select
!
! THE STATEMENT, in the order in which each part reads the ones before:
!
!   THE ASSEMBLY        once: the members, their domains D_m and attachments,
!                         |  the joints, the drivers, the design h;
!                         |  Sigma = (D_1 + ... + D_M) x {h}
!   THE DESIGN          on Sigma: L_design = mass + <mu_bounds, d>
!                         |
!   THE SCENARIOS       each j on Omega_j = [0, T] x Sigma, the assembly shared:
!                         |  L_j = <adjoint, r> + <lambda, g> + <nu, b> + <mu, c>,
!                         |  c the stress bound on the observation window
!   THE DISCRETE IMAGE  Sigma_h by the members' models; the instants between the
!                         |  case's breakpoints; the march
!   THE ZERO            of L = L_design + sum_j L_j
!                         |
!   THE AUDIT           the data at t = 0 and the augmented jacobian there; the
!                       conditions; the stress bound's largest point; the energy
!                       balance and its defect
!
! EACH FACT HAS ONE OWNER:
!   the members, their domains and attachments, the joints,     the assembly,
!   the drivers, the design-dependent properties                   built by the case
!   a member's theory, section, strength and stress               the member's model
!   the controls, the breakpoints, the bounds and the limit       the case
!   the state, the reactions and the efforts of a scenario        the mechanism
!   the objective, the constraints, the audit                     this text
!   the instants, the elements and the march                      the discrete image
!
! A CASE SUPPLIES: machine(), the assembly; controls(t, j), the
! controls of scenario j; bounds(h), the design's bounds; scenarios,
! their number; breakpoints, the instants the time grid keeps, from
! the start to the end; observation, the window of the stress bound and
! of the energy balance; limit, the bound of the ratio of the stress to
! the strength.
!
! NOT YET STATED: members of domains of their own. Sigma is the one
! domain the flexible members share, of one interval and its points,
! times the design; an assembly whose members differ stops the program.
!
! MEASURED on HART II at the step 8e-4, on one core: the start at rest
! finds the kernel of nine at the reference and the droop in fourteen
! iterations of Newton, where its residual stops at 3.2e-5, 5.7e-9 of
! its first: the reaction of the phase condition along the upper
! swashplate's turn. The march reached t = 0.157 s, 198 of its 1395
! blocks, in 55 minutes, some 1.8 s an iteration of Newton: 13
! iterations a block on average in the first 0.04 s, 6 from t = 0.08 s,
! seven in ten of them at the rounding of the residual. Past t = 0.36 s,
! in a run whose Newton stopped at that rounding, the blocks took 16 to
! 20 iterations, up to 38, most of them steps damped from the block's
! estimate. A scenario's march, projected from these rates, takes more
! than ten hours.
!=====================================================================!

program flexible_mechanism

  use util_precision     , only : dp
  use operation_manifold , only : continuous_manifold, discrete_manifold, window, interval, points
  use operation_expansion, only : expansion
  use operation_field    , only : continuous_field, discrete_field, maximum
  use operation_residual , only : continuous_residual, discrete_residual, minimize
  use operation_family   , only : dirk, bdf, adams, chain
  use multibody_dynamics , only : assembly, mechanism
  use hart2_case         , only : machine, controls, bounds, scenarios, breakpoints, observation, limit

  implicit none

  ! the weight of the smooth maximum
  real(dp), parameter :: weight = 1.0e4_dp

  type(assembly)            :: layout
  type(continuous_manifold) :: sigma, omega(scenarios)
  type(discrete_manifold)   :: sigma_h, omega_h(scenarios)
  type(points)              :: instants
  type(window)              :: observed(scenarios)
  type(continuous_field)    :: t, mass, stresses(scenarios)
  type(continuous_field), allocatable :: h(:)
  type(mechanism)           :: mechanisms(scenarios)
  type(continuous_residual) :: d, r, g, b, c, L_design, L(scenarios)
  type(discrete_residual)   :: L_design_h, L_h(scenarios)
  type(chain)               :: march
  type(discrete_field)      :: estimates(scenarios), solutions(scenarios), h_h, mu_h
  real(dp) :: step
  integer  :: num_instants, j, k

  ! the instants: uniform between the case's breakpoints, at most the
  ! step apart, every breakpoint an instant; the march takes three
  step = argument(1, 8.0e-4_dp)
  if (.not. step > 0.0_dp) error stop 'flexible_mechanism: the step is positive'
  instants     = points(breakpoints, step)
  num_instants = instants % count()
  if (num_instants < 3) error stop 'flexible_mechanism: the march takes three instants at least'

  ! THE ASSEMBLY, once; its manifold Sigma = (D_1 + ... + D_M) x {h}
  layout = machine()
  sigma  = layout % manifold()

  ! THE DESIGN on Sigma: the objective, the members' mass in kg - the
  ! rigid members' and each flexible member's integral over its domain -
  ! and the design's bounds
  h        = sigma % coordinate('h')
  mass     = layout % mass()
  d        = continuous_residual(sigma % design(), bounds(h), inequality=.true.)
  L_design = mass + d % multiplier('mu_bounds') * d

  do j = 1, scenarios

     ! THE SCENARIO on Omega_j = [0, T] x Sigma: the assembly in motion,
     ! driven by the scenario's controls, its state its own
     omega(j)      = continuous_manifold(time=interval(breakpoints(1), breakpoints(size(breakpoints))), base=sigma)
     t             = omega(j) % coordinate('t')
     mechanisms(j) = mechanism(omega(j), layout, controls(t, j))

     ! the stress bound, one definition for the constraint and the audit:
     ! the squared ratio of the stress to the strength over the flexible
     ! members, their sections' points and their domains, on the window
     stresses(j) = mechanisms(j) % stress()
     observed(j) = omega(j) % window(t, observation(1), observation(2))

     ! the residuals: on Omega_j the laws; on t = 0 the data at rest, an
     ! equilibrium under the reference phase condition; on each flexible
     ! member's face its attachment; on {h} the stress bound, at most the
     ! limit squared
     r = continuous_residual(omega(j), mechanisms(j) % laws())
     g = continuous_residual(omega(j) % boundary(t, 0.0_dp), mechanisms(j) % initial([1, 2]))
     b = mechanisms(j) % attachments()
     c = continuous_residual(omega(j) % design(), &
          & [maximum(stresses(j), over=observed(j), weight=weight) - limit**2], inequality=.true.)

     ! the scenario's Lagrangian, every residual paired with its multiplier
     L(j) = r % multiplier('adjoint') * r + g % multiplier('lambda') * g + b % multiplier('nu') * b &
          & + c % multiplier('mu') * c

  end do

  ! THE DISCRETE IMAGE: Sigma_h by the members' models, each domain by
  ! its elements at its resolution, along h the expansion of order zero;
  ! each Omega_j,h the instants times Sigma_h; along t dirk(2) at the
  ! first instant and for the startup, then bdf(2), then adams(2), one
  ! block per instant
  sigma_h    = layout % discretize(design=expansion(order=0))
  L_design_h = L_design % discretize(sigma_h)
  march      = chain([dirk(2), dirk(2), bdf(2), (adams(2), k = 4, num_instants)], from=[(k, k = 1, num_instants)])
  do j = 1, scenarios
     omega_h(j)   = omega(j) % discretize(time=instants, base=sigma_h)
     L_h(j)       = L(j) % discretize(omega_h(j), time=march)
     estimates(j) = mechanisms(j) % estimate(omega_h(j))
  end do

  ! THE ZERO of L = L_design + sum_j L_j, the images sharing the design
  call minimize(L_design_h, L_h, estimates, solutions)

  ! the design found, and in each scenario the multiplier of its stress
  ! bound
  h_h = solutions(1) % fields(['h'])
  print '(a, *(f9.5))', 'the design found: ', h_h % value(1, :)
  do j = 1, scenarios
     mu_h = solutions(j) % fields(['mu'])
     print '(a, i0, a, es16.8)', 'scenario ', j, ', mu = ', mu_h % value(1, 1)
  end do

  ! THE AUDIT in each scenario: the controls at rest at t = 0, the
  ! equilibrium there and its augmented jacobian nonsingular; the
  ! joints' conditions at every instant; the
  ! stress bound's largest point on its window against the limit; the
  ! energy balance over the window and its defect
  do j = 1, scenarios
     call mechanisms(j) % audit(solutions(j), stresses(j), observed(j), limit)
  end do

contains

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

end program flexible_mechanism
