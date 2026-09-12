# Accuracy contract

`run.sh` is the required accuracy suite registered in `verify.sh`. It drives
the built application executable with fixed, bounded parameters, reads the
records it prints, measures orders of accuracy, conservation, solver
residuals and transpose identities against analytic references, writes
`results/<set>/summary.json` and returns failure when any declaration is not
met. `cases.py` declares every case with its order or floor and the sentence
that justifies it; `contract.py` is the runner. Nothing in the library is
changed: the application gained the `check = state` word (the last instant's
jet and every block's final and initial imbalance norms) and the order
demonstration gained `ORDER_RECORD`/`ORDER_SUMMARY` lines and an
`acceptance = exploratory | required` setting.

    ./test/accuracy-contract/run.sh                # required + rejection sets
    ./test/accuracy-contract/run.sh --exploratory  # also the exploratory set
    python3 test/accuracy-contract/contract.py --set required --case T03-dirk4

Three sets. **Required** cases must pass. **Rejection** cases are deliberate
failures (wrong order above and below, an under-resolved reference, a
malformed record, a nonfinite value, an unconverged solve, a timeout, a
missing row, a missing column, a process failure, an unmet residual, a slope
at the roundoff floor, a stale limitation); each must be reported with the
status it declares, and the runner fails otherwise. **Exploratory** cases
report measured orders for configurations that are not acceptance criteria;
their statuses never fail `verify.sh`.

A required case that documents a measured shortfall is declared a
*limitation* with its measured value. It runs every time; while the
shortfall persists its status is `limitation` (not a pass); once the
declaration is met the case fails with `unexpected_pass` so the declaration
is updated rather than silently outgrown.

## Analytic references

Van der Pol at design nu = 0 is the harmonic oscillator q'' + q = 0. With
q(0) = 1, q'(0) = 0: q = cos t. Over T = 2:

- energy functional E = int_0^T (q^2 + q'^2)/2 dt = T/2 = 1;
- the tangent state w = dq/dnu at nu = 0 solves w'' + w = (1 - q^2) q',
  w(0) = w'(0) = 0: w = 3 t cos t / 8 + sin^3 t / 8 - 3 sin t / 8;
- dE/dnu = int_0^T (q w + q' w') dt = 3T^2/16 - sin^4 T/16 - 3 sin^2 T/16
  = 0.55224376640617830;
- dissipation D = int nu (1 - q^2) q'^2 dt, dD/dnu at nu = 0 =
  int_0^T sin^4 t dt = 3T/8 - sin 2T/4 + sin 4T/32 = 0.97011806903396274;
- q(T) = cos 2, q'(T) = -sin 2, and (q^2 + q'^2)/2 = 1/2 at every t.

The derivation was checked with sympy (`dsolve` of the forced oscillator
and the two integrals). The algebraic form `vanderpol_algebraic` (fields q
and y = q^2) has the same solution and is the multicomponent temporal case.

Spatial cases use the application's own comparisons: the separated mode
cos(2 pi x) cos(2 pi y) cos(omega t) on the periodic unit box with kappa =
0.1 (the printed semi-discrete error isolates the temporal error on a fixed
16 x 16 mesh, where the mode is an eigenvector of the translation-invariant
discrete operator; the error against the exact mode with the fourth-order
time scheme at ten steps isolates the spatial error, the temporal part being
below 2 % of it at 32 x 32), the discrete Laplacian of the mode against
-kappa |k|^2 times the mode, and the Taylor-Green vortex on the periodic
2 pi box (multicomponent spatial case, velocity error, relative rms).

## Orders and windows

For a discretisation of order p the error against an analytic reference has
the expansion e(h) = C h^p (1 + a h + O(h^2)). Over the finest pair of grids
with step ratio r the observed order is

    p_obs = ln(e(h) / e(h/r)) / ln r = p + ln((1 + a h) / (1 + a h / r)) / ln r.

The asymptotic regime is declared as |a h| <= theta = 1/4 on the coarser
grid of the pair: the next term is at most a quarter of the leading term.
Then |p_obs - p| <= delta(r) = -ln((1 - theta) / (1 - theta / r)) / ln r,
which is 0.2224 at r = 2. Every required refinement uses r = 2 (instants
21, 41, 81, 161, 321 in time; 8, 16, 32 cells in space; 6, 11, 21, 41
instants on the field). The full pairwise sequence, the least-squares slope
over all grids, the Richardson-extrapolated order p_K + (p_K - p_{K-1})/(r-1)
and whether the deviation contracts from the second-finest to the finest
pair are recorded as diagnostics.

Roundoff enters as an absolute floor phi on each value: half a unit of the
last printed significant digit (12 digits in the table, 17 in the state
line, 4 in the check lines) plus gamma_N |reference| with gamma_N = N u /
(1 - N u), u = 2^-53 and N the number of summed terms (instants times
components). A floor perturbs the pair's slope by at most (phi_k / e_k +
phi_{k+1} / e_{k+1}) / ln r; this widening is added to the window, and when
it exceeds delta the pair is `unresolved` and cannot pass. An error of zero
at the printed precision, or an error that does not decrease under
refinement (a sign change of the leading coefficient), is likewise not a
pass (`unresolved`, `not_monotone`). A slope above p + delta + widening is
`exceeds`: the declared order is not the leading order, and the declaration
must state the measured one.

A numerical reference refined rho times beyond the finest grid replaces the
error by |C| (h^p - h_ref^p) and shifts the finest slope upward by

    b(p, r, rho) = ln((r^p - rho^-p) / (1 - rho^-p)) / ln r - p.

It is admitted only when b <= delta, and the window is then widened by b on
the upper side; otherwise the case is refused as `under_resolved_reference`
before any run (rho = 1.5 at p = 2 gives b = 0.68; rho = 4 gives 0.07). Every
required case uses an analytic reference; the exploratory case X01 shows an
admitted numerical reference and the rejection case R03 a refused one.

## Floors

- Conservation: the trapezoidal rule (adams2) and the implicit midpoint
  (dirk2) conserve the quadratic invariant of the linear oscillator; adams2
  also integrates the energy exactly, E = T/2. Contract: |value -
  reference| <= phi.
- Transpose identity: the forward and reverse passes evaluate one bilinear
  form; `check = passes` prints their relative difference, bounded by
  gamma_N.
- Solver residual: every block's final imbalance norm <= 1e-12 x its initial
  residual norm, the configured relative stopping rule, read from the record
  rather than from the converged flag.
- Physics residual: q'' + q at the last instant is a row of the solved
  residual at nu = 0, bounded by tolerance x initial norm + gamma_N |q|.
  Declared for multistep and Newmark rows; the DIRK rows violate it at
  first order (limitation L05).
- Functional quadrature: a functional over one step is integrated on the
  instants the family's rows read (Adams and BDF: `min(order, k)` instants;
  Newmark: the two instants k - 1 and k, the trapezoidal rule), so the
  quadrature is of the order of the scheme; a staged family integrates over
  its stages with the tableau weights.
- Sensitivity demonstration: |tangent - adjoint| <= gamma_N |tangent|; the
  central difference quotient of two functionals each solved to relative
  tolerance 1e-12 satisfies |tangent - quotient| <= (tolerance + u) |F| /
  delta with delta = 1e-6, the delta^2 truncation term being far below it.

## Declared limitations at 3fb9c97

- L01: Adams-Moulton 4 design derivatives and last-instant state converge
  at order 3 while its functional reaches 4.
- L02 (resolved by R07 slice 2): every Newmark scheme's design derivatives
  converged at first order because `family_step_quadrature` integrated the
  functionals of an unstaged one-step family by the right-endpoint rectangle
  rule; the two-instant (trapezoidal) rule of the instants the Newmark rows
  read restores order 2. The case is now required as T11 (beta = 1/4,
  gamma = 1/2: dE, dD, q, q' at order 2, E and the invariant conserved) and
  T12 (Fox-Goodwin: E, dE, dD, q' at order 2), the former exploratory X03.
- L03: Adams-Moulton 3 and 4 on the periodic field converge at first order
  in time while the same schemes reach their order on the ODE.
- L04: the form-degree-4 spatial operator converges at order 2.2 on the
  periodic mode, not 4.
- L05: the DIRK jet at an arriving instant does not satisfy the law:
  |q'' + q| = 0.46 h at the last instant.

## Summary record

`summary.json` carries the set, the provenance (commit, branch, modified
working tree, compiler and version, flags, precision, executable and
source SHA-256, Python, platform, CPU affinity), the constants (theta, u,
delta at r = 2, the tolerance, the analytic references) and one entry per
case with its status, message, elapsed time, every run's argument vector,
exit status, elapsed time and log path, and every check's kind, order or
threshold, values, errors, pairwise slopes, window, roundoff widening,
reference bias, extrapolated order and justification.
