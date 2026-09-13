"""Declared accuracy cases: references, orders, floors and their justification.

Every required case drives the application executable with fixed, bounded
parameters. The temporal cases use van der Pol at design nu = 0, which is the
harmonic oscillator q'' + q = 0 with q(0) = 1, q'(0) = 0, so q = cos t and
every reference below is analytic (derivations in README.md and in the
sympy script recorded there). The spatial cases use the separated mode on
the periodic unit box and the Taylor-Green vortex, whose exact solutions the
application compares against in its own check lines.

A quantity is read from one row of the record:
  E, dE, dD      energy functional, its design derivative, the dissipation's
                 design derivative (table columns)
  q, qd, qdd     the jet at the last instant (check = state)
  invariant      (q^2 + qd^2)/2 at the last instant, exactly 1/2 for q = cos t
  law            qdd + q, the physics residual at the last instant at nu = 0
  transpose      forward against reverse first derivatives (check = passes)
  blocks         every block's final and initial imbalance norms
  semi, mode     error against the semi-discrete and the exact mode
  operator       relative rms error of the discrete Laplacian of the mode
  velocity       Taylor-Green velocity error, relative rms
  estimate:<f>, effectivity:<f>, scale:<f>, transfer:<f>, localization:<f>
                 the functional discretization-error estimate of functional
                 <f> (check = functional_error): eta, eta / (F - F_h) with F
                 the declared reference, the scale S of the relative
                 criterion, the transfer defect, and log2 of the sum of the
                 step indicators |eta_k| over an interval against the same
                 sum on the paired uniform grid (check = indicators)
  derivative_estimate:<f>, derivative_effectivity:<f>
                 the same estimate for the design derivative G_h = dF_h/dnu:
                 eta_G and eta_G / (G - G_h) with G the declared derivative
                 reference (max_derivative_degree >= 1)
"""

import math

APPLICATION = "application/graph_time_integrator"
DURATION = 2.0
INSTANTS = (21, 41, 81, 161, 321)

# analytic references at T = 2 (README.md, section "Analytic references")
E_REF = DURATION / 2
DE_REF = 3 * DURATION**2 / 16 - math.sin(DURATION)**4 / 16 - 3 * math.sin(DURATION)**2 / 16
DD_REF = 3 * DURATION / 8 - math.sin(2 * DURATION) / 4 + math.sin(4 * DURATION) / 32
Q_REF = math.cos(DURATION)
QD_REF = -math.sin(DURATION)
INVARIANT_REF = 0.5
# the exact value of each functional word at nu = 0: over T = 2, E = T/2 and
# D = 0 identically; the mean int q dt over one period T = 2 pi is 0 by
# cancellation, where int |q| dt = 4 keeps the criterion's scale finite
PERIOD = 2 * math.pi
MEAN_REF = 0.0
MEAN_SCALE = 4.0
FUNCTIONAL_REFERENCES = {"energy": E_REF, "dissipation": 0.0, "mean": MEAN_REF}
# the design derivatives at nu = 0: dE/dnu and dD/dnu of README; the mean's
# tangent w'' + w = -sin^3 t, w(0) = w'(0) = 0, is w = 3/8 t cos t - 1/32 sin 3t
# - 9/32 sin t, whose integral over [0, 2 pi] is 0 term by term
FUNCTIONAL_DERIVATIVE_REFERENCES = {"energy": DE_REF, "dissipation": DD_REF, "mean": 0.0}

# THE RADIAL OSCILLATOR q'' + q - nu / q^3 = 0 at nu = L^2, with q(0) = 1,
# q'(0) = 0: q = sqrt(cos^2 t + nu sin^2 t) (with y = q^2 the equation reads
# y y'' - y'^2/2 + 2 y^2 = 2 nu, which a cos^2 t + b sin^2 t satisfies when
# a b = nu). Every reference below follows from that solution.
RADIAL_DESIGN = 2.0


def radial_state(t, nu):
    """q, q' and q'' of the radial oscillator at time t."""
    q = math.sqrt(math.cos(t) ** 2 + nu * math.sin(t) ** 2)
    return q, (nu - 1.0) * math.sin(2 * t) / (2 * q), nu / q ** 3 - q


RADIAL_Q_REF, RADIAL_QD_REF, _ = radial_state(DURATION, RADIAL_DESIGN)
# E = q'^2/2 + q^2/2 + nu/(2 q^2) = (1 + nu)/2 at every instant, so the energy
# functional is T E and its design derivative is T/2
RADIAL_INVARIANT_REF = 0.5 * (1.0 + RADIAL_DESIGN)
RADIAL_E_REF = DURATION * RADIAL_INVARIANT_REF
RADIAL_DE_REF = DURATION / 2
# F_2 = int_0^T q^2 dt = int (cos^2 t + nu sin^2 t) dt
RADIAL_F2_REF = (1.0 + RADIAL_DESIGN) * DURATION / 2 \
    + (1.0 - RADIAL_DESIGN) * math.sin(2 * DURATION) / 4
RADIAL_DF2_REF = DURATION / 2 - math.sin(2 * DURATION) / 4

# the configured Newton stopping rule of every run: relative 1e-12
TOLERANCE = 1.0e-12
STATE_COMPONENTS = 3

THETA_TEXT = ("window from the two-term expansion e = C h^p (1 + a h) with |a h| <= 1/4 "
              "on the coarser grid of the finest pair; roundoff floor = print resolution "
              "+ gamma_N |reference|")


def ode_argv(instants, families="bdf adams dirk", max_order=4, derivative=1,
             physics="vanderpol", chain=None, extra=(), check="state passes",
             duration=DURATION, functionals="energy dissipation"):
    argv = [APPLICATION, f"--physics={physics}", "--design=0.0", "--initial_state=1.0",
            f"--time_duration={duration}", f"--instants={instants}", "--grid=uniform",
            f"--families={families}", f"--max_discretization_order={max_order}",
            f"--max_derivative_degree={derivative}", f"--functionals={functionals}",
            "--combinations=1", f"--check={check}"]
    if chain:
        argv.append(f"--chain={chain}")
    return argv + list(extra)


def radial_argv(instants, families="bdf adams dirk", max_order=4, derivative=1, extra=(),
                check="state passes"):
    return [APPLICATION, "--physics=radial_oscillator", f"--design={RADIAL_DESIGN}",
            "--initial_state=1.0", f"--time_duration={DURATION}", f"--instants={instants}",
            "--grid=uniform", f"--families={families}", f"--max_discretization_order={max_order}",
            f"--max_derivative_degree={derivative}", "--functionals=energy square_integral",
            "--combinations=1", f"--check={check}"] + list(extra)


def ode_grids(instants=INSTANTS, duration=DURATION):
    return [(f"instants={n}", duration / (n - 1), n * STATE_COMPONENTS) for n in instants]


def field_argv(cells, instants, families="dirk bdf", max_order=4, extra=(), check="mode state"):
    return [APPLICATION, "--config=torus2", f"--spatial_counts={cells} {cells}",
            "--spatial_extent=1.0 1.0", f"--instants={instants}", f"--families={families}",
            f"--max_discretization_order={max_order}", "--max_derivative_degree=0",
            f"--check={check}"] + list(extra)


SPARSE = ("--linear_solver=iterative", "--storage=sparse", "--preconditioner=gauss_seidel")


def taylor_green_argv(cells, instants=6):
    return [APPLICATION, "--config=taylor_green", f"--spatial_counts={cells} {cells}",
            f"--instants={instants}", "--export=none", "--max_derivative_degree=0",
            "--families=bdf", "--max_discretization_order=2", "--check=exact state"]


def order_check(quantity, order, reference, digits, justification, design=None):
    return {"kind": "order", "quantity": quantity, "order": order,
            "reference": reference, "digits": digits, "justification": justification,
            "design": design}


def floor_check(quantity, reference, digits, justification, scale=None, design=None):
    return {"kind": "floor", "quantity": quantity, "reference": reference,
            "digits": digits, "scale": scale, "justification": justification,
            "design": design}


def residual_check(tolerance=TOLERANCE):
    return {"kind": "residual", "quantity": "blocks", "tolerance": tolerance,
            "justification": "every block's final imbalance norm <= tolerance x its initial "
                             "residual norm, the configured relative stopping rule, read from "
                             "the record rather than from the converged flag"}


TRANSPOSE = floor_check("transpose", 0.0, 3, "forward and reverse passes of one bilinear form "
                        "differ by at most gamma_N times its magnitude, N = instants x components",
                        scale=1.0)
LAW = floor_check("law", 0.0, 17, "qdd + q at the last instant is a row of the solved residual "
                  "at nu = 0: bounded by tolerance x initial residual norm of its block "
                  "plus gamma_N |q|", scale="law")

TABLE_DIGITS = 12
STATE_DIGITS = 17
CHECK_DIGITS = 4


def temporal_case(identifier, row, order, quantities, description, chain=None,
                  limitation=None, argv_extra=(), instants=INSTANTS, set_name="required"):
    grids = ode_grids(instants)
    if chain:
        runs = {label: ode_argv(n, families="bdf", max_order=1, chain=chain, extra=argv_extra)
                for (label, _, _), n in zip(grids, instants)}
    else:
        runs = {label: ode_argv(n, extra=argv_extra) for (label, _, _), n in zip(grids, instants)}
    checks = []
    for quantity in quantities:
        if quantity == "E":
            checks.append(order_check("E", order, E_REF, TABLE_DIGITS,
                                      f"energy functional against T/2; {THETA_TEXT}"))
        elif quantity == "dE":
            checks.append(order_check("dE", order, DE_REF, TABLE_DIGITS,
                                      f"dE/dnu against the forced-oscillator value; {THETA_TEXT}"))
        elif quantity == "dD":
            checks.append(order_check("dD", order, DD_REF, TABLE_DIGITS,
                                      f"dD/dnu against int sin^4; {THETA_TEXT}"))
        elif quantity == "q":
            checks.append(order_check("q", order, Q_REF, STATE_DIGITS,
                                      f"q(T) against cos T; {THETA_TEXT}"))
        elif quantity == "qd":
            checks.append(order_check("qd", order, QD_REF, STATE_DIGITS,
                                      f"q'(T) against -sin T; {THETA_TEXT}"))
        elif quantity == "invariant":
            checks.append(order_check("invariant", order, INVARIANT_REF, STATE_DIGITS,
                                      f"energy drift (q^2 + q'^2)/2 - 1/2 at T; {THETA_TEXT}"))
        elif quantity == "conserved":
            checks.append(floor_check("E", E_REF, TABLE_DIGITS,
                                      "the trapezoidal discrete energy of the linear oscillator "
                                      "is conserved exactly: |F - T/2| <= print resolution + "
                                      "gamma_N |F|", scale=E_REF))
        elif quantity == "invariant_conserved":
            checks.append(floor_check("invariant", INVARIANT_REF, STATE_DIGITS,
                                      "the scheme conserves the quadratic invariant of the "
                                      "linear oscillator: |(q^2 + q'^2)/2 - 1/2| <= print "
                                      "resolution + gamma_N / 2", scale=INVARIANT_REF))
        else:
            raise ValueError(quantity)
    checks += [LAW, TRANSPOSE, residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": limitation,
            "timeout": 120}


RADIAL_LAW = floor_check("radial_law", 0.0, 17, "q'' + q - nu / q^3 at the last instant is a "
                         "row of the solved residual: bounded by tolerance x initial residual "
                         "norm of its block plus gamma_N |q|", scale="law", design=RADIAL_DESIGN)


def radial_case(identifier, row, order, quantities, description, limitation=None,
                argv_extra=(), instants=INSTANTS, set_name="required"):
    """The radial oscillator at nu = 2 over the refined grids: the two
    functionals, their design derivatives, the state and the law."""
    grids = ode_grids(instants)
    runs = {label: radial_argv(n, extra=argv_extra) for (label, _, _), n in zip(grids, instants)}
    references = {"E": (RADIAL_E_REF, TABLE_DIGITS, "the energy functional against T (1 + nu)/2"),
                  "dE": (RADIAL_DE_REF, TABLE_DIGITS, "dF_E/dnu against T/2"),
                  "F2": (RADIAL_F2_REF, TABLE_DIGITS,
                         "int q^2 against (1 + nu) T/2 + (1 - nu) sin 2T/4"),
                  "dF2": (RADIAL_DF2_REF, TABLE_DIGITS,
                          "dF_2/dnu against T/2 - sin 2T/4"),
                  "q": (RADIAL_Q_REF, STATE_DIGITS,
                        "q(T) against sqrt(cos^2 T + nu sin^2 T)"),
                  "qd": (RADIAL_QD_REF, STATE_DIGITS,
                         "q'(T) against (nu - 1) sin 2T / (2 q(T))"),
                  "radial_invariant": (RADIAL_INVARIANT_REF, STATE_DIGITS,
                                       "the invariant drift q'^2/2 + q^2/2 + nu/(2 q^2) - "
                                       "(1 + nu)/2 at T")}
    checks = []
    for quantity in quantities:
        reference, digits, text = references[quantity]
        checks.append(order_check(quantity, order, reference, digits,
                                  f"{text}; {THETA_TEXT}", design=RADIAL_DESIGN))
    checks += [RADIAL_LAW, TRANSPOSE, residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": limitation,
            "timeout": 180}


def radial_estimator_case(identifier, row, order, description, families="bdf", max_order=4,
                          instants=INSTANTS, set_name="required"):
    """The functional-error estimate on the radial oscillator, the one law whose
    exact functionals are known at every grid: the effectivity of each against
    its own closed form at order 1, the estimate at the order of F_h - F, the
    transfer identity of each enriched block, and the law at the last instant.
    The word energy names a different number here than on van der Pol, so each
    effectivity carries its exact value as functional_reference."""
    grids = ode_grids(instants)
    runs = {label: radial_argv(n, families=families, max_order=max_order,
                               check="state functional_error")
            for (label, _, _), n in zip(grids, instants)}
    checks = []
    for word, reference in (("energy", RADIAL_E_REF), ("square_integral", RADIAL_F2_REF)):
        checks.append(dict(order_check(f"effectivity:{word}", 1, 1.0, TABLE_DIGITS,
                                       EFFECTIVITY_TEXT),
                           functional_reference=reference))
        checks.append(order_check(f"estimate:{word}", order, 0.0, TABLE_DIGITS,
                                  "the estimate eta converges to zero at the order of F_h - F; "
                                  + THETA_TEXT))
        checks.append(transfer_check(word))
    checks += [RADIAL_LAW, residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": None,
            "timeout": 180}


def field_temporal_case(identifier, row, order, description, limitation=None,
                        set_name="required", families="dirk bdf", max_order=4):
    instants = (6, 11, 21, 41)
    grids = [(f"instants={n}", 0.5 / (n - 1), 256 * n * STATE_COMPONENTS) for n in instants]
    runs = {f"instants={n}": field_argv(16, n, families, max_order, SPARSE, check="mode state")
            for n in instants}
    checks = [order_check("semi", order, None, CHECK_DIGITS,
                          "error against the semi-discrete mode isolates the temporal error "
                          f"on the fixed 16 x 16 periodic mesh; {THETA_TEXT}"),
              residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": limitation,
            "timeout": 120}


def spatial_case(identifier, row, quantity, order, description, cells=(8, 16, 32),
                 instants=11, families="dirk", max_order=4, extra=SPARSE, spatial_order=None,
                 limitation=None, set_name="required", timeout=120):
    grids = [(f"cells={n}", 1.0 / n, n * n * instants * STATE_COMPONENTS) for n in cells]
    more = list(extra) + ([f"--spatial_order={spatial_order}"] if spatial_order else [])
    runs = {f"cells={n}": field_argv(n, instants, families, max_order, more,
                                     check=f"{quantity} state") for n in cells}
    if quantity == "operator":
        text = ("relative rms error of the discrete Laplacian applied to the mode against "
                f"-kappa |k|^2 mode, interior cells of the periodic box; {THETA_TEXT}")
    else:
        text = ("error against the exact mode at T = 0.5 with the fourth-order time scheme at "
                "ten steps, whose semi-discrete error is below 2 % of the spatial error at "
                f"32 x 32; {THETA_TEXT}")
    checks = [order_check(quantity, order, None, CHECK_DIGITS, text), residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": limitation,
            "timeout": timeout}


# THE DISC: the polar mesh of n1 radial and n2 angular cells, so the radial
# step is a/n1 and the angular arc 2 pi a/n2, both halving along the sequence.
# The cells are 1 + (n1 - 1) n2, the centre being one polygonal cell.
DISC_COUNTS = ((9, 16), (17, 32), (33, 64))


def disc_argv(counts, instants=3, families="bdf", max_order=1, check="operator",
              derivative=0, spatial_order=2, initial_field="constant", rows=None,
              extra=SPARSE):
    argv = [APPLICATION, "--config=disc", f"--spatial_counts={counts[0]} {counts[1]}",
            f"--instants={instants}", f"--families={families}",
            f"--max_discretization_order={max_order}",
            f"--max_derivative_degree={derivative}",
            f"--spatial_order={spatial_order}", f"--initial_field={initial_field}",
            f"--check={check}"] + list(extra)
    if rows is not None:
        argv.append(f"--rows={rows}")
    return argv


def disc_cells(counts):
    return 1 + (counts[0] - 1) * counts[1]


def disc_operator_case(identifier, quantity, description, order=2, limitation=None,
                       counts=DISC_COUNTS, set_name="required"):
    """One class of cells of the disc against kappa times the laplacian of the
    quartic (r^2 - a^2)^2, whose radial derivative vanishes at r = a."""
    grids = [(f"cells={n1}x{n2}", 1.0 / n1, disc_cells((n1, n2))) for n1, n2 in counts]
    runs = {f"cells={n1}x{n2}": disc_argv((n1, n2)) for n1, n2 in counts}
    text = ("relative rms error of the fitted balance applied to (r^2 - a^2)^2 against "
            "kappa (16 r^2 - 8 a^2), by cell class, at form degree 2 on the polar mesh; "
            f"{THETA_TEXT}")
    return {"id": identifier, "set": set_name, "description": description, "row": "bdf1",
            "grids": grids, "runs": runs,
            "checks": [order_check(quantity, order, None, CHECK_DIGITS, text)],
            "limitation": limitation, "timeout": 600}


def disc_conservation_case(counts=DISC_COUNTS):
    """The discrete divergence theorem on the disc: an interior face is counted
    twice with opposite signs and a boundary face carries the zero Neumann flux,
    so the balance summed over the cells is zero for every field."""
    grids = [(f"cells={n1}x{n2}", 1.0 / n1, 2 * disc_cells((n1, n2))) for n1, n2 in counts]
    runs = {f"cells={n1}x{n2}": disc_argv((n1, n2)) for n1, n2 in counts}
    text = ("the balance summed over every cell, relative to the sum of the magnitudes: "
            "zero for every field, bounded by gamma_N with N twice the cells")
    return {"id": "D04-disc-conservation", "set": "required",
            "description": "the discrete divergence theorem on the disc, for the quartic and "
                           "for a field of the run's own seed",
            "row": "bdf1", "grids": grids, "runs": runs,
            "checks": [floor_check("balance_sum", 0.0, CHECK_DIGITS, text, scale=1.0),
                       floor_check("balance_sum_seeded", 0.0, CHECK_DIGITS, text, scale=1.0)],
            "limitation": None, "timeout": 600}


def disc_mode_case(identifier, description, order=2, rows=None, counts=DISC_COUNTS,
                   limitation=None, set_name="required"):
    """The radially symmetric Neumann mode of the disc, J_0(z_1 r / a) cos(omega t)
    with omega^2 = 1 + kappa (z_1/a)^2, marched by a fourth-order family at ten
    steps so the spatial error is the one measured."""
    instants = 11
    grids = [(f"cells={n1}x{n2}", 1.0 / n1, disc_cells((n1, n2)) * instants * STATE_COMPONENTS)
             for n1, n2 in counts]
    runs = {f"cells={n1}x{n2}": disc_argv((n1, n2), instants=instants, families="dirk",
                                          max_order=4, check="mode state",
                                          initial_field="mode", rows=rows)
            for n1, n2 in counts}
    text = ("error against J_0(z_1 r / a) cos(omega t) at T = 0.5, DIRK-4 at ten steps, "
            "whose temporal error is below 8 % of the spatial error already at five steps; "
            f"{THETA_TEXT}")
    return {"id": identifier, "set": set_name, "description": description, "row": "dirk4",
            "grids": grids, "runs": runs,
            "checks": [order_check("mode", order, None, CHECK_DIGITS, text), residual_check()],
            "limitation": limitation, "timeout": 900}


def disc_transpose_case(identifier, description, rows=None, counts=((9, 16), (17, 32))):
    """One bilinear form evaluated by the forward and the reverse pass on the disc."""
    instants = 6
    grids = [(f"cells={n1}x{n2}", 1.0 / n1, disc_cells((n1, n2)) * instants * STATE_COMPONENTS)
             for n1, n2 in counts]
    runs = {f"cells={n1}x{n2}": disc_argv((n1, n2), instants=instants, families="dirk",
                                          max_order=2, derivative=1,
                                          check="mode state passes", initial_field="mode",
                                          rows=rows)
            for n1, n2 in counts}
    return {"id": identifier, "set": "required", "description": description, "row": "dirk2",
            "grids": grids, "runs": runs, "checks": [TRANSPOSE, residual_check()],
            "limitation": None, "timeout": 900}


def taylor_green_case(identifier, cells, description, set_name="required", timeout=120,
                      instants=6):
    grids = [(f"cells={n}", 2 * math.pi / n, n * n * instants * 6) for n in cells]
    runs = {f"cells={n}": taylor_green_argv(n, instants) for n in cells}
    checks = [order_check("velocity", 2, None, CHECK_DIGITS,
                          "velocity error, relative rms, against the exact vortex at T = 1, "
                          f"nu = 0.01, bdf2 with {instants - 1} steps; {THETA_TEXT}"),
              residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": "bdf2",
            "grids": grids, "runs": runs, "checks": checks, "limitation": None,
            "timeout": timeout}


EFFECTIVITY_TEXT = ("effectivity eta / (F - F_h) of the functional-error estimate against 1: "
                    "enrichment A (the family of order p + 1 on the same grid, the identity "
                    "prolongation) is asymptotically exact, I = 1 + O(h^(min(p+, 2p) - p)); "
                    "the coarse block's first own instant stays fixed in the enriched block "
                    "at its local order p + 1, so |I - 1| converges at order 1; " + THETA_TEXT)


TRANSFER_TEXT = ("the fixed-row residual of every enriched block, Q+_i - Q_j at P Q_h, "
                 "vanishes under the identity prolongation: zero within gamma_N |Q|, "
                 "|Q| <= 1 on the oscillator")

TRANSFER = floor_check("transfer:energy", 0.0, TABLE_DIGITS, TRANSFER_TEXT, scale=1.0)


def transfer_check(word):
    """The transfer identity on the enriched chain of one functional's estimate."""
    return floor_check(f"transfer:{word}", 0.0, TABLE_DIGITS, TRANSFER_TEXT, scale=1.0)


DERIVATIVE_EFFECTIVITY_TEXT = (
    "effectivity eta_G / (G - G_h) of the derivative-functional estimate against 1, "
    "G = dF/dnu: eta_G is the derivative of the estimate, d/dnu[L+(P Q_h, lambda+) - F_h] "
    "= -lambda+'^T R+(P Q_h) + [F+_nu - lambda+^T R+_nu](P Q_h) - G_h, so it is the "
    "derivative of an asymptotically exact estimate and |I - 1| = O(h) where that O(h) "
    "term is smooth in nu; " + THETA_TEXT)


def estimator_case(identifier, row, order, description, families="bdf adams", max_order=4,
                   instants=(21, 41, 81, 161), argv_extra=(), set_name="required",
                   effectivity=True, limitation=None, chain=None):
    """The functional-error estimate of one row: its effectivity at order 1, the
    estimate itself at the order of F_h - F and the transfer identity (README,
    section 'Functional discretization error')."""
    grids = ode_grids(instants)
    if chain:
        runs = {label: ode_argv(n, families="bdf", max_order=1, chain=chain,
                                check="state functional_error", extra=argv_extra)
                for (label, _, _), n in zip(grids, instants)}
    else:
        runs = {label: ode_argv(n, families=families, max_order=max_order,
                                check="state functional_error", extra=argv_extra)
                for (label, _, _), n in zip(grids, instants)}
    checks = []
    if effectivity:
        checks.append(order_check("effectivity:energy", 1, 1.0, TABLE_DIGITS, EFFECTIVITY_TEXT))
    checks.append(order_check("estimate:energy", order, 0.0, TABLE_DIGITS,
                              "the estimate eta converges to zero at the order of F_h - F; "
                              + THETA_TEXT))
    checks += [TRANSFER, LAW, residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": limitation,
            "timeout": 120}


def field_estimator_case(identifier, row, order, description, families="dirk bdf", max_order=3,
                         set_name="required", effectivity=True):
    """Enrichment A in time on the fixed 16 x 16 periodic mesh: the estimate of the
    energy functional against the semi-discrete mode energy printed by the
    application (omega_h from the discrete operator, of which the mode is an
    eigenvector), which isolates the temporal error; effectivity at order 1."""
    instants = (6, 11, 21, 41)
    grids = [(f"instants={n}", 0.5 / (n - 1), 256 * n * STATE_COMPONENTS) for n in instants]
    runs = {f"instants={n}": field_argv(16, n, families, max_order, SPARSE,
                                        check="mode state functional_error") for n in instants}
    checks = []
    if effectivity:
        checks.append(dict(order_check("effectivity:energy", 1, 1.0, TABLE_DIGITS,
                                       "effectivity against the semi-discrete mode energy on the fixed "
                                       "16 x 16 mesh, the temporal part of the error; " + EFFECTIVITY_TEXT),
                           functional_reference="semi_energy"))
    checks.append(order_check("estimate:energy", order, 0.0, TABLE_DIGITS,
                              "the estimate eta converges to zero at the temporal order of E_h - "
                              "E_semi; " + THETA_TEXT))
    checks += [TRANSFER, residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": None, "timeout": 300}


def derivative_estimator_case(identifier, row, word, order, description, families="bdf",
                              max_order=3, instants=(21, 41, 81, 161, 321),
                              set_name="required", effectivity=True):
    """The estimate of the discretization error of a design derivative G_h =
    dF_h/dnu (the order-1 Lagrangian), through the enriched costate rate
    lambda+' of the derivative Lagrangian: its effectivity at order 1 and the
    estimate itself at the order of G_h - G."""
    grids = ode_grids(instants)
    runs = {label: ode_argv(n, families=families, max_order=max_order,
                            check="state functional_error")
            for (label, _, _), n in zip(grids, instants)}
    checks = []
    if effectivity:
        checks.append(order_check(f"derivative_effectivity:{word}", 1, 1.0, TABLE_DIGITS,
                                  DERIVATIVE_EFFECTIVITY_TEXT))
    checks.append(order_check(f"derivative_estimate:{word}", order, 0.0, TABLE_DIGITS,
                              "the estimate eta_G converges to zero at the order of G_h - G; "
                              + THETA_TEXT))
    checks += [transfer_check(word), LAW, residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": None,
            "timeout": 300}


def mean_case(identifier, row, order, description, families="bdf", max_order=2,
              instants=(21, 41, 81, 161), set_name="required"):
    """The near-zero functional: mean = int q dt over one period, whose exact value
    is 0 by cancellation. The criterion divides by the scale S = sum |w_kj f(Q_j)|,
    never by |F_h|; S converges to int |q| dt = 4 at the family's order where the
    quadrature weights are non-negative (the two-instant rule of BDF-2)."""
    grids = ode_grids(instants, duration=PERIOD)
    runs = {label: ode_argv(n, families=families, max_order=max_order, duration=PERIOD,
                            functionals="mean", check="state functional_error")
            for (label, _, _), n in zip(grids, instants)}
    checks = [order_check("estimate:mean", order, MEAN_REF, TABLE_DIGITS,
                          "the estimate eta of a functional whose exact value is zero "
                          "converges to zero at the order of F_h - F, its floor gamma_N S "
                          "rather than gamma_N |F_h|; " + THETA_TEXT),
              order_check("scale:mean", order, MEAN_SCALE, TABLE_DIGITS,
                          "the scale of the relative criterion converges to int |q| dt = 4 "
                          "at the family's order: the two-instant rule of BDF-2 has "
                          "non-negative weights, so sum |w_kj f(Q_j)| tends to int |f|; "
                          + THETA_TEXT),
              transfer_check("mean"), residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": None,
            "timeout": 120}


def localized_case(identifier, row, order, interval, description, instants=(41, 81, 161),
                   families="bdf", max_order=3, set_name="required"):
    """The step indicators on the coarsened grid (steps merged in pairs inside the
    interval) against those of the same interval on the uniform grid of the same
    instant count: log2 of the ratio of the sums of |eta_k| tends to p, the
    indicators being C(t) h_k^(p+1) with the transition steps at the bounds
    contributing O(h) of the sum."""
    grids = ode_grids(instants)
    runs = {}
    for (label, _, _), n in zip(grids, instants):
        runs[label] = ode_argv(n, families=families, max_order=max_order,
                               check="state functional_error indicators",
                               extra=("--grid=coarsened",
                                      f"--coarsened_interval={interval[0]} {interval[1]}"))
        runs[f"uniform={n}"] = ode_argv(n, families=families, max_order=max_order,
                                        check="state functional_error indicators")
    checks = [{"kind": "order", "quantity": "localization:energy", "order": 1, "reference": order,
               "digits": TABLE_DIGITS, "paired": "uniform", "interval": list(interval),
               "justification": "log2(sum of |eta_k| over the coarsened steps of the interval / "
                                "the same sum on the uniform grid) against p: the indicator of "
                                "a step is C(t) h_k^(p+1), so pairs merged into 2h give 2^p, "
                                f"the transition steps at the bounds O(h); {THETA_TEXT}"},
              residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": None, "timeout": 120}


def adaptation_case(identifier, row, tolerance, description, families, max_order, instants=21,
                    limit=None, set_name="required"):
    """The adaptive grid driven by the estimator (adaptive_check = functional_error):
    on the accepted grid |E - E_h| <= tolerance x S, S the scale of the relative
    criterion; the estimator/error ratio is recorded. A limit below the need is
    the rejection case: no acceptance, the process fails with ADAPTATION_UNMET."""
    extra = ["--grid=adaptive", "--adaptive_check=functional_error",
             f"--functional_error_tolerance={tolerance}"]
    if limit:
        extra.append(f"--adaptation_instants={limit}")
    # the accepted grid is not known beforehand: gamma_N is taken at the limit
    terms = (limit if limit else 1281) * STATE_COMPONENTS
    grids = [("adaptive", DURATION / (instants - 1), terms)]
    runs = {"adaptive": ode_argv(instants, families=families, max_order=max_order,
                                 check="state functional_error", extra=extra)}
    checks = [{"kind": "floor", "quantity": "E", "reference": E_REF, "digits": TABLE_DIGITS,
               "scale": ("estimate", tolerance),
               "justification": f"|E - E_h| <= {tolerance:g} x S on the accepted grid, the "
                                "declared criterion of the adaptive loop, S = sum |w_k f(Q_k)|"},
              {"kind": "value", "quantity": "effectivity:energy", "digits": TABLE_DIGITS,
               "justification": "the estimator/error ratio on the accepted grid, recorded"},
              LAW, residual_check()]
    return {"id": identifier, "set": set_name, "description": description, "row": row,
            "grids": grids, "runs": runs, "checks": checks, "limitation": None, "timeout": 120}


def sensitivity_case():
    delta = 1.0e-6
    return {"id": "P02-sensitivity-demo", "set": "required", "record": "sensitivity",
            "description": "the sensitivity demonstration: tangent = adjoint to roundoff and "
                           "the central difference quotient within its declared floor",
            "row": None, "grids": [("demo", None, 21 * STATE_COMPONENTS)],
            "runs": {"demo": [APPLICATION, "--demo=sensitivity"]},
            "checks": [floor_check("tangent_adjoint", 0.0, 3, "|tangent - adjoint| <= gamma_N "
                                   "|tangent|, N = 21 instants x 3 components", scale="tangent"),
                       {"kind": "floor", "quantity": "tangent_difference", "reference": 0.0,
                        "digits": 3, "scale": ("difference", delta),
                        "justification": "central difference of two functionals each solved "
                                         f"to relative tolerance {TOLERANCE:g}: "
                                         "|tangent - quotient| <= (tolerance + u) |F| / delta, "
                                         f"delta = {delta:g}; the delta^2 truncation term is "
                                         "below this floor"}],
            "limitation": None, "timeout": 60}


def required_cases():
    cases = [
        temporal_case("T01-dirk2", "dirk2", 2, ["E", "dE", "dD", "q", "qd", "invariant_conserved"],
                      "implicit midpoint, one stage, order 2; conserves the quadratic invariant"),
        temporal_case("T02-dirk3", "dirk3", 3, ["E", "dE", "dD", "q", "qd", "invariant"],
                      "Crouzeix two-stage DIRK, order 3"),
        temporal_case("T03-dirk4", "dirk4", 4, ["E", "dE", "dD", "q", "qd"],
                      "Crouzeix three-stage DIRK, order 4; its invariant drifts at order 5 "
                      "and is not declared"),
        temporal_case("T04-bdf1", "bdf1", 1, ["E", "dD", "q", "qd", "invariant"],
                      "BDF-1, order 1; dE/dnu changes sign inside these grids and is not "
                      "declared"),
        temporal_case("T05-bdf2", "bdf2", 2, ["dE", "dD", "q", "qd"],
                      "BDF-2, order 2; its energy functional and invariant drift converge at "
                      "order 3 on the linear oscillator and are not declared"),
        temporal_case("T06-bdf3", "bdf3", 3, ["E", "dE", "dD", "q", "qd", "invariant"],
                      "BDF-3, order 3"),
        temporal_case("T07-bdf4", "bdf4", 4, ["dE", "dD", "q", "qd"],
                      "BDF-4, order 4; its energy functional reaches the print resolution "
                      "before the finest grid and its invariant drifts at order 5: neither "
                      "is declared"),
        temporal_case("T08-adams2", "adams2", 2, ["dE", "dD", "q", "qd", "conserved",
                                                  "invariant_conserved"],
                      "Adams-Moulton 2 (trapezoidal), order 2; the discrete energy of the "
                      "linear oscillator is conserved exactly"),
        temporal_case("T09-adams3", "adams3", 3, ["E", "dE", "dD", "q", "qd", "invariant"],
                      "Adams-Moulton 3, order 3"),
        temporal_case("T10-adams4", "adams4", 4, ["dE", "dD", "q", "qd"],
                      "Adams-Moulton 4, order 4 in the design derivatives and the state; its "
                      "energy functional reaches the print resolution before the finest grid "
                      "(X12 states its order 5 above the floor) and its invariant drifts at "
                      "order 5 (T13)"),
        temporal_case("T13-adams4-invariant-drift", "adams4", 5, ["invariant"],
                      "Adams-Moulton 4 invariant drift on the linear oscillator: order 5, one "
                      "above the scheme, as for BDF-4 and DIRK-4"),
        temporal_case("T14-alexander2", "alexander2", 2, ["E", "dE", "dD", "q", "qd"],
                      "Alexander's two-stage L-stable DIRK (gamma = 1 - sqrt 2 / 2), order 2: a "
                      "tableau registered by its data alone in the application's name map",
                      argv_extra=("--families=alexander",)),
        temporal_case("T11-newmark2", "newmark2", 2, ["dE", "dD", "q", "qd", "conserved",
                                                      "invariant_conserved"],
                      "Newmark beta = 1/4, gamma = 1/2 (average acceleration), order 2: the "
                      "trapezoidal rule on the first-order system, so the discrete energy and "
                      "the invariant of the linear oscillator are conserved exactly; functionals "
                      "integrated by the two-instant (trapezoidal) step quadrature",
                      argv_extra=("--families=newmark",)),
        temporal_case("T12-newmark3", "newmark3", 2, ["E", "dE", "dD", "qd"],
                      "Newmark beta = 1/12, gamma = 1/2 (Fox-Goodwin), order 2 in the "
                      "functionals, the design derivatives and the velocity; its position "
                      "reaches order 4 on the linear oscillator and is not declared",
                      argv_extra=("--families=newmark",)),
        temporal_case("C01-dirk3-bdf2", "dirk3-bdf2", 2, ["dE", "dD"],
                      "chain DIRK-3 then BDF-2: order min(3, 2) = 2", chain="dirk:3 bdf:2"),
        temporal_case("C02-bdf2-dirk3", "bdf2-dirk3", 2, ["dE", "dD"],
                      "chain BDF-2 then DIRK-3: order 2", chain="bdf:2 dirk:3"),
        temporal_case("C03-adams3-dirk4", "adams3-dirk4", 3, ["E", "dD"],
                      "chain Adams-3 then DIRK-4: order 3; dE/dnu's next term exceeds the "
                      "declared quarter at these grids and is not declared",
                      chain="adams:3 dirk:4"),
        temporal_case("C04-dirk4-bdf4-adams4", "dirk4-bdf4-adams4", 4, ["E", "dD"],
                      "chain DIRK-4, BDF-4, Adams-4: order 4", chain="dirk:4 bdf:4 adams:4"),
        temporal_case("C05-dirk2-adams2-bdf2", "dirk2-adams2-bdf2", 2, ["dE", "dD"],
                      "chain DIRK-2, Adams-2, BDF-2: order 2", chain="dirk:2 adams:2 bdf:2"),
        temporal_case("M01-algebraic-dirk3", "dirk3", 3, ["E", "dE", "dD"],
                      "two state fields q and y = q^2 tied algebraically, DIRK-3",
                      argv_extra=("--physics=vanderpol_algebraic", "--families=dirk",
                                  "--max_discretization_order=3"),
                      instants=(21, 41, 81, 161)),
        temporal_case("M02-algebraic-bdf2", "bdf2", 2, ["dE", "dD", "q", "qd"],
                      "two state fields q and y = q^2 tied algebraically, BDF-2",
                      argv_extra=("--physics=vanderpol_algebraic", "--families=bdf",
                                  "--max_discretization_order=2")),
        field_temporal_case("F01-field-dirk2", "dirk2", 2,
                            "temporal order on the 16 x 16 periodic mode, DIRK-2"),
        field_temporal_case("F02-field-dirk4", "dirk4", 4,
                            "temporal order on the 16 x 16 periodic mode, DIRK-4"),
        field_temporal_case("F03-field-bdf2", "bdf2", 2,
                            "temporal order on the 16 x 16 periodic mode, BDF-2"),
        field_temporal_case("F04-field-bdf4", "bdf4", 4,
                            "temporal order on the 16 x 16 periodic mode, BDF-4"),
        field_temporal_case("F05-field-adams3", "adams3", 3,
                            "temporal order on the 16 x 16 periodic mode, Adams-Moulton 3 "
                            "from the staged startup", families="adams", max_order=3),
        spatial_case("S01-operator-degree2", "bdf1", "operator", 2,
                     "discrete Laplacian of the mode, form degree 2, cells 8, 16, 32",
                     instants=3, families="bdf", max_order=1, extra=()),
        spatial_case("S02-mode-dirk4", "dirk4", "mode", 2,
                     "marched mode at T = 0.5, DIRK-4 ten steps, cells 8, 16, 32"),
        taylor_green_case("S03-taylor-green", (8, 16),
                          "Taylor-Green vortex, BDF-2, cells 8 and 16: one pair, so the "
                          "asymptotic regime is not verified by a second pair"),
        sensitivity_case(),
        estimator_case("G01-estimate-bdf1", "bdf1", 1,
                       "functional-error estimate of BDF-1 by BDF-2: effectivity at order 1 "
                       "(measured 0.83, 0.91, 0.95, 0.98), estimate at order 1"),
        estimator_case("G02-estimate-bdf3", "bdf3", 3,
                       "functional-error estimate of BDF-3 by BDF-4: effectivity at order 1 "
                       "(measured 0.86, 0.94, 0.97, 0.985), estimate at order 3"),
        estimator_case("G03-estimate-adams3", "adams3", 3,
                       "functional-error estimate of Adams-Moulton 3 by Adams-Moulton 4: "
                       "effectivity at order 1 (measured 0.89, 0.95, 0.975, 0.987), estimate "
                       "at order 3"),
        estimator_case("G04-estimate-newmark3", "newmark3", 2,
                       "functional-error estimate of Newmark beta = 1/12 (Fox-Goodwin) by "
                       "Adams-Moulton 3: effectivity at order 1 (measured 1.06, 1.04, 1.02, "
                       "1.01), estimate at order 2", families="newmark", max_order=3),
        localized_case("G05-localized-bdf3", "bdf3", 3, (0.7, 1.1),
                       "localized error: BDF-3 on the uniform grid with the steps inside "
                       "[0.7, 1.1] merged in pairs; log2 of the indicator sums over the "
                       "interval tends to 3 at order 1 (measured 2.50, 2.77 at 41, 81)"),
        estimator_case("G06-estimate-dirk2", "dirk2", 2,
                       "the grid-stationary counterexample under refinement: implicit midpoint "
                       "estimated by BDF-3 on its arriving instants, whose quadrature part is "
                       "the whole error F+(P Q_h) - F_h; effectivity at order 1 (measured "
                       "1.15, 1.09, 1.05, 1.02), estimate at order 2", families="dirk"),
        estimator_case("G07-estimate-dirk3", "dirk3", 3,
                       "Crouzeix two-stage DIRK estimated by BDF-4 on its arriving instants: "
                       "effectivity at order 1 (measured 0.75, 0.85, 0.92, 0.96), estimate at "
                       "order 3", families="dirk"),
        estimator_case("G08-estimate-dirk4", "dirk4", 4,
                       "Crouzeix three-stage DIRK estimated by BDF-5: effectivity at order 1 "
                       "(measured 1.08, 1.03, 1.01, 1.006, slopes 1.41, 1.24, 1.15), estimate "
                       "at order 4", families="dirk"),
        estimator_case("G09-estimate-adams3-dirk4", "adams3-dirk4", 3,
                       "chain Adams-3 then DIRK-4 by Adams-4 then BDF-5: the junction "
                       "sensitivity flows through the coarse-family block of the enriched "
                       "chain; effectivity at order 1 (measured 0.93, 0.95, 0.97, 0.985), "
                       "estimate at order 3", chain="adams:3 dirk:4"),
        estimator_case("G10-estimate-dirk4-bdf4-adams4", "dirk4-bdf4-adams4", 4,
                       "chain DIRK-4, BDF-4, Adams-4 by BDF-5, BDF-5, Adams-5: effectivity at "
                       "order 1 (measured 1.21, 1.07, 1.03, 1.015), estimate at order 4",
                       chain="dirk:4 bdf:4 adams:4"),
        estimator_case("G11-estimate-dirk3-adams3", "dirk3-adams3", 3,
                       "chain DIRK-3 then Adams-3 by BDF-4 then Adams-4: effectivity at order "
                       "1 (measured 0.61, 0.78, 0.88, 0.94), estimate at order 3",
                       chain="dirk:3 adams:3"),
        estimator_case("G12-estimate-bdf3-dirk3", "bdf3-dirk3", 3,
                       "chain BDF-3 then DIRK-3 by BDF-4 then BDF-4: effectivity at order 1 "
                       "(measured 0.87, 0.95, 0.98, 0.99), estimate at order 3",
                       chain="bdf:3 dirk:3"),
        adaptation_case("A01-adaptive-dirk2", "dirk2", 1.0e-3,
                        "implicit midpoint from 21 instants at the relative tolerance 1e-3: "
                        "E - E_h = 2.5e-3 at h = 0.1 rejects the seed, every step is halved "
                        "once, and the 40-step grid has E - E_h = 6.2e-4 (ratio 1.09)",
                        families="dirk", max_order=2),
        field_estimator_case("G13-field-estimate-dirk3", "dirk3", 3,
                             "enrichment A in time on the 16 x 16 periodic mode: Crouzeix DIRK-3 "
                             "estimated by BDF-4 against the semi-discrete mode energy; "
                             "effectivity at order 1 (measured 0.38, 0.65, 0.81, 0.90 on 6, 11, "
                             "21, 41 instants, slopes 0.82, 0.93, 0.95), estimate at order 3"),
        derivative_estimator_case("G14-derivative-bdf3-energy", "bdf3", "energy", 3,
                                  "the error of the design derivative dE/dnu of BDF-3, "
                                  "estimated through the enriched costate rate of the "
                                  "order-1 Lagrangian: effectivity at order 1 (measured "
                                  "0.951, 0.985, 0.995, 0.998, 0.9991 on 21 to 321 "
                                  "instants, slopes 1.73, 1.53, 1.35, 1.22 falling to 1 "
                                  "from above), estimate at order 3"),
        derivative_estimator_case("G15-derivative-bdf3-dissipation", "bdf3", "dissipation", 3,
                                  "the error of dD/dnu of BDF-3 by the same identity: "
                                  "effectivity at order 1 (measured 0.832, 0.922, 0.963, "
                                  "0.982, 0.991, slopes 1.11, 1.06, 1.03, 1.015), estimate "
                                  "at order 3"),
        mean_case("G16-mean-bdf2", "bdf2", 2,
                  "the near-zero functional mean = int q dt over one period, exact value 0 "
                  "by cancellation: BDF-2 estimated by BDF-3 gives eta at order 2 (measured "
                  "slopes 1.75, 1.94, 1.99) and the scale S at order 2 against int |q| = 4 "
                  "(measured 3.877, 3.970, 3.993, 3.998, slopes 2.02, 2.07, 2.05); its "
                  "effectivity 0.902, 0.978, 0.994, 0.9986 falls at order 2, above the "
                  "order 1 the class declares, and is not declared here"),
        dict(estimator_case("G17-zero-scale-bdf3", "bdf3", 3,
                            "the dissipation functional at nu = 0: its integrand "
                            "nu (1 - q^2) q'^2 vanishes identically, so the scale S of the "
                            "relative criterion is zero, the costate right side is zero and "
                            "the estimate is exactly zero on every grid - the criterion is "
                            "met without dividing by a functional value",
                            effectivity=False),
             checks=[floor_check("estimate:dissipation", 0.0, TABLE_DIGITS,
                                 "a functional whose integrand is identically zero has "
                                 "F = F_h = 0, lambda+ = 0 and both parts of eta exactly "
                                 "zero: within the print resolution plus gamma_N", scale=1.0),
                     floor_check("scale:dissipation", 0.0, TABLE_DIGITS,
                                 "the scale S = sum |w_kj f(Q_j)| of an identically zero "
                                 "integrand is exactly zero", scale=1.0),
                     residual_check()]),
        adaptation_case("A02-adaptive-bdf3", "bdf3", 1.0e-4,
                        "BDF-3 from 21 instants at the relative tolerance 1e-4: eta = -3.3e-4 "
                        "rejects the seed; its indicators vary along t, so the marking divides "
                        "the steps unequally, and the accepted grid has E - E_h = -8.3e-5 "
                        "within 1e-4 x S = 1.17e-4 (ratio 0.61 on this non-uniform grid)",
                        families="bdf", max_order=3),
        # THE RADIAL OSCILLATOR, a law absent from the tested baseline: the
        # residual is nonlinear through a negative integer power and its two
        # functionals have closed forms, the square integral's depending on
        # the trajectory in both its value and its design derivative.
        radial_case("E01-radial-bdf2", "bdf2", 2, ["E", "dE", "F2", "q", "qd"],
                    "radial oscillator q'' + q - nu/q^3 at nu = 2, BDF-2, order 2; its "
                    "dF_2/dnu changes sign inside these grids (X22) and is not declared"),
        radial_case("E02-radial-bdf3", "bdf3", 3, ["E", "dE", "F2", "dF2", "q", "qd"],
                    "radial oscillator, BDF-3, order 3 in both functionals, both design "
                    "derivatives and the state"),
        radial_case("E03-radial-bdf4", "bdf4", 4, ["E", "dE", "q", "qd"],
                    "radial oscillator, BDF-4, order 4 in the energy, its design derivative "
                    "and the state; the square integral and its design derivative carry the "
                    "staged startup's lower order (X23) and are not declared"),
        radial_case("E04-radial-adams2", "adams2", 2,
                    ["E", "dE", "F2", "dF2", "q", "qd"],
                    "radial oscillator, Adams-Moulton 2 (trapezoidal), order 2"),
        radial_case("E05-radial-adams3", "adams3", 3,
                    ["E", "dE", "F2", "dF2", "q", "qd"],
                    "radial oscillator, Adams-Moulton 3, order 3"),
        radial_case("E06-radial-dirk2", "dirk2", 2, ["E", "dE", "F2", "dF2", "q", "qd"],
                    "radial oscillator, implicit midpoint, order 2"),
        radial_case("E07-radial-dirk3", "dirk3", 3, ["E", "dE", "F2", "dF2", "q", "qd"],
                    "radial oscillator, Crouzeix two-stage DIRK, order 3"),
        radial_case("E08-radial-dirk4", "dirk4", 4, ["E", "dE", "F2", "dF2", "q", "qd"],
                    "radial oscillator, Crouzeix three-stage DIRK, order 4"),
        radial_case("E09-radial-newmark2", "newmark2", 2,
                    ["E", "dE", "F2", "dF2", "q", "qd"],
                    "radial oscillator, Newmark beta = 1/4, gamma = 1/2, order 2: the same "
                    "trapezoidal step as Adams-Moulton 2, and the same numbers",
                    argv_extra=("--families=newmark",)),
        # CONSERVATION. The radial oscillator's energy
        # E = q'^2/2 + q^2/2 + nu/(2 q^2) is (1 + nu)/2 at every instant of the
        # continuous flow. None of these families is symplectic, so the defect
        # drifts with the horizon and no exact conservation is claimed; at a
        # fixed horizon it converges at the order of the scheme.
        radial_case("G18-radial-dirk2-invariant", "dirk2", 2, ["radial_invariant"],
                    "the invariant defect of the implicit midpoint at T = 2, order 2"),
        radial_case("G19-radial-dirk3-invariant", "dirk3", 3, ["radial_invariant"],
                    "the invariant defect of the Crouzeix two-stage DIRK, order 3"),
        radial_case("G20-radial-bdf3-invariant", "bdf3", 3, ["radial_invariant"],
                    "the invariant defect of BDF-3, order 3"),
        radial_case("G21-radial-adams2-invariant", "adams2", 2, ["radial_invariant"],
                    "the invariant defect of Adams-Moulton 2, order 2"),
        # THE DISC, a second discretization use on supported geometry: the polar
        # mesh's identified angular seam, curved Neumann boundary, anisotropic
        # cells and polygonal centre cell.
        disc_operator_case("D01-disc-operator-centre", "operator_centre",
                           "the polygonal centre cell of the disc, its own class, at order 2"),
        disc_conservation_case(),
        disc_mode_case("D05-disc-mode-fitted-balance",
                       "the radial Neumann mode marched on the disc, the spatial law "
                       "substituted into the state row as the fitted balance: order 2"),
        disc_transpose_case("D06-disc-transpose",
                            "the disc's forward and reverse passes evaluate one bilinear "
                            "form, the fitted balance"),
        disc_transpose_case("D08-disc-transpose-jet",
                            "the same bilinear form with the spatial derivatives as rows of "
                            "the jet, the second discretization of the same law",
                            rows="states state-time-derivatives state-spatial-derivatives"),
        radial_case("E10-radial-alexander2", "alexander2", 2,
                    ["E", "dE", "F2", "dF2", "q", "qd"],
                    "radial oscillator, Alexander's two-stage L-stable DIRK, order 2: a "
                    "tableau registered by its data alone",
                    argv_extra=("--families=alexander",)),
        # THE FUNCTIONAL DISCRETIZATION-ERROR ESTIMATOR ON THE NEW LAW. The
        # radial oscillator is the one law stated here whose functionals are
        # known in closed form at every grid, so the effectivity is read
        # against an exact F rather than against a refined one.
        radial_estimator_case("G22-radial-estimate-bdf1", "bdf1", 1,
                              "the same estimate of BDF-1 by BDF-2, the pair whose effectivity "
                              "is furthest from 1 on the coarse grids: the energy measured "
                              "0.7098, 0.8495, 0.9243, 0.9622, 0.9811, slopes 0.95, 0.99, "
                              "1.00, 1.00, and the square integral 0.7494, 0.8763, 0.9390, "
                              "0.9698, 0.9850, slopes 1.02, 1.02, 1.01, 1.01; each estimate "
                              "at order 1"),
        radial_estimator_case("G23-radial-estimate-bdf3", "bdf3", 3,
                              "the functional-error estimate on the radial oscillator: BDF-3 "
                              "estimated by BDF-4 on the same grid, each functional against "
                              "its closed form; effectivity at order 1 (the energy measured "
                              "0.8868, 0.9653, 0.9880, 0.9953, slopes 1.70, 1.54, 1.35, and "
                              "the square integral 0.8857, 0.9601, 0.9844, 0.9933, slopes "
                              "1.52, 1.36, 1.21, both falling to 1 from above), each estimate "
                              "at order 3"),
    ]
    # declared limitations: measured below their theoretical order
    cases += [
        # THE SECOND DISCRETIZATION OF THE SAME LAW ON THE DISC. The spatial
        # derivatives are rows of the jet, each tied to the values by the fit's
        # row at form degree 2. That form is the compact one - the powers of one
        # coordinate over the cell and its face neighbours, with no mixed
        # member - and on the polar mesh those neighbours lie along the radial
        # and angular directions, which are the coordinate axes only along two
        # rays. On the box, where they are always axis aligned, the same form
        # reaches order 2 (config/mode_jet.cfg).
        disc_mode_case("D07-disc-mode-jet",
                       "the radial Neumann mode marched on the disc with the spatial "
                       "derivatives as rows of the jet, declared at order 2",
                       rows="states state-time-derivatives state-spatial-derivatives",
                       limitation="the compact form fits the second derivatives from "
                                  "axis-pure members over a neighbourhood the polar mesh "
                                  "does not align with the axes: errors 3.778e-2, 2.190e-2, "
                                  "1.703e-2 at the step ratio 1.9412, slopes 0.822 and 0.379, "
                                  "against the fitted balance's 4.791e-3, 1.162e-3, 3.035e-4 "
                                  "at 2.136 and 2.024"),
        disc_operator_case("D02-disc-operator-interior", "operator",
                           "the interior cells of the disc declared at order 2",
                           limitation="the fitted balance on the polar mesh converges at "
                                      "order 1.74 in the interior: errors 3.653e-2, 1.187e-2, "
                                      "3.734e-3, slopes 1.768 and 1.744 at the step ratio "
                                      "1.9412, and 1.208e-3 over a fourth grid of 65 x 128"),
        disc_operator_case("D03-disc-operator-ring", "operator_boundary",
                           "the boundary ring of the disc declared at order 2",
                           limitation="the one-sided fits against the curved Neumann boundary "
                                      "converge at order 1.33: errors 2.327e-1, 8.220e-2, "
                                      "3.395e-2, slopes 1.636 and 1.333 at the step ratio "
                                      "1.9412, and 1.530e-2 over a fourth grid of 65 x 128"),
        spatial_case("L04-operator-degree4", "bdf1", "operator", 4,
                     "discrete Laplacian of the mode at form degree 4 converges at order 2",
                     instants=3, families="bdf", max_order=1, extra=(), spatial_order=4,
                     limitation="interior relative rms error measured at slopes 2.5, 2.2"),
    ]
    return cases


def exploratory_cases():
    numeric_reference = temporal_case("X01-dirk2-numeric-reference", "dirk2", 2, ["dE"],
                                      "dE/dnu of DIRK-2 against a grid four times finer than "
                                      "the finest: the reference bias is admitted when below "
                                      "the window", instants=(21, 41, 81),
                                      argv_extra=("--families=dirk",))
    numeric_reference["runs"]["reference"] = ode_argv(321, extra=("--families=dirk",))
    numeric_reference["checks"][0]["reference"] = {"run": "reference", "refinement": 4.0}
    numeric_reference["checks"] = numeric_reference["checks"][:1]
    cases = [
        numeric_reference,
        temporal_case("X02-bdf2-energy", "bdf2", 3, ["E"],
                      "BDF-2 energy functional on the linear oscillator: order 3 measured"),
        temporal_case("X12-adams4-energy", "adams4", 5, ["E"],
                      "Adams-Moulton 4 energy functional on the linear oscillator: order 5 "
                      "over the three grids above the 12-digit print floor (the error at 161 "
                      "instants is one unit of the last printed digit)", instants=(21, 41, 81)),
        temporal_case("X08-newmark3-position", "newmark3", 4, ["q"],
                      "Newmark beta = 1/12 (Fox-Goodwin): fourth-order position on the "
                      "linear oscillator", argv_extra=("--families=newmark",)),
        temporal_case("X09-dirk4-invariant-drift", "dirk4", 5, ["invariant"],
                      "DIRK-4 invariant drift at order 5"),
        temporal_case("X10-bdf2-invariant-drift", "bdf2", 3, ["invariant"],
                      "BDF-2 invariant drift at order 3"),
        temporal_case("X11-bdf4-invariant-drift", "bdf4", 5, ["invariant"],
                      "BDF-4 invariant drift at order 5"),
        temporal_case("X04-bdf1-dE-crossing", "bdf1", 1, ["dE"],
                      "BDF-1 dE/dnu changes sign inside these grids"),
        estimator_case("X13-estimate-bdf2-heuristic", "bdf2", 3,
                       "BDF-2 energy (order 3 by superconvergence) estimated by BDF-3 (order "
                       "3): p+ = p, the estimate is a heuristic; its effectivity tends to 2 "
                       "(measured 1.89, 1.95, 1.98, 1.99) and is not declared; the estimate "
                       "itself converges at order 3", effectivity=False),
        estimator_case("X14-estimate-bdf4-heuristic", "bdf4", 5,
                       "BDF-4 energy (order 5 by superconvergence) estimated by BDF-5: p+ = p, "
                       "a heuristic; effectivity measured 1.31, 1.41, 1.46, 1.52 and not "
                       "declared; the estimate converges at order 5", effectivity=False),
        estimator_case("X15-estimate-newmark1", "newmark1", 1,
                       "Newmark beta = gamma = 0 (the explicit Taylor step, order 1) estimated "
                       "by Adams-Moulton 3: effectivity at order 1 (measured 1.02, 1.01, "
                       "1.006, 1.003)", families="newmark", max_order=3),
        field_estimator_case("X19-field-estimate-bdf3", "bdf3", 3,
                             "BDF-3 on the 16 x 16 mode by BDF-4: effectivity 0.37, 0.72, 0.88, "
                             "0.94 (order 1, slopes 1.15, 1.19, 1.10); the estimate is "
                             "declared at the family's order 3 and measured 0.99, 2.29, 2.71 "
                             "below it on these grids (omega_h h from 0.75 to 0.09), the "
                             "recorded shortfall", set_name="exploratory"),
        dict(field_estimator_case("X20-field-estimate-dirk2", "dirk2", 2,
                                  "implicit midpoint on the 16 x 16 mode by BDF-3: effectivity "
                                  "0.79, 0.94, 0.99, 1.0006 converges faster than O(h) and "
                                  "crosses 1 at the finest grid, so no order of |I - 1| is "
                                  "declared; the estimate at order 2", set_name="exploratory",
                                  effectivity=False)),
        dict(field_estimator_case("X21-field-estimate-bdf1", "bdf1", 1,
                                  "BDF-1 on the 16 x 16 mode by BDF-2: effectivity 0.15, 0.45, "
                                  "0.68, 0.83 (order 1, slopes 0.63, 0.79, 0.89); its estimate "
                                  "is not monotone on the two coarsest grids and is not declared",
                                  set_name="exploratory"), checks=[
                 dict(order_check("effectivity:energy", 1, 1.0, TABLE_DIGITS, EFFECTIVITY_TEXT),
                      functional_reference="semi_energy"), TRANSFER, residual_check()]),
        estimator_case("X16-estimate-alexander2", "alexander2", 2,
                       "Alexander's L-stable DIRK estimated by BDF-3: effectivity 1.61, 1.36, "
                       "1.19, 1.10 (order 1 in |I - 1|, far from 1 on the coarse grids), "
                       "estimate at order 2", families="alexander", max_order=2),
        estimator_case("X17-estimate-dirk3-bdf2-heuristic", "dirk3-bdf2", 3,
                       "chain DIRK-3 then BDF-2: the BDF-2 energy is superconvergent (order "
                       "3) and BDF-3 enriches it at the same order, a heuristic; effectivity "
                       "measured 1.05, 1.26, 1.37, 1.43 and not declared; estimate at order 3",
                       chain="dirk:3 bdf:2", effectivity=False),
        estimator_case("X18-estimate-bdf2-dirk3-heuristic", "bdf2-dirk3", 3,
                       "chain BDF-2 then DIRK-3, a heuristic for the same reason; effectivity "
                       "measured 1.60, 1.76, 1.83, 1.86 and not declared; estimate at order 3",
                       chain="bdf:2 dirk:3", effectivity=False),
        radial_case("X22-radial-bdf2-design-derivative", "bdf2", 2, ["dF2"],
                    "radial oscillator, BDF-2: dF_2/dnu changes sign between 21 and 41 "
                    "instants, the pairwise slopes 4.462, 0.595, 1.134, 1.687 contracting "
                    "towards 2", set_name="exploratory"),
        radial_case("X23-radial-bdf4-square-integral", "bdf4", 4, ["F2", "dF2"],
                    "radial oscillator, BDF-4: the square integral converges at the measured "
                    "2.72 and its design derivative changes sign between 161 and 321 instants",
                    set_name="exploratory"),
        field_temporal_case("X05-field-bdf3", "bdf3", 3,
                            "temporal order on the 16 x 16 periodic mode, BDF-3"),
        spatial_case("X06-taylor-green-three-grids", "bdf2", "velocity", 2,
                     "Taylor-Green vortex, BDF-2, cells 8, 16, 32 (about a minute at 32)",
                     cells=(8, 16, 32), instants=6, families="bdf", max_order=2, extra=(),
                     timeout=600),
    ]
    cases[-1]["runs"] = {f"cells={n}": taylor_green_argv(n, 6) for n in (8, 16, 32)}
    cases[-1]["grids"] = [(f"cells={n}", 2 * math.pi / n, n * n * 36) for n in (8, 16, 32)]
    demo = {"id": "X07-order-demo-reduced", "set": "exploratory", "record": "order_demo",
            "description": "the reduced order demonstration's own records: every ORDER_RECORD "
                           "status is reported, unresolved ones included",
            "row": None, "grids": [("demo", None, 0)],
            "runs": {"demo": [APPLICATION, "--demo=order_of_accuracy", "--max_derivative_degree=0",
                              "--coarsest_instants=6", "--refinement_grids=3",
                              "--reference_refinement=1.5", "--time_duration=1"]},
            "checks": [{"kind": "records", "quantity": "order_records",
                        "justification": "a record is exploratory unless every status reaches"}],
            "limitation": None, "timeout": 300}
    cases.append(demo)
    return cases


def rejection_cases(fixtures):
    """Cases that must fail, each with the status the contract must report."""
    dirk = ("--families=dirk",)
    wrong = temporal_case("R01-wrong-order-above", "dirk2", 3, ["E"],
                          "DIRK-2 declared order 3", instants=(21, 41, 81, 161),
                          argv_extra=dirk)
    wrong["expected_status"] = "below"
    low = temporal_case("R02-wrong-order-below", "dirk3", 2, ["E"],
                        "DIRK-3 declared order 2", instants=(21, 41, 81, 161),
                        argv_extra=dirk)
    low["expected_status"] = "exceeds"
    under = temporal_case("R03-under-resolved-reference", "dirk2", 2, ["E"],
                          "reference grid only 1.5 times finer than the finest",
                          instants=(21, 41, 81), argv_extra=dirk)
    under["runs"]["reference"] = ode_argv(121, extra=dirk)
    under["checks"][0]["reference"] = {"run": "reference", "refinement": 1.5}
    under["expected_status"] = "under_resolved_reference"
    malformed = temporal_case("R04-malformed-record", "dirk2", 2, ["E"],
                              "a record whose value column is not a number", instants=(21, 41))
    malformed["runs"] = {"instants=21": ["cat", str(fixtures / "malformed.out")],
                         "instants=41": ["cat", str(fixtures / "malformed.out")]}
    malformed["expected_status"] = "malformed_record"
    nonfinite = temporal_case("R05-nonfinite-value", "dirk2", 2, ["E"],
                              "a record carrying NaN", instants=(21, 41))
    nonfinite["runs"] = {"instants=21": ["cat", str(fixtures / "nonfinite.out")],
                         "instants=41": ["cat", str(fixtures / "nonfinite.out")]}
    nonfinite["expected_status"] = "nonfinite"
    unconverged = temporal_case("R06-solver-failure", "bdf2", 2, ["dE"],
                                "one Newton iteration on van der Pol at nu = 1 leaves the "
                                "rows unconverged", instants=(21, 41),
                                argv_extra=("--design=1.0", "--initial_state=2.0",
                                            "--time_duration=7.0", "--iteration_criterion=by_count",
                                            "--max_iterations=1"))
    unconverged["expected_status"] = "solver_failure"
    slow = spatial_case("R07-timeout", "dirk2", "mode", 2,
                        "the dense direct solver at 32 x 32 exceeds a one second limit",
                        cells=(16, 32), families="dirk", max_order=2, extra=(), timeout=1)
    slow["expected_status"] = "timeout"
    missing = temporal_case("R08-missing-row", "bdf2", 2, ["E"],
                            "the declared row is absent from the record", instants=(21, 41),
                            argv_extra=("--families=dirk",))
    missing["expected_status"] = "missing_result"
    column = temporal_case("R09-missing-column", "dirk2", 2, ["dE"],
                           "the derivative column is not computed", instants=(21, 41),
                           argv_extra=("--max_derivative_degree=0", "--families=dirk"))
    column["expected_status"] = "missing_result"
    process = temporal_case("R10-process-failure", "dirk2", 2, ["E"],
                            "an unknown setting stops the program", instants=(21, 41),
                            argv_extra=("--no_such_setting=1",))
    process["expected_status"] = "process_failure"
    residual = temporal_case("R11-residual-unmet", "dirk2", 2, ["E"],
                             "a stopping tolerance the solved rows do not satisfy",
                             instants=(21, 41), argv_extra=dirk)
    residual["checks"] = [residual_check(tolerance=1.0e-20)]
    residual["expected_status"] = "residual_unmet"
    roundoff = temporal_case("R12-roundoff-unresolved", "bdf4", 5, ["E"],
                             "BDF-4 energy at the print resolution: no slope can be read",
                             instants=(81, 161, 321), argv_extra=("--families=bdf",))
    roundoff["expected_status"] = "unresolved"
    stale = temporal_case("R13-stale-limitation", "dirk2", 2, ["E"],
                          "a declared limitation that is met must be reported",
                          instants=(21, 41, 81, 161), limitation="declared for this test only",
                          argv_extra=dirk)
    stale["expected_status"] = "unexpected_pass"
    unmet = adaptation_case("R14-adaptation-limit", "dirk2", 1.0e-3,
                             "the adaptive loop with adaptation_instants = 30 below the 41 "
                             "instants the tolerance needs: no acceptance, the process reports "
                             "ADAPTATION_UNMET and fails", families="dirk", max_order=2, limit=30)
    unmet["expected_status"] = "process_failure"
    return [wrong, low, under, malformed, nonfinite, unconverged, slow, missing, column,
            process, residual, roundoff, stale, unmet]
