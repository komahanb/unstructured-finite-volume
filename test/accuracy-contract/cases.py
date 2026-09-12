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

# the configured Newton stopping rule of every run: relative 1e-12
TOLERANCE = 1.0e-12
STATE_COMPONENTS = 3

THETA_TEXT = ("window from the two-term expansion e = C h^p (1 + a h) with |a h| <= 1/4 "
              "on the coarser grid of the finest pair; roundoff floor = print resolution "
              "+ gamma_N |reference|")


def ode_argv(instants, families="bdf adams dirk", max_order=4, derivative=1,
             physics="vanderpol", chain=None, extra=()):
    argv = [APPLICATION, f"--physics={physics}", "--design=0.0", "--initial_state=1.0",
            f"--time_duration={DURATION}", f"--instants={instants}", "--grid=uniform",
            f"--families={families}", f"--max_discretization_order={max_order}",
            f"--max_derivative_degree={derivative}", "--functionals=energy dissipation",
            "--combinations=1", "--check=state passes"]
    if chain:
        argv.append(f"--chain={chain}")
    return argv + list(extra)


def ode_grids(instants=INSTANTS):
    return [(f"instants={n}", DURATION / (n - 1), n * STATE_COMPONENTS) for n in instants]


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


def order_check(quantity, order, reference, digits, justification):
    return {"kind": "order", "quantity": quantity, "order": order,
            "reference": reference, "digits": digits, "justification": justification}


def floor_check(quantity, reference, digits, justification, scale=None):
    return {"kind": "floor", "quantity": quantity, "reference": reference,
            "digits": digits, "scale": scale, "justification": justification}


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
    ]
    # declared limitations: measured below their theoretical order
    cases += [
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
    return [wrong, low, under, malformed, nonfinite, unconverged, slow, missing, column,
            process, residual, roundoff, stale]
