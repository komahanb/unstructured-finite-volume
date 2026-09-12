#!/usr/bin/env python3
"""The accuracy contract: run the declared cases, read their records, measure
orders and floors, write a machine-readable summary and return failure when a
declaration is not met.

    contract.py --set required|exploratory|rejection [--case ID ...]
                [--results DIR] [--root DIR]

Check kinds: order, floor, residual, records, value (a quantity recorded
without a declaration, failing only where it is absent or not finite).

Statuses of a case: pass, limitation (a declared limitation that is still
unmet), and the failures below, unresolved, exceeds, not_monotone,
under_resolved_reference, missing_result, malformed_record, nonfinite,
floor_exceeded, solver_failure, residual_unmet, process_failure, timeout,
unexpected_pass.
"""

import argparse
import datetime
import hashlib
import json
import math
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

import cases as declared

UNIT_ROUNDOFF = 2.0 ** -53
THETA = 0.25
DECLARATION_FAILURES = {"below", "exceeds", "unresolved", "not_monotone", "floor_exceeded"}
ACCEPTED = {"pass", "limitation"}


def gamma(terms):
    """First-order bound on the relative error of a sum of `terms` terms."""
    n_u = terms * UNIT_ROUNDOFF
    return n_u / (1.0 - n_u)


def window(ratio, theta=THETA):
    """|p_obs - p| for e = C h^p (1 + a h), |a h| <= theta on the coarser grid."""
    return -math.log((1.0 - theta) / (1.0 - theta / ratio)) / math.log(ratio)


def reference_bias(order, ratio, refinement):
    """Upward shift of the finest slope when the reference is a grid refined
    `refinement` times beyond the finest measured grid."""
    x = refinement ** (-order)
    return math.log((ratio ** order - x) / (1.0 - x)) / math.log(ratio) - order


def print_resolution(value, digits):
    """Half a unit in the last printed significant digit."""
    if value == 0.0 or not math.isfinite(value):
        return 0.0
    return 0.5 * 10.0 ** (math.floor(math.log10(abs(value))) - (digits - 1))


# ----------------------------------------------------------------- records

ROW = re.compile(r"^  (\S+)\s+(\d+)\s+(.*?)\s*$")
FLOAT = re.compile(r"^[-+]?(\d+\.?\d*([eE][-+]?\d+)?|NaN|Infinity|-Infinity)$")
# Fortran ES editing drops the E before a three-digit exponent: 1.0-165
WIDE_EXPONENT = re.compile(r"^([-+]?\d+\.\d+)([-+]\d{3})$")


class MalformedRecord(Exception):
    pass


def number(token):
    if token == "-":
        return None
    if FLOAT.match(token):
        return float(token)
    wide = WIDE_EXPONENT.match(token)
    if wide:
        return float(wide.group(1) + "e" + wide.group(2))
    raise MalformedRecord(f"not a number: {token!r}")


def parse_table(text):
    """Rows of the application table with their check lines."""
    rows = {}
    record = {"rows": rows, "operator": None}
    current = None
    for line in text.splitlines():
        m = re.match(r"^   relative rms error\s+interior\s+(\S+)\s+one boundary\s+(\S+)\s+corner\s+(\S+)", line)
        if m:
            record["operator"] = {"interior": number(m.group(1)), "one_boundary": number(m.group(2)),
                                  "corner": number(m.group(3))}
            continue
        m = ROW.match(line)
        if m and not line.startswith("      "):
            tokens = m.group(3).split()
            status = None
            if tokens and tokens[-1] in ("unconverged", "diverging"):
                status = tokens.pop()
            current = {"label": m.group(1), "solved": int(m.group(2)),
                       "f": [number(t) for t in tokens], "status": status,
                       "dissipation": None, "state": None, "blocks": [], "transpose": None,
                       "mode": None, "semi": None, "velocity": None, "pressure": None,
                       "divergence": None, "functional_error": {}, "indicators": {}}
            rows[current["label"]] = current
            continue
        if current is None:
            continue
        m = re.match(r"^      van der pol dissipation\s+(.*?)\s*$", line)
        if m:
            current["dissipation"] = [number(t) for t in m.group(1).split()]
            continue
        m = re.match(r"^      state at the last instant, node 1:\s*(.*?)\s*$", line)
        if m:
            current["state"] = [number(t) for t in m.group(1).split()]
            continue
        m = re.match(r"^      block (\d+): imbalance\s+(\S+) of initial\s+(\S+) converged ([TF])", line)
        if m:
            current["blocks"].append({"block": int(m.group(1)), "norm": number(m.group(2)),
                                      "initial": number(m.group(3)), "converged": m.group(4) == "T"})
            continue
        m = re.match(r"^      functional error, (\S+): estimate\s+(\S+) residual\s+(\S+) quadrature\s+(\S+)"
                     r" scale\s+(\S+) transfer\s+(\S+) enriched (\S+)", line)
        if m:
            current["functional_error"][m.group(1)] = {
                "estimate": number(m.group(2)), "residual": number(m.group(3)),
                "quadrature": number(m.group(4)), "scale": number(m.group(5)),
                "transfer": number(m.group(6)), "enriched": m.group(7)}
            continue
        m = re.match(r"^      indicator, (\S+), step (\d+): t\s+(\S+) h\s+(\S+) eta\s+(\S+)", line)
        if m:
            current["indicators"].setdefault(m.group(1), []).append(
                {"step": int(m.group(2)), "t": number(m.group(3)), "h": number(m.group(4)),
                 "eta": number(m.group(5))})
            continue
        m = re.match(r"^      tangent against adjoint over the table, relative\s+(\S+)", line)
        if m:
            current["transpose"] = number(m.group(1))
            continue
        m = re.match(r"^      error at the last instant, compared with the mode\s+(\S+)\s+semi-discrete\s+(\S+)", line)
        if m:
            current["mode"] = number(m.group(1))
            current["semi"] = number(m.group(2))
            continue
        m = re.match(r"^      taylor-green at the last instant: velocity error, relative rms\s+(\S+)\s+"
                     r"pressure error, mean removed\s+(\S+)\s+divergence rms\s+(\S+)", line)
        if m:
            current["velocity"] = number(m.group(1))
            current["pressure"] = number(m.group(2))
            current["divergence"] = number(m.group(3))
    return record


def parse_sensitivity(text):
    rows = {}
    current = None
    for line in text.splitlines():
        m = re.match(r"^ (\S.*), van der pol at a design of one$", line)
        if m:
            current = {"label": m.group(1), "status": None, "blocks": []}
            rows[current["label"]] = current
            continue
        if current is None:
            continue
        for key, prefix in (("functional", "the functional"), ("tangent", "sensitivity, tangent"),
                            ("adjoint", "sensitivity, adjoint"), ("difference", "sensitivity, differenced"),
                            ("tangent_adjoint", "tangent against adjoint"),
                            ("tangent_difference", "tangent against difference")):
            m = re.match(r"^   " + re.escape(prefix) + r"\s+(\S+)\s*$", line)
            if m:
                current[key] = number(m.group(1))
    return {"rows": rows, "operator": None}


def parse_order_demo(text):
    records = []
    summary = None
    for line in text.splitlines():
        m = re.match(r'^ ORDER_RECORD chain="(.*)" p=(\d+)(?: r=(\d+))?(?: fitted=\s*(\S+) spread=\s*(\S+))? status=(\S+)', line)
        if m:
            records.append({"chain": m.group(1), "p": int(m.group(2)),
                            "r": int(m.group(3)) if m.group(3) else None,
                            "fitted": number(m.group(4)) if m.group(4) else None,
                            "spread": number(m.group(5)) if m.group(5) else None,
                            "status": m.group(6)})
        m = re.match(r"^ ORDER_SUMMARY (.*)$", line)
        if m:
            summary = dict(item.split("=") for item in m.group(1).split())
    if summary is None:
        raise MalformedRecord("no ORDER_SUMMARY line")
    return {"rows": {"demo": {"label": "demo", "status": None, "blocks": [], "records": records,
                              "summary": summary}}, "operator": None}


PARSERS = {"table": parse_table, "sensitivity": parse_sensitivity, "order_demo": parse_order_demo}


def functional_value(row, word):
    """The value F_h of the functional a configured word names, or None."""
    if word == "energy":
        return row["f"][0] if row.get("f") else None
    if word == "dissipation":
        return row["dissipation"][0] if row.get("dissipation") else None
    values = row.get(word)
    return values[0] if values else None


def interval_sum(indicators, a, b):
    """The sum of |eta_k| over the steps inside [a, b]; an instant within half
    a step of a bound lies on it (the application's coarsening rule)."""
    if not indicators:
        return None
    return sum(abs(i["eta"]) for i in indicators
               if i["t"] - i["h"] >= a - i["h"] / 2 and i["t"] <= b + i["h"] / 2)


def quantity_of(row, record, quantity, check=None, paired=None):
    """The number a quantity names in one row, or None when it is absent.
    The functional-error quantities read `functional error` and `indicator`
    lines: estimate:<word>, scale:<word>, transfer:<word>, effectivity:<word>
    = eta / (F - F_h) with F the declared reference of the word (infinite
    where F_h = F), and localization:<word> = log2 of the sum of |eta_k|
    over the check's interval against the same sum on the paired row."""
    if ":" in quantity:
        kind, word = quantity.split(":", 1)
        estimate = row.get("functional_error", {}).get(word)
        if estimate is None:
            return None
        if kind in ("estimate", "scale", "transfer"):
            return estimate[kind]
        if kind == "effectivity":
            value = functional_value(row, word)
            if value is None:
                return None
            error = declared.FUNCTIONAL_REFERENCES[word] - value
            return estimate["estimate"] / error if error != 0.0 else math.inf
        if kind == "localization":
            if paired is None or check is None:
                return None
            a, b = check["interval"]
            own = interval_sum(row["indicators"].get(word), a, b)
            other = interval_sum(paired["indicators"].get(word), a, b)
            if own is None or other is None or own <= 0.0 or other <= 0.0:
                return None
            return math.log2(own / other)
        raise ValueError(quantity)
    if quantity == "E":
        return row["f"][0] if row.get("f") else None
    if quantity == "dE":
        return row["f"][1] if row.get("f") and len(row["f"]) > 1 else None
    if quantity == "dD":
        return row["dissipation"][1] if row.get("dissipation") and len(row["dissipation"]) > 1 else None
    if quantity in ("q", "qd", "qdd"):
        index = ("q", "qd", "qdd").index(quantity)
        state = row.get("state")
        return state[index] if state and len(state) > index else None
    if quantity == "invariant":
        state = row.get("state")
        if not state or len(state) < 2 or state[0] is None or state[1] is None:
            return None
        return 0.5 * (state[0] ** 2 + state[1] ** 2)
    if quantity == "law":
        state = row.get("state")
        if not state or len(state) < 3 or state[0] is None or state[2] is None:
            return None
        return state[2] + state[0]
    if quantity == "operator":
        return record["operator"]["interior"] if record.get("operator") else None
    return row.get(quantity)


# ----------------------------------------------------------------- execution

class Executor:
    def __init__(self, root, results):
        self.root = root
        self.results = results
        self.completed = {}
        (results / "runs").mkdir(parents=True, exist_ok=True)

    def run(self, case_id, label, argv, timeout):
        key = tuple(argv)
        if key in self.completed:
            return self.completed[key]
        digest = hashlib.sha256(" ".join(argv).encode()).hexdigest()[:12]
        log = self.results / "runs" / f"{case_id}-{label}-{digest}.out"
        command = list(argv)
        directory = self.root
        if command[0] == declared.APPLICATION:
            # the application reads config/*.cfg relative to its own directory
            directory = self.root / Path(declared.APPLICATION).parent
            command[0] = "./" + Path(declared.APPLICATION).name
        start = time.perf_counter()
        outcome = {"label": label, "argv": argv, "log": str(log), "timed_out": False,
                   "exit_status": None, "elapsed_seconds": None, "text": ""}
        try:
            completed = subprocess.run(command, cwd=directory, capture_output=True, text=True,
                                       timeout=timeout, check=False)
            outcome["exit_status"] = completed.returncode
            outcome["text"] = completed.stdout
            log.write_text(completed.stdout + ("\n--- stderr ---\n" + completed.stderr
                                               if completed.stderr else ""))
        except subprocess.TimeoutExpired as expired:
            outcome["timed_out"] = True
            outcome["text"] = expired.stdout.decode() if isinstance(expired.stdout, bytes) else (expired.stdout or "")
            log.write_text(outcome["text"] + f"\n--- timeout after {timeout} s ---\n")
        except OSError as failure:
            outcome["exit_status"] = -1
            log.write_text(f"execution failure: {failure}\n")
        outcome["elapsed_seconds"] = time.perf_counter() - start
        self.completed[key] = outcome
        return outcome


# ----------------------------------------------------------------- evaluation

def resolution_of(check, row, value):
    """Half a unit of the last printed digit of a quantity: propagated through
    the quotient for an effectivity, eta / (F - F_h), whose two operands are
    printed; the quantity's own print resolution otherwise."""
    quantity = check["quantity"]
    if quantity.startswith("effectivity:"):
        word = quantity.split(":", 1)[1]
        estimate = row["functional_error"][word]["estimate"]
        f = functional_value(row, word)
        error = abs(declared.FUNCTIONAL_REFERENCES[word] - f)
        if error == 0.0:
            return math.inf
        return (print_resolution(estimate, check["digits"])
                + abs(value) * print_resolution(f, check["digits"])) / error
    return print_resolution(value, check["digits"])


def order_of(check, grids, values, references, terms, scale, resolutions):
    """Status and diagnostics of an order declaration over the grids."""
    p = check["order"]
    result = {"errors": [], "pairwise_slopes": [], "status": None}
    errors = []
    for value, reference in zip(values, references):
        errors.append(abs(value - reference) if reference is not None else value)
    result["errors"] = errors
    if any(e is None or not math.isfinite(e) for e in errors):
        result["status"] = "nonfinite"
        return result
    hs = [h for _, h, _ in grids]
    if any(e == 0.0 for e in errors):
        result["status"] = "unresolved"
        result["message"] = "an error is zero at the printed precision"
        return result
    slopes = [math.log(errors[k] / errors[k + 1]) / math.log(hs[k] / hs[k + 1])
              for k in range(len(errors) - 1)]
    result["pairwise_slopes"] = slopes
    if any(errors[k + 1] >= errors[k] for k in range(len(errors) - 1)):
        result["status"] = "not_monotone"
        result["message"] = "the error does not decrease under refinement"
        return result
    ratio = hs[-2] / hs[-1]
    delta = window(ratio)
    floors = [r + gamma(n) * s for r, n, s in zip(resolutions, terms, scale)]
    widening = (floors[-2] / errors[-2] + floors[-1] / errors[-1]) / math.log(ratio)
    bias = 0.0
    if isinstance(check["reference"], dict):
        bias = reference_bias(p, ratio, check["reference"]["refinement"])
    observed = slopes[-1]
    x = [math.log(h) for h in hs]
    y = [math.log(e) for e in errors]
    mx, my = sum(x) / len(x), sum(y) / len(y)
    least_squares = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sum((a - mx) ** 2 for a in x)
    result.update({"ratio": ratio, "window": delta, "roundoff_widening": widening,
                   "reference_bias": bias, "observed_order": observed,
                   "least_squares_order": least_squares, "floors": floors,
                   "deviation": observed - p,
                   "extrapolated_order": observed + (observed - slopes[-2]) / (ratio - 1.0)
                   if len(slopes) > 1 else None,
                   "contracting": abs(slopes[-1] - p) <= abs(slopes[-2] - p)
                   if len(slopes) > 1 else None,
                   "pairs": len(slopes)})
    if widening > delta:
        result["status"] = "unresolved"
        result["message"] = (f"roundoff widens the slope by {widening:.3f}, beyond the window "
                             f"{delta:.3f}: the finest errors are at the floor")
    elif observed - p < -(delta + widening):
        result["status"] = "below"
        result["message"] = f"observed order {observed:.3f} below {p} - {delta + widening:.3f}"
    elif observed - p > delta + widening + bias:
        result["status"] = "exceeds"
        result["message"] = (f"observed order {observed:.3f} above {p} + {delta + widening + bias:.3f}: "
                             "the declared order is not the leading order")
    else:
        result["status"] = "pass"
        result["message"] = f"observed order {observed:.3f} within {p} -/+ {delta + widening:.3f}"
    return result


def floor_threshold(check, value, row, terms, case):
    scale = check.get("scale")
    if scale == "law":
        initial = max((b["initial"] for b in row["blocks"]), default=0.0)
        q = row["state"][0] if row.get("state") else 0.0
        return declared.TOLERANCE * initial + gamma(terms) * abs(q)
    if scale == "tangent":
        return gamma(terms) * abs(row.get("tangent", 0.0))
    if isinstance(scale, (list, tuple)) and scale[0] == "difference":
        return (declared.TOLERANCE + UNIT_ROUNDOFF) * abs(row.get("functional", 0.0)) / scale[1]
    if isinstance(scale, (list, tuple)) and scale[0] == "estimate":
        # the adaptive criterion: tolerance x S, S the scale of the estimate
        word = check["quantity"] if ":" in check["quantity"] else "energy"
        estimate = row["functional_error"][word.split(":")[-1]]
        return scale[1] * estimate["scale"] + print_resolution(value, check["digits"])
    magnitude = abs(check["reference"]) if scale is None else scale
    return print_resolution(check["reference"] if check["reference"] else value, check["digits"]) \
        + gamma(terms) * magnitude


def evaluate_check(check, case, rows, records):
    """One check over every grid; rows[label] is the case's row in that grid's record."""
    grids = case["grids"]
    kind = check["kind"]
    outcome = {"quantity": check["quantity"], "kind": kind, "justification": check["justification"]}
    if kind == "records":
        statuses = [r["status"] for r in rows[grids[0][0]]["records"]]
        outcome["records"] = rows[grids[0][0]]["records"]
        outcome["summary"] = rows[grids[0][0]]["summary"]
        outcome["status"] = "pass" if statuses and all(s == "reaches" for s in statuses) else "unresolved"
        outcome["message"] = f"{statuses.count('reaches')} of {len(statuses)} records reach p"
        return outcome
    if kind == "residual":
        outcome["tolerance"] = check["tolerance"]
        outcome["blocks"] = {}
        for label, _, _ in grids:
            blocks = rows[label]["blocks"]
            if not blocks:
                outcome["status"] = "missing_result"
                outcome["message"] = f"no block imbalance lines in {label}"
                return outcome
            outcome["blocks"][label] = blocks
            for b in blocks:
                if b["norm"] is None or b["initial"] is None or not math.isfinite(b["norm"]):
                    outcome["status"] = "nonfinite"
                    outcome["message"] = f"block {b['block']} of {label}"
                    return outcome
                limit = check["tolerance"] * b["initial"]
                if b["norm"] > limit:
                    outcome["status"] = "residual_unmet"
                    outcome["message"] = (f"{label} block {b['block']}: imbalance {b['norm']:.3e} above "
                                          f"{check['tolerance']:.1e} x {b['initial']:.3e}")
                    return outcome
        outcome["status"] = "pass"
        outcome["message"] = "every block's imbalance is within tolerance x initial norm"
        return outcome
    values, terms, scales, references, resolutions = [], [], [], [], []
    for label, _, n in grids:
        row = rows[label]
        paired = None
        if check.get("paired"):
            paired = rows.get(check["paired"] + label[label.index("="):])
            if paired is None:
                outcome["status"] = "missing_result"
                outcome["message"] = f"no paired {check['paired']} run for {label}"
                return outcome
        value = quantity_of(row, records[label], check["quantity"], check, paired)
        if value is None:
            outcome["status"] = "missing_result"
            outcome["message"] = f"{check['quantity']} absent from {label}"
            return outcome
        if not math.isfinite(value):
            outcome["status"] = "nonfinite"
            outcome["message"] = f"{check['quantity']} is {value} in {label}"
            return outcome
        values.append(value)
        terms.append(n)
        resolutions.append(resolution_of(check, row, value))
        if kind == "order":
            reference = check["reference"]
            if isinstance(reference, dict):
                reference = quantity_of(rows[reference["run"]], records[reference["run"]], check["quantity"])
                if reference is None:
                    outcome["status"] = "missing_result"
                    outcome["message"] = "the reference run lacks the quantity"
                    return outcome
            references.append(reference)
            scales.append(abs(reference) if reference is not None else 1.0)
    outcome["values"] = values
    if kind == "order":
        outcome["order"] = check["order"]
        outcome["references"] = references
        outcome.update(order_of(check, grids, values, references, terms, scales, resolutions))
        return outcome
    if kind == "value":
        outcome["status"] = "pass"
        outcome["message"] = "recorded: " + ", ".join(f"{v:.6g}" for v in values)
        return outcome
    if kind == "floor":
        outcome["thresholds"] = []
        outcome["errors"] = []
        for (label, _, n), value in zip(grids, values):
            threshold = floor_threshold(check, value, rows[label], n, case)
            error = abs(value - check["reference"])
            outcome["thresholds"].append(threshold)
            outcome["errors"].append(error)
            if error > threshold:
                outcome["status"] = "floor_exceeded"
                outcome["message"] = f"{label}: |value - reference| = {error:.3e} above the floor {threshold:.3e}"
                return outcome
        outcome["status"] = "pass"
        outcome["message"] = "within the roundoff floor on every grid"
        return outcome
    raise ValueError(kind)


def evaluate_case(case, executor):
    """Run a case's records and evaluate every check."""
    start = time.perf_counter()
    result = {"id": case["id"], "set": case["set"], "description": case["description"],
              "row": case["row"], "limitation": case.get("limitation"),
              "expected_status": case.get("expected_status"), "runs": [], "checks": [],
              "grids": [{"label": l, "h": h, "terms": n} for l, h, n in case["grids"]]}

    def close(status, message):
        result["status"] = status
        result["message"] = message
        result["elapsed_seconds"] = time.perf_counter() - start
        return result

    for check in case["checks"]:
        reference = check.get("reference")
        if isinstance(reference, dict):
            ratio = case["grids"][-2][1] / case["grids"][-1][1]
            bias = reference_bias(check["order"], ratio, reference["refinement"])
            if bias > window(ratio):
                return close("under_resolved_reference",
                             f"a reference {reference['refinement']} times finer than the finest grid "
                             f"shifts the slope by {bias:.3f}, beyond the window {window(ratio):.3f}")
    parser = PARSERS[case.get("record", "table")]
    records, rows = {}, {}
    for label, argv in case["runs"].items():
        outcome = executor.run(case["id"], label, argv, case["timeout"])
        result["runs"].append({k: v for k, v in outcome.items() if k != "text"})
        if outcome["timed_out"]:
            return close("timeout", f"{label} exceeded {case['timeout']} s")
        if outcome["exit_status"] != 0:
            return close("process_failure", f"{label} returned status {outcome['exit_status']}")
        try:
            record = parser(outcome["text"])
        except MalformedRecord as failure:
            return close("malformed_record", f"{label}: {failure}")
        records[label] = record
        if case["row"] is None:
            continue
        if case["row"] not in record["rows"]:
            return close("missing_result", f"row {case['row']} absent from {label}")
        row = record["rows"][case["row"]]
        if row.get("status"):
            return close("solver_failure", f"{label}: row {case['row']} is {row['status']}")
        rows[label] = row
    if case["row"] is None:
        for label in case["runs"]:
            rows[label] = None
    status, message = "pass", "every declaration met"
    for check in case["checks"]:
        if case["row"] is None and case.get("record") == "table":
            raise ValueError("a table case names its row")
        if case["row"] is None:
            for record_label, record in records.items():
                for row_label, row in record["rows"].items():
                    one = evaluate_check(check, case, {record_label: row}, records)
                    one["row"] = row_label
                    result["checks"].append(one)
                    if one["status"] != "pass" and status == "pass":
                        status, message = one["status"], f"{row_label}: {one.get('message', '')}"
            continue
        one = evaluate_check(check, case, rows, records)
        result["checks"].append(one)
        if one["status"] != "pass" and status == "pass":
            status, message = one["status"], f"{check['quantity']}: {one.get('message', '')}"
    if case.get("limitation"):
        if status == "pass":
            return close("unexpected_pass", "the declared limitation is met: update the declaration")
        if status in DECLARATION_FAILURES:
            return close("limitation", f"declared: {case['limitation']}; measured: {message}")
    return close(status, message)


# ----------------------------------------------------------------- provenance

def provenance(root):
    def git(*args):
        try:
            return subprocess.run(["git", "-C", str(root)] + list(args), capture_output=True,
                                  text=True, check=False).stdout.strip()
        except OSError:
            return ""

    def sha256(path):
        try:
            return hashlib.sha256(path.read_bytes()).hexdigest()
        except OSError:
            return None

    flags = re.search(r'^FLAGS="([^"]*)"', (root / "application/build.sh").read_text(), re.M)
    compiler = re.search(r"^F90=\$\{F90:-([^}]*)\}", (root / "application/build.sh").read_text(), re.M)
    compiler_name = os.environ.get("F90", compiler.group(1) if compiler else "gfortran")
    try:
        version = subprocess.run([compiler_name, "--version"], capture_output=True, text=True,
                                 check=False).stdout.splitlines()[0]
    except (OSError, IndexError):
        version = None
    return {"commit": git("rev-parse", "HEAD"), "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "working_tree_modified": bool(git("status", "--porcelain", "--untracked-files=no")),
            "compiler": compiler_name, "compiler_version": version,
            "flags": flags.group(1) if flags else None, "precision": "real64",
            "executable": declared.APPLICATION,
            "executable_sha256": sha256(root / declared.APPLICATION),
            "source_sha256": sha256(root / "application/module_graph_time_integrator.f90"),
            "python": platform.python_version(), "platform": platform.platform(),
            "hostname": platform.node(), "affinity": sorted(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity") else None}


# ----------------------------------------------------------------- main

def report_line(result):
    parts = []
    for check in result.get("checks", []):
        if check["kind"] == "order" and "observed_order" in check:
            parts.append(f"{check['quantity']} {check['observed_order']:.2f}/{check['order']}")
        elif check["kind"] == "order":
            parts.append(f"{check['quantity']} {check['status']}")
    return " ".join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--set", default="required", choices=["required", "exploratory", "rejection"])
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--results", default=None)
    parser.add_argument("--root", default=None)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = Path(args.root).resolve() if args.root else here.parents[1]
    results = Path(args.results).resolve() if args.results else here / "results" / args.set
    results.mkdir(parents=True, exist_ok=True)
    if args.set == "required":
        selected = declared.required_cases()
    elif args.set == "exploratory":
        selected = declared.exploratory_cases()
    else:
        selected = declared.rejection_cases(here / "fixtures")
    if args.case:
        selected = [c for c in selected if c["id"] in args.case]
        if not selected:
            print(f" no case named {args.case} in the {args.set} set")
            return 2
    if not (root / declared.APPLICATION).exists():
        print(f" FAIL : {declared.APPLICATION} is not built under {root}")
        return 2
    executor = Executor(root, results)
    start = time.perf_counter()
    summary = {"contract": "accuracy-contract", "set": args.set,
               "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "provenance": provenance(root),
               "constants": {"theta": THETA, "unit_roundoff": UNIT_ROUNDOFF,
                             "window_ratio_2": window(2.0), "tolerance": declared.TOLERANCE,
                             "references": {"E": declared.E_REF, "dE": declared.DE_REF,
                                            "dD": declared.DD_REF, "q": declared.Q_REF,
                                            "qd": declared.QD_REF}},
               "cases": []}
    print(f" ACCURACY CONTRACT ({args.set} set), window at ratio 2: {window(2.0):.4f}")
    counts = {}
    accepted_all = True
    for case in selected:
        result = evaluate_case(case, executor)
        summary["cases"].append(result)
        status = result["status"]
        counts[status] = counts.get(status, 0) + 1
        if args.set == "rejection":
            accepted = status == case["expected_status"]
            mark = "PASS" if accepted else "FAIL"
            print(f" {mark} : {case['id']:32s} rejected as {status:26s} expected {case['expected_status']}"
                  f"   {result['elapsed_seconds']:6.1f} s")
        else:
            accepted = status in ACCEPTED
            mark = "PASS" if status == "pass" else ("LIMIT" if status == "limitation" else "FAIL")
            print(f" {mark:5s}: {case['id']:32s} {status:14s} {report_line(result):44s} {result['elapsed_seconds']:6.1f} s")
            if status != "pass":
                print(f"          {result['message']}")
        if not accepted:
            accepted_all = False
    if args.set == "exploratory":
        accepted_all = all(c["status"] not in ("process_failure", "timeout", "malformed_record")
                 for c in summary["cases"])
    summary["counts"] = counts
    summary["elapsed_seconds"] = time.perf_counter() - start
    summary["status"] = "pass" if accepted_all else "fail"
    path = results / "summary.json"
    path.write_text(json.dumps(summary, indent=1, default=str))
    print(f" ACCURACY_CONTRACT set={args.set} status={summary['status']} "
          + " ".join(f"{k}={v}" for k, v in sorted(counts.items()))
          + f" elapsed={summary['elapsed_seconds']:.1f}s summary={path}")
    return 0 if accepted_all else 1


if __name__ == "__main__":
    sys.exit(main())
