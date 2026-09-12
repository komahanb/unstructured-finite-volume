#!/usr/bin/env python3
"""Scaling measurements: schedule construction, listed subsets, elimination fill, a GTI horizon.

Every case runs as its own process, pinned when --cpus is given, repeated
--samples times; with several --build label=directory pairs the builds are
interleaved (A B, B A, ...) so machine drift affects each equally. Records
are the key=value lines the programs print; the peak resident size and
user time come from the kernel's rusage of each child. Output: records.json
(raw) and tables.md (medians, ranges, local log-log slopes).
"""

import argparse
import datetime
import json
import math
import os
import platform
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

TOKEN = re.compile(r"(\w+)=\s*(\S+)")


def cases(quick, suites):
    """The measured cases: (suite, program, arguments, series) tuples."""
    chosen = []
    geometric = lambda first, count, ratio=2: [first * ratio ** k for k in range(count)]
    if "schedule" in suites:
        sizes = [200, 400] if quick else geometric(500, 7)
        for n in sizes:
            chosen.append(("schedule", "schedule", ["chain", str(n), "1"], "chain n"))
            chosen.append(("schedule", "schedule", ["tree", str(n), "1"], "tree n"))
        for reach in ([1, 4] if quick else [1, 2, 4, 8, 16, 32]):
            chosen.append(("schedule", "schedule", ["history", "400" if quick else "4000", str(reach)], "history reach"))
        for n in ([200, 400] if quick else geometric(500, 6)):
            chosen.append(("schedule", "schedule", ["history", str(n), "4"], "history n"))
    if "subset" in suites:
        sizes = [1000, 4000] if quick else geometric(1000, 6, 4)
        for m in sizes:
            reps = max(1, 400000 // m)
            chosen.append(("subset", "subset", ["construct", str(4 * m), str(m), "random", "0", "1", str(reps)], "construct M"))
        m = 4000 if quick else 64000
        for order in ["ascending", "descending", "random"]:
            for duplicated in ["0", "1"]:
                chosen.append(("subset", "subset", ["construct", str(4 * m), str(m), order, duplicated, "1", str(max(1, 400000 // m))], "construct order"))
        for n in ([4000, 16000] if quick else geometric(10000, 4, 4)):
            chosen.append(("subset", "subset", ["transport", str(n), str(n // 4), "random", "0", "4", "1"], "transport N"))
        n = 16000 if quick else 160000
        for parts in ([1, 4] if quick else [1, 4, 16, 64]):
            chosen.append(("subset", "subset", ["transport", str(n), str(n // 4), "random", "0", str(parts), "1"], "transport parts"))
        for divisor in ([4, 1] if quick else [64, 16, 4, 1]):
            chosen.append(("subset", "subset", ["transport", str(n), str(n // divisor), "random", "0", "4", "1"], "transport M"))
    if "schur" in suites:
        sizes = [40, 80] if quick else geometric(40, 7)
        for ne in sizes:
            chosen.append(("schur", "schur", [str(ne), "8", "2", "2", "2", "1e-9", "400"], "ne"))
        ne = "80" if quick else "640"
        for p in ([1, 2] if quick else [1, 2, 4, 8, 16]):
            chosen.append(("schur", "schur", [ne, "8", str(p), "2", "2", "1e-9", "400"], "p"))
        for nk in ([4, 16] if quick else [4, 16, 64, 256]):
            chosen.append(("schur", "schur", [ne, str(nk), "2", "2", "2", "1e-9", "400"], "nk"))
        for c in ([1, 4] if quick else [1, 2, 4, 8]):
            chosen.append(("schur", "schur", [ne, "64", "2", str(c), "2", "1e-9", "400"], "c"))
    if "horizon" in suites:
        for n in ([10, 20] if quick else geometric(20, 5)):
            chosen.append(("horizon", "horizon", [str(n), "primal", "4", "1e-12"], "primal n"))
        for n in ([10, 20] if quick else geometric(20, 4)):
            chosen.append(("horizon", "horizon", [str(n), "taylor", "4", "1e-12"], "taylor n"))
    if "application" in suites and not quick:
        for cells in ["4 4", "8 8", "16 16"]:
            for rows in ["states", "states state-time-derivatives state-spatial-derivatives"]:
                chosen.append(("application", "graph_time_integrator",
                               ["--config=taylor_green", "spatial_counts=" + cells, "instants=3", "time_duration=0.25",
                                "max_derivative_degree=0", "export=none", "rows=" + rows], "mesh " + rows.split()[0]
                               + ("_all" if " " in rows else "")))
    return chosen


def execute(command, cwd, cpus, priority):
    """Run one process; its rusage (peak resident size, user time) is read at wait4."""
    prefix = []
    if cpus:
        prefix += ["taskset", "-c", cpus]
    if priority is not None:
        prefix += ["nice", "-n", str(priority)]
    with tempfile.TemporaryFile(mode="w+") as out, tempfile.TemporaryFile(mode="w+") as err:
        started = time.monotonic()
        process = subprocess.Popen(prefix + command, cwd=cwd, stdout=out, stderr=err)
        _, status, usage = os.wait4(process.pid, 0)
        elapsed = time.monotonic() - started
        process.returncode = os.waitstatus_to_exitcode(status)
        out.seek(0)
        err.seek(0)
        stdout, stderr = out.read(), err.read()
    # child_maxrss includes the pre-exec image of this interpreter, a floor of
    # about 15 MB; the programs report their own peak from getrusage instead.
    return {"returncode": process.returncode, "seconds": elapsed, "user_seconds": usage.ru_utime,
            "system_seconds": usage.ru_stime, "child_maxrss_kilobytes": usage.ru_maxrss, "stdout": stdout, "stderr": stderr}


def parse(stdout):
    phases, verifications, summary = {}, {}, {}
    for line in stdout.splitlines():
        words = line.split()
        if not words:
            continue
        tokens = dict(TOKEN.findall(line))
        if words[0] == "record":
            phases[tokens["phase"]] = {"seconds": float(tokens["seconds"]), "calls": int(tokens["calls"]), "bytes": int(tokens["bytes"])}
        elif words[0] == "verification":
            verifications[tokens["name"]] = tokens["satisfied"] == "T"
        elif words[0] in ("summary", "numbers"):
            for key, value in tokens.items():
                if key in ("suite", "shape", "mode", "order"):
                    continue
                try:
                    summary[key] = int(value)
                except ValueError:
                    try:
                        summary[key] = float(value)
                    except ValueError:
                        summary[key] = value
    return phases, verifications, summary


def application_numbers(stdout):
    """Printed functional and exact-solution diagnostics of an application run, for equality across variants."""
    numbers = re.findall(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?", stdout)
    return numbers


def provenance(builds, cpus, priority):
    def output(command, cwd=ROOT):
        try:
            return subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=True).stdout.strip()
        except Exception as error:  # noqa: BLE001 - provenance is best effort
            return "unavailable: %s" % error
    return {
        "date": datetime.datetime.now().isoformat(timespec="seconds"),
        "commit": output(["git", "rev-parse", "HEAD"]),
        "branch": output(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "working_tree": output(["git", "status", "--short"]),
        "compiler": output(["gfortran-15", "--version"]).splitlines()[0],
        "library_flags": "build.sh: make OPTIMIZE=yes with DEBUG=yes: -O3 -g -Wall -pedantic -fbounds-check -fbacktrace "
                         "-ffpe-trap=invalid,overflow,underflow -std=f2023 -fcoarray=single -cpp -fPIC",
        "benchmark_flags": "test/graph-benchmark/Makefile: the same DEBUG flags without -O3, plus -O2",
        "application_flags": "application/build.sh: -std=f2023 -fcoarray=single -cpp -Wall -fbounds-check -O2",
        "precision": "real64 (util_precision dp)",
        "hardware": output(["sh", "-c", "lscpu | grep -E 'Model name|^CPU\\(s\\)|L2|L3' | sed 's/  */ /g'"]),
        "memory": output(["sh", "-c", "free -m | head -2"]),
        "kernel": platform.release(),
        "affinity": cpus or "none",
        "scheduling_priority": priority,
        "builds": {label: str(path) for label, path in builds.items()},
        "concurrent_activity": "other agents were building and testing on the same 4-core machine during these measurements",
    }


def measure(arguments):
    builds = {}
    for pair in arguments.build:
        label, _, path = pair.partition("=")
        builds[label] = Path(path).resolve()
    suites = arguments.suites.split(",")
    chosen = cases(arguments.quick, suites)
    output_dir = Path(arguments.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = {"provenance": provenance(builds, arguments.cpus, arguments.priority), "samples": arguments.samples,
               "quick": arguments.quick, "cases": []}
    failures = 0
    for suite, program, args, series in chosen:
        entry = {"suite": suite, "program": program, "arguments": args, "series": series, "runs": {label: [] for label in builds}}
        records["cases"].append(entry)
        for sample in range(arguments.samples):
            labels = list(builds)
            if sample % 2 == 1:
                labels.reverse()
            for label in labels:
                executable = builds[label] / program
                cwd = ROOT / "application" if suite == "application" else builds[label]
                run = execute([str(executable)] + args, cwd, arguments.cpus, arguments.priority)
                phases, verifications, summary = parse(run["stdout"])
                run.update({"phases": phases, "verifications": verifications, "summary": summary})
                if suite == "application":
                    run["printed_numbers"] = application_numbers(run["stdout"])
                    run["phases"] = {"process": {"seconds": run["seconds"], "calls": -1, "bytes": -1}}
                if run["returncode"] != 0 or not all(verifications.values()):
                    failures += 1
                    (output_dir / ("failure-%s-%s-%d.log" % (suite, "_".join(a.replace(" ", "") for a in args), sample))).write_text(run["stdout"] + run["stderr"])
                run.pop("stdout")
                run.pop("stderr")
                entry["runs"][label].append(run)
                print("%s %s %s sample %d: %s seconds=%.4f maxrss=%d KB %s" % (
                    label, suite, " ".join(args), sample, "ok" if run["returncode"] == 0 else "FAILED", run["seconds"],
                    run["summary"].get("peak_rss_kilobytes", run["child_maxrss_kilobytes"]), " ".join("%s=%.3es" % (k, v["seconds"]) for k, v in run["phases"].items())), flush=True)
    records["failures"] = failures
    (output_dir / "records.json").write_text(json.dumps(records, indent=1))
    (output_dir / "tables.md").write_text(tables(records))
    print("failures: %d; records in %s" % (failures, output_dir))
    return failures


def median_range(values):
    if not values:
        return float("nan"), float("nan"), float("nan")
    return statistics.median(values), min(values), max(values)


def scientific(value):
    if isinstance(value, str):
        return value
    if value != value:
        return "nan"
    if abs(value) >= 1000 or (abs(value) < 0.01 and value != 0):
        return "%.3e" % value
    return "%.4g" % value


def tables(records):
    """Markdown tables per suite and series: one row per case, phase medians, ranges and slopes."""
    lines = ["# Scaling tables", "", "Commit %s; %d samples; affinity %s; %s" % (
        records["provenance"]["commit"], records["samples"], records["provenance"]["affinity"],
        records["provenance"]["concurrent_activity"]), ""]
    builds = list(records["cases"][0]["runs"]) if records["cases"] else []
    grouped = {}
    for case in records["cases"]:
        grouped.setdefault((case["suite"], case["series"]), []).append(case)
    for (suite, series), group in grouped.items():
        phases = []
        for case in group:
            for label in builds:
                for run in case["runs"][label]:
                    for name in run["phases"]:
                        if name not in phases:
                            phases.append(name)
        lines.append("## %s: %s" % (suite, series))
        lines.append("")
        header = ["arguments", "build"] + ["%s median s (min..max)" % p for p in phases] + ["slope"] + ["calls " + p for p in phases] + ["peak RSS KB"]
        lines.append("| " + " | ".join(header) + " |")
        lines.append("|" + "---|" * len(header))
        previous = {}
        for case in group:
            for label in builds:
                runs = [r for r in case["runs"][label] if r["returncode"] == 0]
                cells = [" ".join(case["arguments"]), label]
                medians = {}
                for p in phases:
                    values = [r["phases"][p]["seconds"] for r in runs if p in r["phases"]]
                    m, low, high = median_range(values)
                    medians[p] = m
                    cells.append("%s (%s..%s)" % (scientific(m), scientific(low), scientific(high)) if values else "-")
                size = case_size(case)
                slope = "-"
                if label in previous and previous[label][0] and size and phases:
                    p = phases[1] if suite in ("schedule",) and len(phases) > 1 else phases[0]
                    p = {"schedule": "schedule", "subset": phases[0], "schur": "state", "horizon": "initialize"}.get(suite, p)
                    if p in medians and p in previous[label][1] and previous[label][1][p] > 0 and medians[p] > 0 and size != previous[label][0]:
                        slope = "%s %.2f" % (p, math.log(medians[p] / previous[label][1][p]) / math.log(size / previous[label][0]))
                previous[label] = (size, medians)
                cells.append(slope)
                for p in phases:
                    values = [r["phases"][p]["calls"] for r in runs if p in r["phases"]]
                    cells.append(str(int(statistics.median(values))) if values else "-")
                values = [r["summary"].get("peak_rss_kilobytes", r["child_maxrss_kilobytes"]) for r in runs]
                cells.append(str(int(statistics.median(values))) if values else "-")
                lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
        summaries = {}
        for case in group:
            for label in builds:
                for run in case["runs"][label]:
                    if run["summary"]:
                        summaries[(" ".join(case["arguments"]), label)] = run["summary"]
        if summaries:
            keys = []
            for summary in summaries.values():
                for key in summary:
                    if key not in keys and key != "peak_rss_kilobytes":
                        keys.append(key)
            lines.append("| arguments | build | " + " | ".join(keys) + " |")
            lines.append("|---|---|" + "---|" * len(keys))
            for (args, label), summary in summaries.items():
                lines.append("| %s | %s | %s |" % (args, label, " | ".join(scientific(summary.get(k, "-")) for k in keys)))
            lines.append("")
        if len(builds) > 1:
            lines.append("Ratios of medians, %s over %s (slowdowns above 1.10 marked):" % (builds[1], builds[0]))
            lines.append("")
            for case in group:
                base = case["runs"][builds[0]]
                cand = case["runs"][builds[1]]
                parts = []
                for p in phases:
                    b = [r["phases"][p]["seconds"] for r in base if p in r["phases"] and r["returncode"] == 0]
                    c = [r["phases"][p]["seconds"] for r in cand if p in r["phases"] and r["returncode"] == 0]
                    if b and c and statistics.median(b) > 0:
                        ratio = statistics.median(c) / statistics.median(b)
                        parts.append("%s %.3f%s" % (p, ratio, " SLOWER" if ratio > 1.10 else ""))
                lines.append("- %s: %s" % (" ".join(case["arguments"]), "; ".join(parts)))
            lines.append("")
    return "\n".join(lines)


def case_size(case):
    a = case["arguments"]
    suite = case["suite"]
    series = case["series"]
    try:
        if suite == "schedule":
            return int(a[2]) if series.endswith("reach") else int(a[1])
        if suite == "subset":
            return {"construct M": int(a[2]), "transport N": int(a[1]), "transport parts": int(a[5]), "transport M": int(a[2])}.get(series)
        if suite == "schur":
            return {"ne": int(a[0]), "p": int(a[2]), "nk": int(a[1]), "c": int(a[3])}.get(series)
        if suite == "horizon":
            return int(a[0])
        if suite == "application":
            return int(a[1].split("=")[1].split()[0])
    except (ValueError, IndexError):
        return None
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", default=str(HERE / "scaling-records"))
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--cpus", default=None, help="processor list for taskset, e.g. 2,3")
    parser.add_argument("--priority", type=int, default=None, help="nice increment for every measured process")
    parser.add_argument("--build", action="append", default=[], help="label=directory of built benchmark programs; repeatable")
    parser.add_argument("--suites", default="schedule,subset,schur,horizon", help="comma-separated subset of schedule,subset,schur,horizon,application")
    parser.add_argument("--quick", action="store_true", help="small sizes, one sample: the oracle check for run.sh")
    parser.add_argument("--tables", default=None, help="rewrite tables.md from an existing records.json and exit")
    arguments = parser.parse_args()
    if arguments.tables:
        records = json.loads(Path(arguments.tables).read_text())
        print(tables(records))
        return 0
    if not arguments.build:
        arguments.build = ["current=" + str(HERE)]
    if arguments.quick:
        arguments.samples = 1
    return 1 if measure(arguments) else 0


if __name__ == "__main__":
    sys.exit(main())
