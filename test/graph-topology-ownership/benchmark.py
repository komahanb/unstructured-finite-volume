"""Compare traversal against a preserved library from commit 8d59154."""

import argparse
import json
from pathlib import Path
import re
import statistics
import subprocess
import tempfile


def execute(command):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stdout + result.stderr)
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-library", type=Path, required=True)
    parser.add_argument("--candidate-library", type=Path, default=Path(__file__).resolve().parents[2] / "lib")
    parser.add_argument("--compiler", default="gfortran-15")
    parser.add_argument("--samples", type=int, default=7)
    parser.add_argument("--vertices", type=int, default=20000)
    parser.add_argument("--repetitions", type=int, default=1000)
    parser.add_argument("--cpu", type=int)
    arguments = parser.parse_args()
    directory = Path(__file__).resolve().parent
    measurements = {"baseline": {}, "candidate": {}}
    with tempfile.TemporaryDirectory(prefix="ufvm-topology-benchmark-") as temporary:
        temporary_path = Path(temporary)
        counter = temporary_path / "allocation_counter.o"
        execute(["cc", "-O2", "-c", str(directory / "allocation_counter.c"), "-o", str(counter)])
        for name, library in (("baseline", arguments.baseline_library), ("candidate", arguments.candidate_library)):
            command = [arguments.compiler, "-fcoarray=single", "-cpp", "-std=f2023", "-O3"]
            if name == "baseline":
                command.append("-DBASELINE")
            execute(command + ["-J" + str(temporary_path), "-I" + str(library.resolve()), str(directory / "traversal.F90"), str(counter),
                               str(library.resolve() / "libufvm.a"), "-Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc",
                               "-o", str(temporary_path / name)])
        for sample in range(arguments.samples):
            order = ("baseline", "candidate") if sample % 2 == 0 else ("candidate", "baseline")
            for name in order:
                command = [str(temporary_path / name), str(arguments.vertices), str(arguments.repetitions)]
                if arguments.cpu is not None:
                    command = ["taskset", "-c", str(arguments.cpu)] + command
                output = execute(command)
                for operation, seconds, allocations, checksum in re.findall(
                        r"(\w+) seconds=\s*([\d.]+) allocations=(\d+) checksum=(\d+)", output):
                    measurements[name].setdefault(operation, []).append({
                        "seconds": float(seconds), "allocations": int(allocations), "checksum": int(checksum)})
    report = {"vertices": arguments.vertices, "degree": 6, "repetitions": arguments.repetitions,
              "samples": arguments.samples, "cpu": arguments.cpu, "operations": {}}
    for operation in ("fibre_index", "fibre_read", "fibre_context", "incoming_read"):
        baseline = measurements["baseline"][operation]
        candidate = measurements["candidate"][operation]
        if len(baseline) != arguments.samples or len(candidate) != arguments.samples:
            raise RuntimeError("incomplete traversal measurements")
        if len({measurement["checksum"] for measurement in baseline + candidate}) != 1:
            raise RuntimeError("baseline and candidate incidence sums differ")
        baseline_median = statistics.median(measurement["seconds"] for measurement in baseline)
        candidate_median = statistics.median(measurement["seconds"] for measurement in candidate)
        report["operations"][operation] = {
            "baseline_median_seconds": baseline_median, "candidate_median_seconds": candidate_median,
            "candidate_over_baseline": candidate_median / baseline_median,
            "allocations": max(measurement["allocations"] for measurement in baseline + candidate),
            "checksum": baseline[0]["checksum"], "baseline": baseline, "candidate": candidate}
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
