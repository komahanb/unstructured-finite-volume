#!/bin/bash
# Independent concurrent executions: the heterogeneous set run serially and over
# an OpenMP parallel do, per-execution results compared at tolerance zero, each
# execution's own output file compared with its serial file, and the standard
# output of the runs at different thread counts compared line for line.
#
# OPENMP=yes builds the program with -fopenmp against a library built with the
# same setting (OPENMP=yes ./build.sh); the default build is the serial fallback,
# in which the requested thread count is ignored and every run is serial.
set -euo pipefail
suite_dir="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$suite_dir/../.." && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    (cd "$root" && ./build.sh)
fi
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
compiler=${F90:-gfortran-15}
flags="-std=f2023 -fcoarray=single -cpp -fbounds-check -O2"
openmp=${OPENMP:-no}
if [ "$openmp" = yes ]; then
    # -fopenmp also places every local array on the thread's stack
    # (-frecursive), so the application modules are compiled here with it
    flags="$flags -fopenmp"
    module_dir=$work
else
    module_dir=${UFVM_GTI_MODULE_DIR:-$work}
fi
if [ ! -f "$module_dir/gti_demos.mod" ]; then
    sed '/^program graph_time_integrator/,$d' "$root/application/module_graph_time_integrator.f90" > "$module_dir/modules.f90"
    "$compiler" $flags -I"$root/lib" -J"$module_dir" -c "$module_dir/modules.f90" -o "$module_dir/modules.o"
fi
"$compiler" $flags -I"$root/lib" -I"$module_dir" -J"$work" \
    "$suite_dir/test.f90" "$module_dir/modules.o" "$root/lib/libufvm.a" -o "$work/run"
if [ "$openmp" = yes ]; then
    counts="1 2 4"
else
    counts="1 4"
fi
for threads in $counts; do
    mkdir -p "$work/out_$threads"
    "$work/run" "$work/out_$threads" "$threads" | tee "$work/stdout_$threads.log"
done
first=$(echo $counts | cut -d' ' -f1)
for threads in $counts; do
    for f in "$work/out_$threads"/serial/execution_*.txt; do
        name=$(basename "$f")
        cmp "$f" "$work/out_$first/serial/$name"
        for d in "$work/out_$threads"/threads_*_repetition_*; do
            cmp "$f" "$d/$name"
        done
    done
    if ! diff <(grep -v '^ threads requested' "$work/stdout_$first.log") \
              <(grep -v '^ threads requested' "$work/stdout_$threads.log") > /dev/null; then
        echo "FAIL: standard output at $threads threads differs from $first"
        exit 1
    fi
done
echo "PASS: every execution's own file equals its serial file at every thread count ($counts)"
echo "PASS: standard output is independent of the thread count ($counts)"
