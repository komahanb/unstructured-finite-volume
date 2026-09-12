#!/bin/bash
set -euo pipefail
suite_dir="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$suite_dir/../.." && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    (cd "$root" && ./build.sh)
fi
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
module_dir=${UFVM_GTI_MODULE_DIR:-$work}
compiler=${F90:-gfortran-15}
if [ ! -f "$module_dir/gti_demos.mod" ]; then
    sed '/^program graph_time_integrator/,$d' "$root/application/module_graph_time_integrator.f90" > "$module_dir/modules.f90"
    "$compiler" -std=f2023 -fcoarray=single -cpp -fbounds-check -O2 -I"$root/lib" -J"$module_dir" \
        -c "$module_dir/modules.f90" -o "$module_dir/modules.o"
fi
"$compiler" -std=f2023 -fcoarray=single -fbounds-check -O2 -I"$root/lib" -I"$module_dir" -J"$work" \
    "$suite_dir/test.f90" "$module_dir/modules.o" "$root/lib/libufvm.a" -o "$work/run"
"$work/run"
for mode in source_twin advance_uninitialized results_incomplete derivative_incomplete derivative_streamed; do
    case "$mode" in
        source_twin) diagnostic="view_level: this storage's hierarchy has been released" ;;
        advance_uninitialized) diagnostic='gti_chain: initialize an execution before advancing it' ;;
        results_incomplete) diagnostic='gti_chain: results require a completed execution' ;;
        derivative_incomplete) diagnostic='gti_chain: a derivative requires a completed primal execution' ;;
        derivative_streamed) diagnostic='gti_chain: a streamed Taylor execution has released its primal state' ;;
    esac
    if "$work/run" "$mode" > "$work/$mode.log" 2>&1; then
        echo "FAIL: execution accepted $mode"
        exit 1
    fi
    if ! grep -Fq "$diagnostic" "$work/$mode.log"; then
        cat "$work/$mode.log"
        exit 1
    fi
    echo "PASS: execution refused $mode"
done
