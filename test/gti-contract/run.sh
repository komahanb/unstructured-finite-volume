#!/bin/bash
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    (cd "$root" && ./build.sh)
fi
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
# Compile the application's modules independently of its main program.
sed '/^program graph_time_integrator/,$d' "$root/application/module_graph_time_integrator.f90" > "$work/modules.f90"
compiler=${F90:-gfortran-15}
"$compiler" -std=f2023 -fcoarray=single -cpp -fbounds-check -O2 -I"$root/lib" -J"$work" \
    -c "$work/modules.f90" -o "$work/modules.o"
"$compiler" -std=f2023 -fcoarray=single -fbounds-check -O2 -I"$root/lib" -I"$work" -J"$work" \
    "$here/test.f90" "$work/modules.o" "$root/lib/libufvm.a" -o "$work/run"
"$work/run" accuracy
"$work/run" status
for mode in adaptive_failure minimum_step forward reverse linear_forward linear_reverse; do
    case "$mode" in
        adaptive_failure) expected='step solve did not converge' ;;
        minimum_step) expected='minimum step cannot meet the tolerance' ;;
        forward|reverse) expected='primal march must converge' ;;
        linear_forward|linear_reverse) expected='derivative solve did not converge' ;;
    esac
    if "$work/run" "$mode" > "$work/$mode.log" 2>&1; then
        echo "FAIL: $mode was accepted"
        exit 1
    fi
    if ! grep -q "$expected" "$work/$mode.log"; then
        cat "$work/$mode.log"
        exit 1
    fi
    echo "PASS: $mode reports its failure"
done
