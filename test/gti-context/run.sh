#!/bin/bash
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    (cd "$root" && ./build.sh)
    (cd "$root" && ./application/build.sh)
fi
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
module_dir=${UFVM_GTI_MODULE_DIR:-$work}
compiler=${F90:-gfortran-15}
if [ ! -f "$module_dir/gti_march.mod" ]; then
    # Standalone use needs only the prefix; verify.sh supplies the full module build.
    sed '/^end module gti_march/q' "$root/application/module_graph_time_integrator.f90" > "$module_dir/modules.f90"
    "$compiler" -std=f2023 -fcoarray=single -cpp -fbounds-check -O2 -I"$root/lib" -J"$module_dir" \
        -c "$module_dir/modules.f90" -o "$module_dir/modules.o"
fi
"$compiler" -std=f2023 -fcoarray=single -fbounds-check -O2 -I"$root/lib" -I"$module_dir" -J"$work" \
    "$here/test.f90" "$module_dir/modules.o" "$root/lib/libufvm.a" -o "$work/run"
"$work/run"
# The initial consistency solve uses the configured stopping rule too.
(cd "$root/application" && ./graph_time_integrator --families=bdf \
    --max_discretization_order=1 --max_derivative_degree=0 --instants=3 \
    --time_duration=0.1 --grid=uniform --tolerance=2.0 \
    --tolerance_criterion=absolute --max_iterations=1 --iteration_criterion=by_count) \
    > "$work/initial.log" 2>&1
if ! awk '/initial state, consistent/ {matched = ($4 == 1 && $5 == 0 && $6 == 0)} END {exit !matched}' \
    "$work/initial.log"; then
    cat "$work/initial.log"
    echo 'FAIL: the initial consistency solve ignored its configured absolute tolerance'
    exit 1
fi
echo 'PASS: the configured stopping rule applies before the initial consistency solve'
# The reverse storage limit (reverse_entries): a limit below one entry is
# refused where it is set; a limit below the working set of one
# recomputation is refused with the accounts before any tower is solved;
# a limit at retention runs the reverse pass, agreeing with the forward.
reverse_case() {
    (cd "$root/application" && ./graph_time_integrator --families=dirk --max_discretization_order=3 \
        --max_derivative_degree=1 --instants=13 --time_duration=1 --grid=uniform \
        --chain="dirk:3 dirk:3 dirk:3 dirk:3" --functionals="energy dissipation" --check=passes "$@")
}
if reverse_case --reverse_entries=0 > "$work/reverse-zero.log" 2>&1 || \
    ! grep -Fq 'gti_sweeps: the reverse storage limit must be one entry at least' "$work/reverse-zero.log"; then
    cat "$work/reverse-zero.log"
    echo 'FAIL: a reverse storage limit below one entry was accepted'
    exit 1
fi
if reverse_case --reverse_entries=1 > "$work/reverse-one.log" 2>&1 || \
    ! grep -Fq 'gti_chain: the reverse storage limit must admit the working set of one recomputation at least' "$work/reverse-one.log" || \
    ! grep -Eq 'reverse storage limit must admit the working set of one recomputation at least; limit 1 entries, recomputation requires [0-9]+ \(restart state [0-9]+, leaf [0-9]+, costate window [0-9]+, Lagrangian terms [0-9]+\); retention requires [0-9]+' "$work/reverse-one.log"; then
    cat "$work/reverse-one.log"
    echo 'FAIL: an insufficient reverse storage limit was accepted or reported without its accounts'
    exit 1
fi
# above the retention of every row (the rows differ in size; the smallest
# admissible limit of one chain is exercised by test/gti-execution): the
# limit is shown among the settings and every row's passes agree
if ! reverse_case --reverse_entries=100000 > "$work/reverse-retained.log" 2>&1 || \
    ! grep -Fq '   reverse entries          100000' "$work/reverse-retained.log" || \
    ! awk '/tangent against adjoint over the table, relative/ {rows++; if ($8 + 0.0 > 1.0e-12) failed++} \
        END {exit !(rows >= 2 && failed == 0)}' "$work/reverse-retained.log"; then
    cat "$work/reverse-retained.log"
    echo 'FAIL: the reverse pass under a limit above retention did not agree with the forward pass'
    exit 1
fi
echo 'PASS: the reverse storage limit is refused below one entry and below one recomputation, and applies above retention'
