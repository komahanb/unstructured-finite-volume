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
