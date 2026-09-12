#!/bin/bash
# build the library and the traversal benchmark, run it; then build the
# scaling programs and verify their oracles at small sizes.
set -e
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$root" && ./build.sh >/dev/null )
fi
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run

# The section-66 reference: printed beside every run so a regression
# reads as a before/after table.
echo ""
echo " --- baseline (test/graph-benchmark/baseline) ---"
cat "$here/baseline"

# The scaling programs: schedule construction, listed subsets, elimination
# fill and one GTI horizon. The horizon program links the packed
# application modules, compiled once into UFVM_GTI_MODULE_DIR when set.
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
module_dir=${UFVM_GTI_MODULE_DIR:-$work}
compiler=${F90:-gfortran-15}
if [ ! -f "$module_dir/gti_chain.mod" ]; then
    sed '/^program graph_time_integrator/,$d' "$root/application/module_graph_time_integrator.f90" > "$module_dir/modules.f90"
    "$compiler" -std=f2023 -fcoarray=single -cpp -fbounds-check -O2 -I"$root/lib" -J"$module_dir" \
        -c "$module_dir/modules.f90" -o "$module_dir/modules.o"
fi
make -C "$here" horizon GTI_MODULE_DIR="$module_dir" >/dev/null
echo ""
echo " --- scaling oracles at small sizes (python3 scaling.py --quick) ---"
python3 "$here/scaling.py" --quick --output-dir "$work/scaling"
