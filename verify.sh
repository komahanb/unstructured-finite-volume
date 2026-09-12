#!/bin/bash
# Build once, then verify solver, domain, scheduling and derivative contracts.
set -euo pipefail
cd "$(dirname "$0")"
./check_naming.sh
python3 test/naming/test.py
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ./build.sh
    ./application/build.sh
fi
export UFVM_SKIP_LIBRARY_BUILD=1
gti_module_dir=$(mktemp -d)
trap 'rm -rf "$gti_module_dir"' EXIT
export UFVM_GTI_MODULE_DIR="$gti_module_dir"
for suite in graph-contract graph-robustness graph-minimization graph-multigrid \
    graph-dense-direct graph-partition graph-field-transport graph-set-view \
    graph-marching graph-execution graph-elimination gti-contract gti-context gti-execution graph-benchmark time-integration-tower \
    adjoint-tower derivative-action-tower fractal-graph fractal-map \
    graph-algebra graph-algorithms graph-binary graph-topology-ownership graph-change graph-characterization \
    graph-constitution graph-differentiation graph-field graph-identity-map \
    graph-inclusion graph-mesh graph-ordinary graph-relation graph-relational \
    graph-sequence graph-state calculator-tower learning-tower visualization-tower \
    partitioned-implicit-pde-tower; do
    echo "Verifying $suite"
    "./test/$suite/run.sh"
done
for demo in function_identities scheme_weights randomized_checks lagrangian_expansion taylor_state; do
    (cd application && ./graph_time_integrator "--demo=$demo")
done
echo 'All numerical and contract checks passed.'
