#!/bin/bash
# The partitioned PDE tower's import gate: the dependency-
# stratification law, made mechanical - the fifth client of the
# calculator gate's philosophy, with its own allowlists. Every
# source may `use` only the framework modules its GATE has been
# granted; a directory with sources but no allowlist fails closed.
#
# The gate grouping does not weaken the law. In particular Gate A
# may not touch a differential operator or a solver at all:
# structure, ownership and transport must stand before any operator
# consumes them.
#
# This is deliberately a NON-DERIVATIVE-FAMILY client: the
# derivative-action and adjoint fixtures, and class_graph_linearization,
# are forbidden everywhere.

here="$(cd "$(dirname "$0")" && pwd)"

intrinsics="iso_fortran_env iso_c_binding ieee_arithmetic ieee_exceptions ieee_features"

# Allowlists are keyed by "<dir>/<file>" first, then by "<dir>".
# The per-file form exists so that common/ can hold BOTH the
# dependency-free assert module and a fixture that legitimately
# imports production - without the fixture's permissions silently
# extending to the assert module.
allowed_for() {
    case "$1" in
        # The assert module imports nothing but intrinsics, and the
        # gate proves it rather than trusting the comment.
        common/partitioned_pde_assert.f90) echo "" ;;
        # The shifted adapter delegates topology to production; it
        # may see the differential operator and the field, and it
        # must NOT see the partitioner, the assembler or a solver.
        common/shifted_laplacian_fixture.f90) echo "graph_carrier graph_grammar class_graph_field class_graph_differential_operator" ;;
        common)            echo "__no_allowlist__" ;;
        # Gate A: structure, ownership, visibility, transport,
        # reconstruction. No operator, no solver.
        gate-a-partition)  echo "partitioned_pde_assert graph_carrier graph_grammar class_graph class_graph_field class_graph_partitioner class_graph_assembler" ;;
        # Gate B: the topology-consuming action. Adds the production
        # differential operator and the GMRES citizen, plus the
        # tower's own shifted adapter. Still no derivative-family
        # fixture and no class_graph_linearization anywhere.
        gate-b-operator)   echo "partitioned_pde_assert shifted_laplacian_fixture graph_carrier graph_grammar class_graph class_graph_field class_graph_partitioner class_graph_assembler class_graph_differential_operator class_graph_gmres" ;;
        # The composite is built from ONE decomposition: it needs the
        # partition/assembly stack and the Gate-B local action, and
        # nothing from any other tower.
        common/partitioned_shifted_laplacian_fixture.f90) echo "graph_carrier graph_grammar class_graph class_graph_field class_graph_partitioner class_graph_assembler shifted_laplacian_fixture" ;;
        # Gate C: the statement. Composition only - the two fixtures,
        # the transport stack and the solver. No linearization, no
        # derivative or adjoint fixture, no relational-graph machinery.
        gate-c-statement)  echo "partitioned_pde_assert shifted_laplacian_fixture partitioned_shifted_laplacian_fixture graph_carrier graph_grammar class_graph class_graph_field class_graph_partitioner class_graph_assembler class_graph_gmres" ;;
        *)                 echo "__no_allowlist__" ;;
    esac
}

violation=0

for dir in "$here"/common "$here"/gate-*; do
    [ -d "$dir" ] || continue
    name="$(basename "$dir")"
    sources=$(ls "$dir"/*.f90 2>/dev/null)
    [ -n "$sources" ] || continue

    for src in $sources; do
        # per-file allowlist first, then the directory's
        allow="$(allowed_for "$name/$(basename "$src")")"
        if [ "$allow" = "__no_allowlist__" ]; then
            allow="$(allowed_for "$name")"
        fi
        if [ "$allow" = "__no_allowlist__" ]; then
            echo "IMPORT GATE: $(basename "$src") in $name has no declared allowlist"
            violation=1
            continue
        fi
        uses=$(grep -ihE '^[[:space:]]*use[[:space:]]' "$src" \
               | sed -E 's/^[[:space:]]*[uU][sS][eE][[:space:]]*(,[[:space:]]*[iI][nN][tT][rR][iI][nN][sS][iI][cC][[:space:]]*::)?[[:space:]]*([a-zA-Z][a-zA-Z0-9_]*).*/\2/' \
               | tr 'A-Z' 'a-z' | sort -u)
        for mod in $uses; do
            ok=0
            for a in $allow $intrinsics; do
                [ "$mod" = "$a" ] && ok=1 && break
            done
            if [ "$ok" -eq 0 ]; then
                echo "IMPORT GATE: $(basename "$src") in $name uses '$mod' - not on the gate's allowlist"
                violation=1
            fi
        done
    done
done

if [ "$violation" -ne 0 ]; then
    echo "IMPORT GATE: the tower layering is violated"
    exit 1
fi
echo "import gate: every source imports only its gate and below"
