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

allowed_for() {
    case "$1" in
        common)            echo "" ;;
        # Gate A: structure, ownership, visibility, transport,
        # reconstruction. No operator, no solver.
        gate-a-partition)  echo "partitioned_pde_assert graph_carrier graph_grammar class_graph class_graph_field class_graph_partitioner class_graph_assembler" ;;
        *)                 echo "__no_allowlist__" ;;
    esac
}

violation=0

for dir in "$here"/common "$here"/gate-*; do
    [ -d "$dir" ] || continue
    name="$(basename "$dir")"
    sources=$(ls "$dir"/*.f90 2>/dev/null)
    [ -n "$sources" ] || continue

    allow="$(allowed_for "$name")"
    if [ "$allow" = "__no_allowlist__" ]; then
        echo "IMPORT GATE: $name has sources but no declared allowlist"
        violation=1
        continue
    fi

    for src in $sources; do
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
