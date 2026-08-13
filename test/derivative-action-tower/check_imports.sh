#!/bin/bash
# The derivative action tower's import gate: the dependency-
# stratification law, made mechanical - the third client of the
# calculator gate's philosophy, with its own allowlists. Every
# derivative source may `use` only the framework modules its level
# has been explicitly granted; a directory with sources but no
# allowlist fails closed. Gate A especially forbids everywhere:
# graph_minimization, class_graph_gmres, any legacy tangent/adjoint
# or linearization machinery - structure first, numbers never, at
# this gate. This gate audits the DERIVATIVE TESTS' imports only.

here="$(cd "$(dirname "$0")" && pwd)"

intrinsics="iso_fortran_env iso_c_binding ieee_arithmetic ieee_exceptions ieee_features"

allowed_for() {
    case "$1" in
        common)                   echo "" ;;
        level-0-carrier)          echo "derivative_assert graph_carrier" ;;
        level-1-relation)         echo "derivative_assert graph_carrier graph_relation" ;;
        level-2-relation-algebra) echo "derivative_assert graph_carrier graph_relation graph_relation_algebra" ;;
        level-3-graph)            echo "derivative_assert graph_carrier graph_relation graph_relation_algebra graph_structure" ;;
        level-4-graph-calculus)   echo "derivative_assert graph_carrier graph_relation graph_relation_algebra graph_structure graph_profile graph_algorithms" ;;
        # level 5: primal values need domains, not graphs - and no
        # tangent/cotangent/seed type exists to be imported.
        level-5-field-calculus)   echo "derivative_assert graph_carrier class_graph_field" ;;
        # level 6: the structural rung - the algebra derives value
        # dependency, the profile walks it, and the binary citizen
        # materializes J_ZX and answers its reverse as a view.
        # Fields stay forbidden: the pattern needs no numbers.
        level-6-derivative-structure) echo "derivative_assert graph_carrier graph_relation graph_relation_algebra graph_binary_relation graph_structure graph_profile graph_algorithms" ;;
        *)                        echo "__no_allowlist__" ;;
    esac
}

violation=0

for dir in "$here"/common "$here"/level-*; do
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
                echo "IMPORT GATE: $(basename "$src") in $name uses '$mod' - not on the level's allowlist"
                violation=1
            fi
        done
    done
done

if [ "$violation" -ne 0 ]; then
    echo "IMPORT GATE: the tower layering is violated"
    exit 1
fi
echo "import gate: every derivative source imports only its level and below"
