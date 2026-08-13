#!/bin/bash
# The adjoint tower's import gate: the dependency-stratification
# law, made mechanical - the fourth client of the calculator gate's
# philosophy, with its own allowlists. Every adjoint source may
# `use` only the framework modules its LEVEL has been granted; a
# directory with sources but no allowlist fails closed.
#
# Gate grouping does not weaken level checking. In particular Gate
# A (levels 0-6) may not touch minimization or a solver at all, and
# levels 0-4 may not touch the field: structure precedes numbers,
# and no gate boundary excuses a rung from its own stratum.

here="$(cd "$(dirname "$0")" && pwd)"

intrinsics="iso_fortran_env iso_c_binding ieee_arithmetic ieee_exceptions ieee_features"

allowed_for() {
    case "$1" in
        common)                   echo "" ;;
        level-0-carrier)          echo "adjoint_assert graph_carrier" ;;
        level-1-relation)         echo "adjoint_assert graph_carrier graph_relation" ;;
        # level 2: the algebra, and the binary citizen for the
        # subobjects' own inclusions and their transposed use.
        level-2-relation-algebra) echo "adjoint_assert graph_carrier graph_relation graph_relation_algebra graph_binary_relation" ;;
        level-3-graph)            echo "adjoint_assert graph_carrier graph_relation graph_relation_algebra graph_binary_relation graph_structure" ;;
        # level 4: the profile and the algorithms - which are here to
        # REFUSE an order, not to produce one.
        level-4-graph-calculus)   echo "adjoint_assert graph_carrier graph_relation graph_relation_algebra graph_binary_relation graph_structure graph_profile graph_algorithms" ;;
        # level 5: values need domains, not graphs - the smallest
        # allowlist above the ground.
        level-5-field-calculus)   echo "adjoint_assert graph_carrier class_graph_field" ;;
        # level 6: support and orientation; still no field, because
        # where an operator stands is not what it multiplies.
        level-6-discretization)   echo "adjoint_assert graph_carrier graph_relation graph_relation_algebra graph_binary_relation" ;;
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
echo "import gate: every adjoint source imports only its level and below"
