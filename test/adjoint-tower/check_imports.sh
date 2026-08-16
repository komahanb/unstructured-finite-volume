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
        # 2026-08-16: the relational container is retired. A level reading
        # (S, P) is granted fractal_graph and graph_relational_view,
        # and builds the representation itself. Granted per level, in
        # review; the list is an assertion, not a history.
        common)                   echo "" ;;
        level-0-carrier)          echo "adjoint_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map" ;;
        level-1-relation)         echo "adjoint_assert fractal_graph graph_set_representation graph_set_map graph_relation" ;;
        # level 2: the algebra, and the binary citizen for the
        # subobjects' own inclusions and their transposed use.
        level-2-relation-algebra) echo "adjoint_assert fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation" ;;
        level-3-graph)            echo "adjoint_assert fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation fractal_graph graph_relational_view" ;;
        # level 4: the profile and the algorithms - which are here to
        # REFUSE an order, not to produce one.
        level-4-graph-calculus)   echo "adjoint_assert fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation graph_profile graph_algorithms fractal_graph graph_relational_view" ;;
        # level 5: values need domains, not graphs - the smallest
        # allowlist above the ground.
        level-5-field-calculus)   echo "adjoint_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map class_graph_field" ;;
        # level 6: support and orientation; still no field, because
        # where an operator stands is not what it multiplies.
        level-6-discretization)   echo "adjoint_assert fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation" ;;
        # level 7: the solver rung. gmres inherits attach/constant/
        # apply from the minimizer base, so graph_minimization is not
        # imported directly; the equations are SUPPLIED by the level's
        # own fixture, and no Level-8 constitution may be reached.
        level-7-minimization)     echo "adjoint_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_operation_view graph_ordinary_view graph_field_calculus class_graph_field class_graph class_graph_gmres opaque_equation_fixture" ;;
        # level 8: one constitution. It may see everything legitimately
        # below it - the relation stack that carries the structural
        # supports, the field, and the solver - plus its own fixture.
        # It may NOT reach back into Level 7's supplied equations.
        level-8-constitution)     echo "adjoint_assert fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation graph_operation_view graph_ordinary_view graph_field_calculus class_graph_field class_graph class_graph_gmres adjoint_constitution_fixture" ;;
        # level 9: the statement - the composition rung. It may see
        # the relation stack that carries the supports, the model
        # graph that owns them, the field, the legacy host and the
        # solver, plus the REUSED level-8 constitution. It may NOT
        # import class_graph_linearization (the specialized
        # same-domain path this tower cannot use), nor the graph
        # algorithms: an implicit system does not become a DAG at
        # the statement.
        level-9-statement)        echo "adjoint_assert fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_gmres adjoint_constitution_fixture fractal_graph graph_relational_view" ;;
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
