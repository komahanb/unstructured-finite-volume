#!/bin/bash
# The partitioned PDE tower's import gate, keyed PER LEVEL.
#
# Levels are the implementation architecture, so the dependency
# ceiling rises level by level - never gate by gate. A level source
# may `use` only the framework modules its own rung has been
# granted, and a directory with sources but no allowlist fails
# closed.
#
# common/ is NOT a hole in the stratification. Each shared fixture
# is keyed by FILE and classified by the first level that earns it:
#
#     partitioned_pde_assert            below everything (nothing)
#     chain_relations_fixture           earned at Level 1
#     chain_algebra_fixture             earned at Level 2
#     shifted_laplacian_fixture         earned at Level 6
#     partitioned_shifted_laplacian     earned at Level 8
#
# so Level 5 cannot reach the Level-6 operator, and Level 7 cannot
# reach the Level-8 composite. The per-file keying also proves the
# assert module's freedom from framework imports mechanically
# rather than by comment.
#
# This is deliberately a NON-DERIVATIVE-FAMILY client: the
# derivative and adjoint fixtures, and class_graph_linearization,
# are forbidden at every level.

here="$(cd "$(dirname "$0")" && pwd)"

intrinsics="iso_fortran_env iso_c_binding ieee_arithmetic ieee_exceptions ieee_features"

allowed_for() {
    case "$1" in
        # ---- shared fixtures, keyed by the level that earns them
        common/partitioned_pde_assert.f90) echo "" ;;
        common/chain_relations_fixture.f90) echo "graph_carrier graph_relation graph_binary_relation" ;;
        common/chain_algebra_fixture.f90) echo "graph_carrier graph_relation graph_relation_algebra graph_binary_relation" ;;
        common/shifted_laplacian_fixture.f90) echo "graph_carrier graph_grammar class_graph_field class_graph_differential_operator" ;;
        common/partitioned_shifted_laplacian_fixture.f90) echo "graph_carrier graph_grammar class_graph class_graph_field class_graph_partitioner class_graph_assembler shifted_laplacian_fixture" ;;
        common)            echo "__no_allowlist__" ;;

        # ---- L0: carriers only. No relation, no graph, no field.
        level-0-carrier)   echo "partitioned_pde_assert chain_relations_fixture graph_carrier" ;;
        # ---- L1: + the relation nucleus
        level-1-relation)  echo "partitioned_pde_assert chain_relations_fixture graph_carrier graph_relation graph_binary_relation" ;;
        # ---- L2: + relation algebra
        level-2-relation-algebra) echo "partitioned_pde_assert chain_relations_fixture chain_algebra_fixture graph_carrier graph_relation graph_relation_algebra graph_binary_relation" ;;
        # ---- L3: + the ordinary graph realization. No partitioner.
        level-3-graph)     echo "partitioned_pde_assert chain_relations_fixture graph_carrier graph_relation graph_binary_relation class_graph" ;;
        # ---- L4: + the partitioner. Graph-to-graph only; no field.
        level-4-graph-calculus) echo "partitioned_pde_assert chain_relations_fixture chain_algebra_fixture graph_carrier graph_grammar graph_relation graph_relation_algebra graph_binary_relation class_graph class_graph_partitioner" ;;
        # ---- L5: + fields and transport. NO differential operator.
        level-5-field-calculus) echo "partitioned_pde_assert graph_carrier graph_grammar class_graph class_graph_field class_graph_partitioner class_graph_assembler" ;;
        # ---- L6: + the differential operator and its fixture.
        #          Still no solver.
        level-6-discretization) echo "partitioned_pde_assert shifted_laplacian_fixture graph_carrier graph_grammar class_graph class_graph_field class_graph_partitioner class_graph_assembler class_graph_differential_operator" ;;
        # ---- L7: + minimization. Consumes Level 6, not Level 8.
        level-7-minimization) echo "partitioned_pde_assert shifted_laplacian_fixture graph_carrier graph_grammar class_graph class_graph_field class_graph_gmres" ;;
        # ---- L8: + the partitioned constitution, which consumes L6.
        level-8-constitution) echo "partitioned_pde_assert shifted_laplacian_fixture partitioned_shifted_laplacian_fixture graph_carrier graph_grammar class_graph class_graph_field class_graph_gmres" ;;
        # ---- L9: the statement, consuming Level 8 and the solver.
        level-9-statement) echo "partitioned_pde_assert shifted_laplacian_fixture partitioned_shifted_laplacian_fixture graph_carrier graph_grammar class_graph class_graph_field class_graph_gmres" ;;

        *)                 echo "__no_allowlist__" ;;
    esac
}

violation=0

for dir in "$here"/common "$here"/level-*; do
    [ -d "$dir" ] || continue
    name="$(basename "$dir")"
    sources=$(ls "$dir"/*.f90 2>/dev/null)
    [ -n "$sources" ] || continue

    for src in $sources; do
        # per-file allowlist first, then the level's
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
echo "import gate: every source imports only its level and below"
