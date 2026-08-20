#!/bin/bash
# The learning tower's import gate: the dependency-stratification
# law, made mechanical - a clone of the calculator gate's
# philosophy, with its own allowlists. Every learning source may
# `use` only the framework modules its level has been explicitly
# granted; a directory with sources but no allowlist fails closed.
# This gate audits the LEARNING TESTS' imports only.

here="$(cd "$(dirname "$0")" && pwd)"

intrinsics="iso_fortran_env iso_c_binding ieee_arithmetic ieee_exceptions ieee_features"

allowed_for() {
    case "$1" in
        # 2026-08-16: the relational container is retired. A level reading
        # (S, P) is granted graph_fractal and view_relational,
        # and builds the representation itself. Granted per level, in
        # review; the list is an assertion, not a history.
        common)                  echo "" ;;
        level-0-carrier)         echo "learning_assert graph_fractal map_set_representation map_set" ;;
        level-1-relation)        echo "learning_assert graph_fractal map_set_representation map_set relation_finitary" ;;
        # level 2: the algebra, and NOT the binary storage - D is
        # held as class(relation).
        level-2-relation-algebra) echo "learning_assert graph_fractal map_set_representation map_set map_inclusion relation_finitary relation_algebra" ;;
        # level 5: values need domains, not graphs - the smallest
        # allowlist of any rung above the ground.
        level-5-field-calculus)  echo "learning_assert graph_fractal map_set_representation map_set map_inclusion field_stored" ;;
        # level 4: the profile's interpretation and the algorithms
        # that walk it - the binary storage the view leans on stays
        # forbidden to the learning client.
        level-4-graph-calculus)  echo "learning_assert graph_fractal map_set_representation map_set map_label map_inclusion relation_finitary relation_algebra relation_algorithms graph_fractal view_relational relation_binary" ;;
        # level 3: the container; relation_binary is granted
        # for the view refusal ONLY - the production path never
        # touches it.
        level-3-graph)           echo "learning_assert graph_fractal map_set_representation map_set map_inclusion relation_finitary relation_algebra relation_binary graph_fractal view_relational" ;;
        # level 6: the structural-Jacobian rung. The binary citizen
        # is earned at last - J_Theta must materialize its generated
        # pairs and answer its reverse as a transpose view. Fields
        # stay forbidden: dependency structure belongs to the model,
        # not to one data instance.
        level-6-discretization)  echo "learning_assert graph_fractal map_set_representation map_set map_inclusion relation_finitary relation_algebra relation_binary relation_algorithms graph_fractal view_relational" ;;
        # level 7: the fitting rung. The GMRES citizen inherits
        # attach/constant/solve from the minimizer base, so
        # graph_minimization is not directly imported and stays off
        # the list; the relation/algebra/profile stack is not needed
        # at all - the solver sees an opaque R : Theta -> Y. The
        # named fixture is the level's own test-local oracle.
        level-7-minimization)    echo "learning_assert graph_fractal map_set_representation map_set map_inclusion operation_action view_directed field_calculus field_stored view_directed_stored class_graph_gmres learning_residual_fixture" ;;
        # level 8: the meaning rung. Structure derives the order,
        # the test-local constitution supplies the laws, L supplies
        # the home - no solver, no legacy host, no minimization, no
        # binary storage. The named fixture is the level's own.
        level-8-constitution)    echo "learning_assert graph_fractal map_set_representation map_set map_inclusion relation_finitary relation_algebra relation_algorithms field_stored learning_constitution_fixture graph_fractal view_relational" ;;
        # level 9: the statement - the composition rung. The REUSED
        # level-8 constitution and the level's own adapter are the
        # named fixtures; Level 7's affine oracle is deliberately
        # absent, and graph_minimization stays off: gmres inherits
        # the minimizer face. No new mathematics enters here.
        level-9-statement)       echo "learning_assert graph_fractal map_set_representation map_set map_inclusion relation_finitary relation_algebra relation_algorithms operation_action view_directed field_calculus view_directed_stored field_stored class_graph_gmres learning_constitution_fixture constituted_residual_fixture graph_fractal view_relational" ;;
        *)                       echo "__no_allowlist__" ;;
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
echo "import gate: every learning source imports only its level and below"
