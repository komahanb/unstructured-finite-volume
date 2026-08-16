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
        # (S, P) is granted fractal_graph and graph_relational_view,
        # and builds the representation itself. Granted per level, in
        # review; the list is an assertion, not a history.
        common)                  echo "" ;;
        level-0-carrier)         echo "learning_assert fractal_graph graph_set_representation graph_set_map" ;;
        level-1-relation)        echo "learning_assert fractal_graph graph_set_representation graph_set_map graph_relation" ;;
        # level 2: the algebra, and NOT the binary storage - D is
        # held as class(relation).
        level-2-relation-algebra) echo "learning_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra" ;;
        # level 5: values need domains, not graphs - the smallest
        # allowlist of any rung above the ground.
        level-5-field-calculus)  echo "learning_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map class_graph_field" ;;
        # level 4: the profile's interpretation and the algorithms
        # that walk it - the binary storage the view leans on stays
        # forbidden to the learning client.
        level-4-graph-calculus)  echo "learning_assert fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_relation graph_relation_algebra graph_profile graph_algorithms fractal_graph graph_relational_view" ;;
        # level 3: the container; graph_binary_relation is granted
        # for the view refusal ONLY - the production path never
        # touches it.
        level-3-graph)           echo "learning_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation fractal_graph graph_relational_view" ;;
        # level 6: the structural-Jacobian rung. The binary citizen
        # is earned at last - J_Theta must materialize its generated
        # pairs and answer its reverse as a transpose view. Fields
        # stay forbidden: dependency structure belongs to the model,
        # not to one data instance.
        level-6-discretization)  echo "learning_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation graph_profile graph_algorithms fractal_graph graph_relational_view" ;;
        # level 7: the fitting rung. The GMRES citizen inherits
        # attach/constant/solve from the minimizer base, so
        # graph_minimization is not directly imported and stays off
        # the list; the relation/algebra/profile stack is not needed
        # at all - the solver sees an opaque R : Theta -> Y. The
        # named fixture is the level's own test-local oracle.
        level-7-minimization)    echo "learning_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_grammar graph_field_calculus class_graph_field class_graph class_graph_gmres learning_residual_fixture" ;;
        # level 8: the meaning rung. Structure derives the order,
        # the test-local constitution supplies the laws, L supplies
        # the home - no solver, no legacy host, no minimization, no
        # binary storage. The named fixture is the level's own.
        level-8-constitution)    echo "learning_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra graph_profile graph_algorithms class_graph_field learning_constitution_fixture fractal_graph graph_relational_view" ;;
        # level 9: the statement - the composition rung. The REUSED
        # level-8 constitution and the level's own adapter are the
        # named fixtures; Level 7's affine oracle is deliberately
        # absent, and graph_minimization stays off: gmres inherits
        # the minimizer face. No new mathematics enters here.
        level-9-statement)       echo "learning_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra graph_profile graph_algorithms graph_grammar graph_field_calculus class_graph class_graph_field class_graph_gmres learning_constitution_fixture constituted_residual_fixture fractal_graph graph_relational_view" ;;
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
