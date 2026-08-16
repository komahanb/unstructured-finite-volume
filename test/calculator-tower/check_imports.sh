#!/bin/bash
# The import gate: the dependency-stratification law, made
# mechanical. Every calculator source may `use` only the framework
# modules its level has been explicitly granted; the allowlist
# grows one level at a time, in review, never by accident. A level
# directory with sources but no allowlist is itself a violation -
# the list evolves explicitly or not at all.
#
# This gate audits the CALCULATOR TESTS' imports only. It does not
# analyze the production dependency graph.

here="$(cd "$(dirname "$0")" && pwd)"

# Fortran intrinsic modules are always permitted.
intrinsics="iso_fortran_env iso_c_binding ieee_arithmetic ieee_exceptions ieee_features"

allowed_for() {
    case "$1" in
        # 2026-08-16: the relational container is retired. A level reading
        # (S, P) is granted fractal_graph and graph_relational_view,
        # and builds the representation itself. Granted per level, in
        # review; the list is an assertion, not a history.
        common)                  echo "" ;;
        level-0-carrier)         echo "calculator_assert fractal_graph graph_set_representation graph_set_map" ;;
        level-1-relation)        echo "calculator_assert fractal_graph graph_set_representation graph_set_map graph_relation" ;;
        level-2-relation-algebra) echo "calculator_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra" ;;
        # level 3: graph_binary_relation is granted for the view
        # refusal alone (transpose_of is level-1 relation
        # infrastructure); the profile and everything above stay
        # forbidden.
        level-3-graph)           echo "calculator_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation fractal_graph graph_relational_view" ;;
        # level 4: the profile's interpretation and the algorithms
        # that walk it; the legacy graph stack, fields and all above
        # stay forbidden.
        level-4-graph-calculus)  echo "calculator_assert fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_relation graph_relation_algebra graph_profile graph_algorithms fractal_graph graph_relational_view" ;;
        # level 6: sparsity is structure - carriers, relations, the
        # algebra, and the binary transpose view; no field, no graph
        # container, no value ever needed.
        level-6-discretization)  echo "calculator_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation" ;;
        # level 7: the ordinary minimization machinery and the legacy
        # operation-host compatibility it still rides; no constitution,
        # no statement, no flow relation, no operator meanings.
        level-7-minimization)    echo "calculator_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_grammar graph_ordinary_view graph_field_calculus class_graph_field class_graph graph_minimization class_graph_gmres affine_residual_fixture" ;;
        # level 8: constitution binds meaning to domains, relations
        # and fields - no solver, no graph host, no statement.
        level-8-constitution)    echo "calculator_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra class_graph_field arithmetic_constitution_fixture" ;;
        # level 9: the composition rung - everything the statement
        # SELECTS, and nothing it would have to invent.
        level-9-statement)       echo "calculator_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra graph_field_calculus class_graph_field graph_grammar graph_ordinary_view class_graph graph_minimization class_graph_gmres arithmetic_constitution_fixture constituted_residual_fixture fractal_graph graph_relational_view" ;;
        # level 5: a field needs a domain, never a graph container.
        level-5-field-calculus)  echo "calculator_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_field_calculus class_graph_field" ;;
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
        # Every `use` statement, tolerant of `use, intrinsic ::` and
        # `use module, only : ...` spellings; case-insensitive.
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
echo "import gate: every calculator source imports only its level and below"
