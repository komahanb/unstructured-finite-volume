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
        # 2026-08-16: the relational container is retired. A level reading
        # (S, P) is granted fractal_graph and graph_relational_view,
        # and builds the representation itself. Granted per level, in
        # review; the list is an assertion, not a history.
        common)                   echo "" ;;
        level-0-carrier)          echo "derivative_assert fractal_graph graph_set_representation graph_set_map" ;;
        level-1-relation)         echo "derivative_assert fractal_graph graph_set_representation graph_set_map graph_relation" ;;
        level-2-relation-algebra) echo "derivative_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra" ;;
        level-3-graph)            echo "derivative_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra fractal_graph graph_relational_view" ;;
        level-4-graph-calculus)   echo "derivative_assert fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_relation graph_relation_algebra graph_algorithms fractal_graph graph_relational_view graph_binary_relation" ;;
        # level 5: primal values need domains, not graphs - and no
        # tangent/cotangent/seed type exists to be imported.
        level-5-field-calculus)   echo "derivative_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map class_graph_field" ;;
        # level 6: the structural rung - the algebra derives value
        # dependency, the profile walks it, and the binary citizen
        # materializes J_ZX and answers its reverse as a view.
        # Fields stay forbidden: the pattern needs no numbers.
        level-6-derivative-structure) echo "derivative_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation graph_algorithms fractal_graph graph_relational_view" ;;
        # level 8 (Gate B): the numerical action rung. Structure
        # derives the order, the test-local constitution supplies
        # primal laws and ONE local linearization per operation,
        # fields carry seeds and results. NO binary storage - the
        # J-pattern is support metadata, never the propagation
        # itinerary - and no solver, ever, at this gate.
        level-8-derivative-constitution) echo "derivative_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra graph_algorithms class_graph_field derivative_constitution_fixture fractal_graph graph_relational_view" ;;
        # level 9 (Gate C): the statement - the composition rung.
        # The REUSED level-8 constitution is the only fixture; no
        # adapter exists because nothing here must satisfy a legacy
        # operation face. No new mathematics, and still no solver.
        level-9-statement)        echo "derivative_assert fractal_graph graph_set_representation graph_set_map graph_inclusion_map graph_relation graph_relation_algebra graph_algorithms class_graph_field derivative_constitution_fixture fractal_graph graph_relational_view" ;;
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
