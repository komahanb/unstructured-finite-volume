#!/bin/bash
# The time integration tower's import gate, keyed PER LEVEL.
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
#     time_assert                below everything (nothing)
#     time_sets_fixture      earned at Level 0
#     time_relations_fixture     earned at Level 1
#     time_algebra_fixture       earned at Level 2
#     time_fields_fixture        earned at Level 5
#     triangular_decay_fixture   earned at Level 6
#
# The fixture ladder is the tower's own stratification applied to
# itself. Level 0 declares sets and may reach ONLY the set
# fixture; the relation fixture is one rung above it. A Level-0
# source that says
#
#     use time_relations_fixture
#
# is a layering violation and this gate must refuse it. That exact
# leak was found and closed in the partitioned tower; --selftest
# asserts it here from the first commit rather than after a review.
#
#                  THE FRONTIER, HELD SHUT
#
# Two kinds of refusal live here, and the difference matters.
#
# STAGED refusals are ceilings that RISE. Levels 0-4 established
# what time IS independently of how this repository marches through
# it, so they see no machinery at all; the machinery then arrives
# one rung at a time, each where its level earns it:
#
#     class_graph_field       earned at Level 5
#     class_graph_step        earned at Level 6, refused at 0-5
#     graph_minimization      earned at Level 7, refused at 0-6
#     class_graph_gmres       earned at Level 7, refused at 0-6
#     class_graph_marcher     earned at Level 8, refused at 0-7
#     class_graph_newton      earned at Level 8, refused at 0-7
#
# A level that could reach the step operator before Level 6 could
# have redescribed production instead of establishing the meaning
# production is measured against - which is exactly what Level 6's
# RED depended on not having happened.
#
# UNIVERSAL refusals never lift, at ANY level of this tower:
#
#     class_graph_linearization
#     the derivative and adjoint fixtures
#
# THE LINEARIZATION REFUSAL IS THE SHARPEST ASSERTION IN THIS GATE,
# and it is load-bearing evidence rather than hygiene.
#
# Level 8 marches implicitly through
#
#     marcher -> newton -> difference_linearization -> gmres
#
# and difference_linearization FAILED there, on the state domain,
# exactly as the reverse review's Class-2 seam predicted. Because
# no level of this tower may name that module, the failure cannot
# have been manufactured: it was reached by the production call
# chain an implicit march requires, and by nothing else. The gate
# is what makes that claim checkable instead of promised.
#
# The derivative and adjoint fixtures stay refused for the older
# reason: this is not a derivative-family client, and its
# independence as evidence depends on that being mechanical.

here="$(cd "$(dirname "$0")" && pwd)"

intrinsics="iso_fortran_env iso_c_binding ieee_arithmetic ieee_exceptions ieee_features"

allowed_for() {
    case "$1" in
        # ---- shared fixtures, keyed by the level that earns them
        common/time_assert.f90) echo "" ;;
        common/time_sets_fixture.f90) echo "graph_set" ;;
        common/time_relations_fixture.f90) echo "graph_set graph_relation graph_binary_relation" ;;
        common/time_algebra_fixture.f90) echo "graph_set graph_relation graph_relation_algebra graph_binary_relation" ;;
        common/time_fields_fixture.f90) echo "graph_set graph_field_calculus class_graph_field time_assert" ;;
        common/triangular_decay_fixture.f90) echo "graph_set graph_grammar graph_field_calculus class_graph_field" ;;
        common)            echo "__no_allowlist__" ;;

        # ---- L0: sets only. NOTHING relational - not the
        #          relation nucleus, and not the Level-1 fixture.
        level-0-set)   echo "time_assert time_sets_fixture graph_set" ;;
        # ---- L1: + the relation nucleus
        level-1-relation)  echo "time_assert time_sets_fixture time_relations_fixture graph_set graph_relation graph_binary_relation" ;;
        # ---- L2: + relation algebra
        level-2-relation-algebra) echo "time_assert time_sets_fixture time_relations_fixture time_algebra_fixture graph_set graph_relation graph_relation_algebra graph_binary_relation" ;;
        # ---- L3: + the related graph container
        level-3-graph)     echo "time_assert time_sets_fixture time_relations_fixture time_algebra_fixture graph_set graph_relation graph_relation_algebra graph_binary_relation graph_structure" ;;
        # ---- L4: + the interpretation and its algorithms. No marcher.
        level-4-graph-calculus) echo "time_assert time_sets_fixture time_relations_fixture time_algebra_fixture graph_set graph_relation graph_relation_algebra graph_binary_relation graph_structure graph_interpretation graph_algorithms" ;;

        # ===== REVIEW GATE A =====

        # ---- L5: + fields. Values, and nothing that steps or solves.
        level-5-field-calculus) echo "time_assert time_sets_fixture time_relations_fixture time_algebra_fixture time_fields_fixture graph_set graph_relation graph_relation_algebra graph_binary_relation graph_field_calculus class_graph_field" ;;
        # ---- L6: + the ordinary graph (the compatibility host), the
        #          operation contract, and the step operators. NO
        #          minimizer: the scheme is tested before the solve.
        level-6-discretization) echo "time_assert time_sets_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture graph_set graph_grammar graph_relation graph_relation_algebra graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_step" ;;
        # ---- L7: + minimization and its gmres concretion. Still no
        #          marcher.
        level-7-minimization) echo "time_assert time_sets_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture graph_set graph_grammar graph_relation graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_step graph_minimization class_graph_gmres" ;;

        # ===== REVIEW GATE B =====

        # ---- L8: + the marcher and newton, the constituted citizens
        #          under test. NOT class_graph_linearization: newton
        #          reaches it, and the tower may not.
        level-8-constitution) echo "time_assert time_sets_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture graph_set graph_grammar graph_relation graph_relation_algebra graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_step graph_minimization class_graph_gmres class_graph_newton class_graph_marcher" ;;
        # ---- L9: the statement, on the same constitution.
        level-9-statement) echo "time_assert time_sets_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture graph_set graph_grammar graph_relation graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_step graph_minimization class_graph_gmres class_graph_newton class_graph_marcher" ;;

        *)                 echo "__no_allowlist__" ;;
    esac
}

# The per-file key wins; the level's is the fallback. A source whose
# level has no allowlist either fails closed.
allowlist_for() {
    local allow
    allow="$(allowed_for "$1")"
    if [ "$allow" = "__no_allowlist__" ]; then
        allow="$(allowed_for "${1%%/*}")"
    fi
    echo "$allow"
}

# allows <level[/file]> <module>
#     0  the module is on that source's ceiling
#     1  it is not - a layering violation
#     2  the source has no declared allowlist at all
allows() {
    local allow a
    allow="$(allowlist_for "$1")"
    [ "$allow" = "__no_allowlist__" ] && return 2
    for a in $allow $intrinsics; do
        [ "$2" = "$a" ] && return 0
    done
    return 1
}

#---------------------------------------------------------------------
# The gate's own test: the decision function, exercised on the two
# questions this tower must never get wrong - the fixture ladder,
# and the frontier. This proves the ALLOWLISTS say what they must;
# the bare scan below proves the scanner acts on them.
#---------------------------------------------------------------------

if [ "$1" = "--selftest" ]; then
    fail=0
    levels="level-0-set level-1-relation level-2-relation-algebra level-3-graph level-4-graph-calculus level-5-field-calculus level-6-discretization level-7-minimization level-8-constitution level-9-statement"
    before_six="level-0-set level-1-relation level-2-relation-algebra level-3-graph level-4-graph-calculus level-5-field-calculus"
    before_seven="$before_six level-6-discretization"
    before_eight="$before_seven level-7-minimization"

    permits() {
        if allows "$1" "$2"; then :; else
            echo " FAIL : the import gate refused '$2' at $1"
            fail=1
        fi
    }
    refuses() {
        if allows "$1" "$2"; then
            echo " FAIL : the import gate permitted '$2' at $1"
            fail=1
        fi
    }

    # L0 earns sets and the set fixture, and NOTHING relational.
    permits level-0-set time_sets_fixture
    permits level-0-set time_assert
    permits level-0-set graph_set
    permits level-0-set iso_fortran_env
    refuses level-0-set time_relations_fixture   # the fixture ladder
    refuses level-0-set time_algebra_fixture
    refuses level-0-set graph_relation
    refuses level-0-set graph_binary_relation
    refuses level-0-set graph_structure

    # L1 stands on L0's fixture and adds its own; L2's is still above it.
    permits level-1-relation time_sets_fixture
    permits level-1-relation time_relations_fixture
    refuses level-1-relation time_algebra_fixture
    refuses level-1-relation graph_structure

    # L2 earns the algebra; the container is still one rung up.
    permits level-2-relation-algebra time_algebra_fixture
    permits level-2-relation-algebra graph_relation_algebra
    refuses level-2-relation-algebra graph_structure

    # L3 earns the container; the interpretation and its algorithms are not
    # its business.
    permits level-3-graph graph_structure
    refuses level-3-graph graph_interpretation
    refuses level-3-graph graph_algorithms

    # L4 earns the interpretation.
    permits level-4-graph-calculus graph_interpretation
    permits level-4-graph-calculus graph_algorithms

    # ---- STAGED: ceilings that rise, each at the level that earns it.
    permits level-5-field-calculus time_fields_fixture
    permits level-5-field-calculus class_graph_field
    refuses level-4-graph-calculus class_graph_field
    refuses level-4-graph-calculus time_fields_fixture

    permits level-6-discretization class_graph_step
    permits level-6-discretization triangular_decay_fixture
    for lvl in $before_six; do
        refuses "$lvl" class_graph_step
        refuses "$lvl" triangular_decay_fixture
    done

    permits level-7-minimization graph_minimization
    permits level-7-minimization class_graph_gmres
    for lvl in $before_seven; do
        refuses "$lvl" graph_minimization
        refuses "$lvl" class_graph_gmres
    done

    permits level-8-constitution class_graph_marcher
    permits level-8-constitution class_graph_newton
    permits level-9-statement class_graph_marcher
    for lvl in $before_eight; do
        refuses "$lvl" class_graph_marcher
        refuses "$lvl" class_graph_newton
    done

    # ---- UNIVERSAL: refusals that never lift at ANY level.
    #      class_graph_linearization is the load-bearing one: Level 8
    #      reached its Class-2 defect through
    #      marcher -> newton -> difference_linearization, and because
    #      no level may name that module, the failure cannot have been
    #      manufactured.
    for lvl in $levels; do
        refuses "$lvl" class_graph_linearization
        refuses "$lvl" derivative_fixture
        refuses "$lvl" adjoint_fixture
    done

    # The fixtures themselves are keyed per file.
    permits common/time_sets_fixture.f90 graph_set
    refuses common/time_sets_fixture.f90 graph_binary_relation
    refuses common/time_assert.f90 graph_set
    permits common/triangular_decay_fixture.f90 graph_grammar
    refuses common/triangular_decay_fixture.f90 class_graph_step
    refuses common/triangular_decay_fixture.f90 class_graph_marcher

    # An unclassified source still fails closed rather than silently
    # open - the ten levels are named, and nothing else is.
    allows level-10-nowhere graph_set
    if [ "$?" -ne 2 ]; then
        echo " FAIL : an undeclared level did not fail closed"
        fail=1
    fi

    if [ "$fail" -ne 0 ]; then
        echo "IMPORT GATE: the layering decision is wrong"
        exit 1
    fi
    echo "import gate: the ladder rises one rung at a time, and the linearization is refused at every level"
    exit 0
fi

violation=0

for dir in "$here"/common "$here"/level-*; do
    [ -d "$dir" ] || continue
    name="$(basename "$dir")"
    sources=$(ls "$dir"/*.f90 2>/dev/null)
    [ -n "$sources" ] || continue

    for src in $sources; do
        key="$name/$(basename "$src")"
        if [ "$(allowlist_for "$key")" = "__no_allowlist__" ]; then
            echo "IMPORT GATE: $(basename "$src") in $name has no declared allowlist"
            violation=1
            continue
        fi

        uses=$(grep -ihE '^[[:space:]]*use[[:space:]]' "$src" \
               | sed -E 's/^[[:space:]]*[uU][sS][eE][[:space:]]*(,[[:space:]]*[iI][nN][tT][rR][iI][nN][sS][iI][cC][[:space:]]*::)?[[:space:]]*([a-zA-Z][a-zA-Z0-9_]*).*/\2/' \
               | tr 'A-Z' 'a-z' | sort -u)
        for mod in $uses; do
            if ! allows "$key" "$mod"; then
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
