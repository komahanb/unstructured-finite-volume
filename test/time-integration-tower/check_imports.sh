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
#     time_carriers_fixture      earned at Level 0
#     time_relations_fixture     earned at Level 1
#     time_algebra_fixture       earned at Level 2
#     time_fields_fixture        earned at Level 5
#     triangular_decay_fixture   earned at Level 6
#
# The fixture ladder is the tower's own stratification applied to
# itself. Level 0 declares carriers and may reach ONLY the carrier
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
#
# A level that could reach the step operator before Level 6 could
# have redescribed production instead of establishing the meaning
# production is measured against - which is exactly what Level 6's
# RED depended on not having happened.
#
# UNIVERSAL refusals never lift, at any built level:
#
#     class_graph_marcher     no level has earned it yet
#     class_graph_newton      no level has earned it
#     class_graph_linearization + the derivative/adjoint fixtures
#                             this is not a derivative-family client
#
# The marcher is the sharpest of these. Levels 5-7 test one step,
# by hand, twice. Whether the machinery that stamps a step along a
# chain can carry a state domain that is not its host's vertex set
# is a question for a level that does not exist, and importing the
# module before then would answer it by assumption.
#
# Levels 8-9 are unbuilt. Their allowlists do not exist yet, and a
# level directory that appears without one fails closed - which is
# the correct behaviour for a frontier.

here="$(cd "$(dirname "$0")" && pwd)"

intrinsics="iso_fortran_env iso_c_binding ieee_arithmetic ieee_exceptions ieee_features"

allowed_for() {
    case "$1" in
        # ---- shared fixtures, keyed by the level that earns them
        common/time_assert.f90) echo "" ;;
        common/time_carriers_fixture.f90) echo "graph_carrier" ;;
        common/time_relations_fixture.f90) echo "graph_carrier graph_relation graph_binary_relation" ;;
        common/time_algebra_fixture.f90) echo "graph_carrier graph_relation graph_relation_algebra graph_binary_relation" ;;
        common/time_fields_fixture.f90) echo "graph_carrier graph_field_calculus class_graph_field time_assert" ;;
        common/triangular_decay_fixture.f90) echo "graph_carrier graph_grammar graph_field_calculus class_graph_field" ;;
        common)            echo "__no_allowlist__" ;;

        # ---- L0: carriers only. NOTHING relational - not the
        #          relation nucleus, and not the Level-1 fixture.
        level-0-carrier)   echo "time_assert time_carriers_fixture graph_carrier" ;;
        # ---- L1: + the relation nucleus
        level-1-relation)  echo "time_assert time_carriers_fixture time_relations_fixture graph_carrier graph_relation graph_binary_relation" ;;
        # ---- L2: + relation algebra
        level-2-relation-algebra) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture graph_carrier graph_relation graph_relation_algebra graph_binary_relation" ;;
        # ---- L3: + the relational graph container
        level-3-graph)     echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture graph_carrier graph_relation graph_relation_algebra graph_binary_relation graph_structure" ;;
        # ---- L4: + the profile and its algorithms. No marcher.
        level-4-graph-calculus) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture graph_carrier graph_relation graph_relation_algebra graph_binary_relation graph_structure graph_profile graph_algorithms" ;;

        # ===== REVIEW GATE A =====

        # ---- L5: + fields. Values, and nothing that steps or solves.
        level-5-field-calculus) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture graph_carrier graph_relation graph_relation_algebra graph_binary_relation graph_field_calculus class_graph_field" ;;
        # ---- L6: + the ordinary graph (the compatibility host), the
        #          operation contract, and the step operators. NO
        #          minimizer: the scheme is tested before the solve.
        level-6-discretization) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture graph_carrier graph_grammar graph_relation graph_relation_algebra graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_step" ;;
        # ---- L7: + minimization and its gmres concretion. Still no
        #          marcher.
        level-7-minimization) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture graph_carrier graph_grammar graph_relation graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_step graph_minimization class_graph_gmres" ;;

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
    levels="level-0-carrier level-1-relation level-2-relation-algebra level-3-graph level-4-graph-calculus level-5-field-calculus level-6-discretization level-7-minimization"
    before_six="level-0-carrier level-1-relation level-2-relation-algebra level-3-graph level-4-graph-calculus level-5-field-calculus"
    before_seven="$before_six level-6-discretization"

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

    # L0 earns carriers and the carrier fixture, and NOTHING relational.
    permits level-0-carrier time_carriers_fixture
    permits level-0-carrier time_assert
    permits level-0-carrier graph_carrier
    permits level-0-carrier iso_fortran_env
    refuses level-0-carrier time_relations_fixture   # the fixture ladder
    refuses level-0-carrier time_algebra_fixture
    refuses level-0-carrier graph_relation
    refuses level-0-carrier graph_binary_relation
    refuses level-0-carrier graph_structure

    # L1 stands on L0's fixture and adds its own; L2's is still above it.
    permits level-1-relation time_carriers_fixture
    permits level-1-relation time_relations_fixture
    refuses level-1-relation time_algebra_fixture
    refuses level-1-relation graph_structure

    # L2 earns the algebra; the container is still one rung up.
    permits level-2-relation-algebra time_algebra_fixture
    permits level-2-relation-algebra graph_relation_algebra
    refuses level-2-relation-algebra graph_structure

    # L3 earns the container; the profile and its algorithms are not
    # its business.
    permits level-3-graph graph_structure
    refuses level-3-graph graph_profile
    refuses level-3-graph graph_algorithms

    # L4 earns the interpretation.
    permits level-4-graph-calculus graph_profile
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

    # ---- UNIVERSAL: refusals that never lift at any built level.
    #      The marcher is the sharpest: Levels 5-7 solve ONE step by
    #      hand, and whether the machinery that stamps a step along a
    #      chain can carry a non-vertex state domain is a question no
    #      built level has asked.
    for lvl in $levels; do
        refuses "$lvl" class_graph_marcher
        refuses "$lvl" class_graph_newton
        refuses "$lvl" class_graph_linearization
        refuses "$lvl" derivative_fixture
        refuses "$lvl" adjoint_fixture
    done

    # The fixtures themselves are keyed per file.
    permits common/time_carriers_fixture.f90 graph_carrier
    refuses common/time_carriers_fixture.f90 graph_binary_relation
    refuses common/time_assert.f90 graph_carrier
    permits common/triangular_decay_fixture.f90 graph_grammar
    refuses common/triangular_decay_fixture.f90 class_graph_step
    refuses common/triangular_decay_fixture.f90 class_graph_marcher

    # An unclassified source fails closed rather than silently open -
    # which is what the unbuilt Level 8 meets.
    allows level-8-constitution graph_carrier
    if [ "$?" -ne 2 ]; then
        echo " FAIL : an unbuilt level did not fail closed"
        fail=1
    fi

    if [ "$fail" -ne 0 ]; then
        echo "IMPORT GATE: the layering decision is wrong"
        exit 1
    fi
    echo "import gate: the ladder rises one rung at a time, and the marcher is refused at every level"
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
