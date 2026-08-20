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
#     field_stored       earned at Level 5
#     operation_step        earned at Level 6, refused at 0-5
#     operation_minimization      earned at Level 7, refused at 0-6
#     operation_gmres       earned at Level 7, refused at 0-6
#     operation_marching     earned at Level 8, refused at 0-7
#     operation_newton      earned at Level 8, refused at 0-7
#
# A level that could reach the step operator before Level 6 could
# have redescribed production instead of establishing the meaning
# production is measured against - which is exactly what Level 6's
# RED depended on not having happened.
#
# UNIVERSAL refusals never lift, at ANY level of this tower:
#
#     operation_linearization
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
        # 2026-08-16: the relational container is retired. A level reading
        # (S, P) is granted graph_fractal and view_relational,
        # and builds the representation itself. Granted per level, in
        # review; the list is an assertion, not a history.
        # ---- shared fixtures, keyed by the level that earns them
        #
        # 2026-08-16: graph_carrier is retired as a domain source. A
        # domain is now a set GRAPH plus a representation held in a
        # map, so the single carrier grant splits into three, and each
        # source gets only the part it uses:
        #
        #     graph_fractal            WHICH set - identity, always
        #     map_set_representation HOW its members are stored -
        #                              only where one is CONSTRUCTED
        #     map_set            the association - only where a
        #                              domain is described or queried
        #     map_inclusion      provenance - only where an
        #                              embedding is ASSERTED
        #
        # No level constructs a representation: they all take their
        # domains from the carrier fixture, so map_set_representation
        # is granted THERE and nowhere else. And no label map appears
        # anywhere in this tower, because nothing in it asks a domain
        # what it is called.
        common/time_assert.f90) echo "" ;;
        common/time_carriers_fixture.f90) echo "graph_fractal map_set_representation map_set" ;;
        common/time_relations_fixture.f90) echo "graph_fractal map_set relation_finitary relation_binary" ;;
        common/time_algebra_fixture.f90) echo "graph_fractal map_set relation_finitary relation_algebra relation_binary" ;;
        common/time_fields_fixture.f90) echo "graph_fractal field_calculus field_stored time_assert" ;;
        # The action stores an identity and a count, so it needs
        # neither a representation nor a map. Each name comes from
        # its owner: graph from the KERNEL, the directed graph
        # from view_directed, the field from the field
        # calculus, and the operation contract from
        # operation_action. graph_grammar, which once lent all
        # four, is deleted.
        common/triangular_decay_fixture.f90) echo "graph_fractal operation_action view_directed field_calculus field_stored" ;;
        common)            echo "__no_allowlist__" ;;

        # ---- L0: sets only. NOTHING relational - not the relation
        #          nucleus, and not the Level-1 fixture.
        level-0-carrier)   echo "time_assert time_carriers_fixture graph_fractal map_set" ;;
        # ---- L1: + the relation nucleus
        level-1-relation)  echo "time_assert time_carriers_fixture time_relations_fixture graph_fractal map_set relation_finitary relation_binary" ;;
        # ---- L2: + relation algebra
        level-2-relation-algebra) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture graph_fractal map_set relation_finitary relation_algebra relation_binary" ;;
        # ---- L3: + the relational graph container
        level-3-graph)     echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture graph_fractal map_set relation_finitary relation_algebra relation_binary view_relational" ;;
        # ---- L4: + the profile and its algorithms. No marcher.
        #
        #          This is the ONE level granted map_inclusion:
        #          it is the only level that ASSERTS an embedding, and
        #          it asserts it of the sets sources/sinks carve.
        #
        #          It is also the one level granted map_label,
        #          and NOT because it asks a domain its name - nothing
        #          in this tower ever does. sources/sinks CARVE, and a
        #          carve binds extension, label and embedding TOGETHER
        #          so that no half-described set escapes. A caller of
        #          a carving operation must therefore HOLD a label map
        #          even where it will never read one. The grant records
        #          that cost rather than hiding it: held, not queried.
        level-4-graph-calculus) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture graph_fractal map_set map_label map_inclusion relation_finitary relation_algebra relation_binary relation_algorithms view_relational" ;;

        # ===== REVIEW GATE A =====

        # ---- L5: + fields. Values, and nothing that steps or solves.
        level-5-field-calculus) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture graph_fractal map_set relation_finitary relation_algebra relation_binary field_calculus field_stored" ;;
        # ---- L6: + the directed graph (the compatibility host), the
        #          operation contract, and the step operators. NO
        #          minimizer: the scheme is tested before the solve.
        level-6-discretization) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture graph_fractal map_set view_directed relation_finitary relation_algebra relation_binary field_calculus view_directed_stored field_stored operation_step" ;;
        # ---- L7: + minimization and its gmres concretion. Still no
        #          marcher.
        level-7-minimization) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture graph_fractal map_set view_directed relation_finitary relation_binary field_calculus view_directed_stored field_stored operation_step operation_minimization operation_gmres" ;;

        # ===== REVIEW GATE B =====

        # ---- L8: + the marcher and newton, the constituted citizens
        #          under test. NOT operation_linearization: newton
        #          reaches it, and the tower may not.
        level-8-constitution) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture graph_fractal map_set relation_finitary relation_algebra relation_binary field_calculus view_directed_stored field_stored operation_step operation_minimization operation_gmres operation_newton operation_marching" ;;
        # ---- L9: the statement, on the same constitution.
        level-9-statement) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture graph_fractal map_set relation_finitary relation_binary field_calculus view_directed_stored field_stored operation_step operation_minimization operation_gmres operation_newton operation_marching" ;;

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
    levels="level-0-carrier level-1-relation level-2-relation-algebra level-3-graph level-4-graph-calculus level-5-field-calculus level-6-discretization level-7-minimization level-8-constitution level-9-statement"
    before_six="level-0-carrier level-1-relation level-2-relation-algebra level-3-graph level-4-graph-calculus level-5-field-calculus"
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

    # L0 earns carriers and the carrier fixture, and NOTHING relational.
    permits level-0-carrier time_carriers_fixture
    permits level-0-carrier time_assert
    permits level-0-carrier graph_fractal
    permits level-0-carrier map_set
    permits level-0-carrier iso_fortran_env
    refuses level-0-carrier graph_carrier            # the retired source
    # A level takes its domains from the fixture, so it never builds a
    # representation and is not granted the module that makes one.
    refuses level-0-carrier map_set_representation
    # Provenance is asserted at Level 4 alone.
    refuses level-0-carrier map_inclusion
    refuses level-0-carrier time_relations_fixture   # the fixture ladder
    refuses level-0-carrier time_algebra_fixture
    refuses level-0-carrier relation_finitary
    refuses level-0-carrier relation_binary
    refuses level-0-carrier view_relational

    # L1 stands on L0's fixture and adds its own; L2's is still above it.
    permits level-1-relation time_carriers_fixture
    permits level-1-relation time_relations_fixture
    refuses level-1-relation time_algebra_fixture
    refuses level-1-relation view_relational

    # L2 earns the algebra; the container is still one rung up.
    permits level-2-relation-algebra time_algebra_fixture
    permits level-2-relation-algebra relation_algebra
    refuses level-2-relation-algebra view_relational

    # L3 earns the container; the profile and its algorithms are not
    # its business.
    permits level-3-graph graph_fractal
    permits level-3-graph view_relational
    refuses level-3-graph graph_profile
    refuses level-3-graph relation_algorithms

    # L4 earns the interpretation.
    refuses level-4-graph-calculus graph_profile
    permits level-4-graph-calculus relation_algorithms

    # ---- STAGED: ceilings that rise, each at the level that earns it.
    permits level-5-field-calculus time_fields_fixture
    permits level-5-field-calculus field_stored
    refuses level-4-graph-calculus field_stored
    refuses level-4-graph-calculus time_fields_fixture

    permits level-6-discretization operation_step
    permits level-6-discretization triangular_decay_fixture
    for lvl in $before_six; do
        refuses "$lvl" operation_step
        refuses "$lvl" triangular_decay_fixture
    done

    permits level-7-minimization operation_minimization
    permits level-7-minimization operation_gmres
    for lvl in $before_seven; do
        refuses "$lvl" operation_minimization
        refuses "$lvl" operation_gmres
    done

    permits level-8-constitution operation_marching
    permits level-8-constitution operation_newton
    permits level-9-statement operation_marching
    for lvl in $before_eight; do
        refuses "$lvl" operation_marching
        refuses "$lvl" operation_newton
    done

    # ---- UNIVERSAL: refusals that never lift at ANY level.
    #      operation_linearization is the load-bearing one: Level 8
    #      reached its Class-2 defect through
    #      marcher -> newton -> difference_linearization, and because
    #      no level may name that module, the failure cannot have been
    #      manufactured.
    for lvl in $levels; do
        refuses "$lvl" operation_linearization
        refuses "$lvl" derivative_fixture
        refuses "$lvl" adjoint_fixture
    done

    # The fixtures themselves are keyed per file. The carrier fixture
    # is the ONE source that builds a representation, so it is the one
    # source granted the module that makes them.
    permits common/time_carriers_fixture.f90 graph_fractal
    permits common/time_carriers_fixture.f90 map_set_representation
    permits common/time_carriers_fixture.f90 map_set
    refuses common/time_carriers_fixture.f90 relation_binary
    refuses common/time_carriers_fixture.f90 map_inclusion
    refuses common/time_assert.f90 graph_fractal
    # The action holds an identity and a count, and asks its domain
    # nothing else - so no map reaches it, at any rung.
    refuses common/triangular_decay_fixture.f90 map_set
    refuses common/triangular_decay_fixture.f90 map_set_representation
    # Provenance is Level 4's alone: no other level may assert an
    # embedding, and no fixture may.
    refuses level-5-field-calculus map_inclusion
    refuses level-9-statement map_inclusion
    # The label map is Level 4's alone, and only because carving
    # compels holding one - no other level carves, so no other level
    # may hold it.
    permits level-4-graph-calculus map_label
    refuses level-0-carrier map_label
    refuses level-5-field-calculus map_label
    refuses level-9-statement map_label
    permits common/triangular_decay_fixture.f90 operation_action
    refuses common/triangular_decay_fixture.f90 operation_step
    refuses common/triangular_decay_fixture.f90 operation_marching

    # An unclassified source still fails closed rather than silently
    # open - the ten levels are named, and nothing else is.
    allows level-10-nowhere graph_fractal
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
