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
        # 2026-08-16: the relational container is retired. A level reading
        # (S, P) is granted fractal_graph and graph_relational_view,
        # and builds the representation itself. Granted per level, in
        # review; the list is an assertion, not a history.
        # ---- shared fixtures, keyed by the level that earns them
        #
        # 2026-08-16: graph_carrier is retired as a domain source. A
        # domain is now a set GRAPH plus a representation held in a
        # map, so the single carrier grant splits into three, and each
        # source gets only the part it uses:
        #
        #     fractal_graph            WHICH set - identity, always
        #     graph_set_representation HOW its members are stored -
        #                              only where one is CONSTRUCTED
        #     graph_set_map            the association - only where a
        #                              domain is described or queried
        #     graph_inclusion_map      provenance - only where an
        #                              embedding is ASSERTED
        #
        # No level constructs a representation: they all take their
        # domains from the carrier fixture, so graph_set_representation
        # is granted THERE and nowhere else. And no label map appears
        # anywhere in this tower, because nothing in it asks a domain
        # what it is called.
        common/time_assert.f90) echo "" ;;
        common/time_carriers_fixture.f90) echo "fractal_graph graph_set_representation graph_set_map" ;;
        common/time_relations_fixture.f90) echo "fractal_graph graph_set_map graph_relation graph_binary_relation" ;;
        common/time_algebra_fixture.f90) echo "fractal_graph graph_set_map graph_relation graph_relation_algebra graph_binary_relation" ;;
        common/time_fields_fixture.f90) echo "fractal_graph graph_field_calculus class_graph_field time_assert" ;;
        # The action stores an identity and a count, so it needs
        # neither a representation nor a map. Each name comes from
        # its owner: set_graph from the KERNEL, the directed graph
        # from graph_directed_view, the field from the field
        # calculus, and the operation contract from
        # graph_operation_view. graph_grammar, which once lent all
        # four, is deleted.
        common/triangular_decay_fixture.f90) echo "fractal_graph graph_operation_view graph_directed_view graph_field_calculus class_graph_field" ;;
        common)            echo "__no_allowlist__" ;;

        # ---- L0: sets only. NOTHING relational - not the relation
        #          nucleus, and not the Level-1 fixture.
        level-0-carrier)   echo "time_assert time_carriers_fixture fractal_graph graph_set_map" ;;
        # ---- L1: + the relation nucleus
        level-1-relation)  echo "time_assert time_carriers_fixture time_relations_fixture fractal_graph graph_set_map graph_relation graph_binary_relation" ;;
        # ---- L2: + relation algebra
        level-2-relation-algebra) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture fractal_graph graph_set_map graph_relation graph_relation_algebra graph_binary_relation" ;;
        # ---- L3: + the relational graph container
        level-3-graph)     echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture fractal_graph graph_set_map graph_relation graph_relation_algebra graph_binary_relation graph_relational_view" ;;
        # ---- L4: + the profile and its algorithms. No marcher.
        #
        #          This is the ONE level granted graph_inclusion_map:
        #          it is the only level that ASSERTS an embedding, and
        #          it asserts it of the sets sources/sinks carve.
        #
        #          It is also the one level granted graph_label_map,
        #          and NOT because it asks a domain its name - nothing
        #          in this tower ever does. sources/sinks CARVE, and a
        #          carve binds extension, label and embedding TOGETHER
        #          so that no half-described set escapes. A caller of
        #          a carving operation must therefore HOLD a label map
        #          even where it will never read one. The grant records
        #          that cost rather than hiding it: held, not queried.
        level-4-graph-calculus) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture fractal_graph graph_set_map graph_label_map graph_inclusion_map graph_relation graph_relation_algebra graph_binary_relation graph_algorithms graph_relational_view" ;;

        # ===== REVIEW GATE A =====

        # ---- L5: + fields. Values, and nothing that steps or solves.
        level-5-field-calculus) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture fractal_graph graph_set_map graph_relation graph_relation_algebra graph_binary_relation graph_field_calculus class_graph_field" ;;
        # ---- L6: + the directed graph (the compatibility host), the
        #          operation contract, and the step operators. NO
        #          minimizer: the scheme is tested before the solve.
        level-6-discretization) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture fractal_graph graph_set_map graph_directed_view graph_relation graph_relation_algebra graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_step" ;;
        # ---- L7: + minimization and its gmres concretion. Still no
        #          marcher.
        level-7-minimization) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture fractal_graph graph_set_map graph_directed_view graph_relation graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_step graph_minimization class_graph_gmres" ;;

        # ===== REVIEW GATE B =====

        # ---- L8: + the marcher and newton, the constituted citizens
        #          under test. NOT class_graph_linearization: newton
        #          reaches it, and the tower may not.
        level-8-constitution) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture fractal_graph graph_set_map graph_relation graph_relation_algebra graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_step graph_minimization class_graph_gmres class_graph_newton class_graph_marcher" ;;
        # ---- L9: the statement, on the same constitution.
        level-9-statement) echo "time_assert time_carriers_fixture time_relations_fixture time_algebra_fixture time_fields_fixture triangular_decay_fixture fractal_graph graph_set_map graph_relation graph_binary_relation graph_field_calculus class_graph class_graph_field class_graph_step graph_minimization class_graph_gmres class_graph_newton class_graph_marcher" ;;

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
    permits level-0-carrier fractal_graph
    permits level-0-carrier graph_set_map
    permits level-0-carrier iso_fortran_env
    refuses level-0-carrier graph_carrier            # the retired source
    # A level takes its domains from the fixture, so it never builds a
    # representation and is not granted the module that makes one.
    refuses level-0-carrier graph_set_representation
    # Provenance is asserted at Level 4 alone.
    refuses level-0-carrier graph_inclusion_map
    refuses level-0-carrier time_relations_fixture   # the fixture ladder
    refuses level-0-carrier time_algebra_fixture
    refuses level-0-carrier graph_relation
    refuses level-0-carrier graph_binary_relation
    refuses level-0-carrier graph_relational_view

    # L1 stands on L0's fixture and adds its own; L2's is still above it.
    permits level-1-relation time_carriers_fixture
    permits level-1-relation time_relations_fixture
    refuses level-1-relation time_algebra_fixture
    refuses level-1-relation graph_relational_view

    # L2 earns the algebra; the container is still one rung up.
    permits level-2-relation-algebra time_algebra_fixture
    permits level-2-relation-algebra graph_relation_algebra
    refuses level-2-relation-algebra graph_relational_view

    # L3 earns the container; the profile and its algorithms are not
    # its business.
    permits level-3-graph fractal_graph
    permits level-3-graph graph_relational_view
    refuses level-3-graph
    refuses level-3-graph graph_algorithms

    # L4 earns the interpretation.
    refuses level-4-graph-calculus graph_profile
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

    # The fixtures themselves are keyed per file. The carrier fixture
    # is the ONE source that builds a representation, so it is the one
    # source granted the module that makes them.
    permits common/time_carriers_fixture.f90 fractal_graph
    permits common/time_carriers_fixture.f90 graph_set_representation
    permits common/time_carriers_fixture.f90 graph_set_map
    refuses common/time_carriers_fixture.f90 graph_binary_relation
    refuses common/time_carriers_fixture.f90 graph_inclusion_map
    refuses common/time_assert.f90 fractal_graph
    # The action holds an identity and a count, and asks its domain
    # nothing else - so no map reaches it, at any rung.
    refuses common/triangular_decay_fixture.f90 graph_set_map
    refuses common/triangular_decay_fixture.f90 graph_set_representation
    # Provenance is Level 4's alone: no other level may assert an
    # embedding, and no fixture may.
    refuses level-5-field-calculus graph_inclusion_map
    refuses level-9-statement graph_inclusion_map
    # The label map is Level 4's alone, and only because carving
    # compels holding one - no other level carves, so no other level
    # may hold it.
    permits level-4-graph-calculus graph_label_map
    refuses level-0-carrier graph_label_map
    refuses level-5-field-calculus graph_label_map
    refuses level-9-statement graph_label_map
    permits common/triangular_decay_fixture.f90 graph_operation_view
    refuses common/triangular_decay_fixture.f90 class_graph_step
    refuses common/triangular_decay_fixture.f90 class_graph_marcher

    # An unclassified source still fails closed rather than silently
    # open - the ten levels are named, and nothing else is.
    allows level-10-nowhere fractal_graph
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
