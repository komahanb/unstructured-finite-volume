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
#     chain_carriers_fixture            earned at Level 0
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
# The fixture ladder is the tower's own stratification applied to
# itself. Level 0 declares carriers and may reach ONLY the carrier
# fixture; the relation fixture is one rung above it and is out of
# reach. A Level-0 source that says
#
#     use chain_relations_fixture
#
# is a layering violation and this gate must refuse it - which
# --selftest asserts directly on the allowlists, and the tower's
# history once got wrong.
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
        common/chain_carriers_fixture.f90) echo "graph_carrier" ;;
        common/chain_relations_fixture.f90) echo "graph_carrier graph_relation graph_binary_relation" ;;
        common/chain_algebra_fixture.f90) echo "graph_carrier graph_relation graph_relation_algebra graph_binary_relation" ;;
        common/shifted_laplacian_fixture.f90) echo "graph_carrier graph_grammar class_graph_field class_graph_differential_operator" ;;
        common/partitioned_shifted_laplacian_fixture.f90) echo "graph_carrier graph_grammar class_graph class_graph_field class_graph_partitioner class_graph_assembler shifted_laplacian_fixture" ;;
        common)            echo "__no_allowlist__" ;;

        # ---- L0: carriers only. NOTHING relational - not the
        #          relation nucleus, and not the Level-1 fixture.
        level-0-carrier)   echo "partitioned_pde_assert chain_carriers_fixture graph_carrier" ;;
        # ---- L1: + the relation nucleus
        level-1-relation)  echo "partitioned_pde_assert chain_carriers_fixture chain_relations_fixture graph_carrier graph_relation graph_binary_relation" ;;
        # ---- L2: + relation algebra
        level-2-relation-algebra) echo "partitioned_pde_assert chain_carriers_fixture chain_relations_fixture chain_algebra_fixture graph_carrier graph_relation graph_relation_algebra graph_binary_relation" ;;
        # ---- L3: + the ordinary graph realization. No partitioner.
        level-3-graph)     echo "partitioned_pde_assert chain_carriers_fixture chain_relations_fixture graph_carrier graph_relation graph_binary_relation class_graph" ;;
        # ---- L4: + the partitioner. Graph-to-graph only; no field.
        level-4-graph-calculus) echo "partitioned_pde_assert chain_carriers_fixture chain_relations_fixture chain_algebra_fixture graph_carrier graph_grammar graph_relation graph_relation_algebra graph_binary_relation class_graph class_graph_partitioner" ;;
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
# The gate's own test: the decision function, exercised on the
# question the tower once answered wrongly. This proves the
# ALLOWLISTS say what they must; the bare scan below proves the
# scanner acts on them.
#---------------------------------------------------------------------

if [ "$1" = "--selftest" ]; then
    fail=0

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
    permits level-0-carrier chain_carriers_fixture
    permits level-0-carrier partitioned_pde_assert
    permits level-0-carrier graph_carrier
    permits level-0-carrier iso_fortran_env
    refuses level-0-carrier chain_relations_fixture     # the C1 regression
    refuses level-0-carrier chain_algebra_fixture
    refuses level-0-carrier graph_relation
    refuses level-0-carrier graph_binary_relation
    refuses level-0-carrier class_graph

    # L1 stands on L0's fixture and adds its own; L2's is still above it.
    permits level-1-relation chain_carriers_fixture
    permits level-1-relation chain_relations_fixture
    refuses level-1-relation chain_algebra_fixture

    # The fixtures themselves are keyed per file, and the carrier
    # fixture is a Level-0 file: carriers only.
    permits common/chain_carriers_fixture.f90 graph_carrier
    refuses common/chain_carriers_fixture.f90 graph_binary_relation
    refuses common/partitioned_pde_assert.f90 graph_carrier

    # An unclassified source fails closed rather than silently open.
    allows common/not_a_real_fixture.f90 graph_carrier
    if [ "$?" -ne 2 ]; then
        echo " FAIL : an unclassified source did not fail closed"
        fail=1
    fi

    if [ "$fail" -ne 0 ]; then
        echo "IMPORT GATE: the layering decision is wrong"
        exit 1
    fi
    echo "import gate: L0 admits the carrier fixture and refuses the relation fixture"
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
