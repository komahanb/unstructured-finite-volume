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
        #
        # 2026-08-16: graph_carrier is retired as a domain source, and
        # the one grant it used to carry splits into four. Each source
        # gets ONLY the part it uses, and the four rules are these:
        #
        #   fractal_graph            identity - granted broadly, because
        #                            WHICH set is a question everything
        #                            asks
        #   graph_set_representation granted only where a representation
        #                            is CONSTRUCTED
        #   graph_set_map            granted only where one is BOUND or
        #                            QUERIED
        #   graph_inclusion_map      granted only where provenance is
        #                            ASSERTED, or a carve is CALLED
        #   graph_label_map          granted only where a label is
        #                            QUERIED, or a carve is CALLED
        #
        # The last two say "or a carve is called" because carving binds
        # extension, label and embedding TOGETHER so no half-described
        # set escapes. A caller of a carving operation must therefore
        # HOLD both maps even where it reads neither. Nothing in this
        # tower ever asks a domain its name; every label grant below is
        # held-not-queried, and says so.
        common/partitioned_pde_assert.f90) echo "" ;;
        # Constructs the three chain carriers: the only source here
        # that builds a representation.
        common/chain_carriers_fixture.f90) echo "fractal_graph graph_set_representation graph_set_map" ;;
        common/chain_relations_fixture.f90) echo "fractal_graph graph_set_map graph_relation graph_binary_relation" ;;
        common/chain_algebra_fixture.f90) echo "fractal_graph graph_set_map graph_relation graph_relation_algebra graph_binary_relation" ;;
        # The operation holds an identity and a count. No map.
        common/shifted_laplacian_fixture.f90) echo "fractal_graph graph_operation_view graph_ordinary_view graph_field_calculus class_graph_field class_graph_differential_operator" ;;
        # Calls partition_data and assemble_data - both carve - and
        # queries the domains it transports. It also DESCRIBES the
        # part carriers it transports onto, so it constructs
        # representations and earns the module that makes them.
        common/partitioned_shifted_laplacian_fixture.f90) echo "fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_operation_view graph_ordinary_view graph_field_calculus class_graph class_graph_field class_graph_partitioner class_graph_assembler shifted_laplacian_fixture" ;;
        common)            echo "__no_allowlist__" ;;

        # ---- L0: sets only. NOTHING relational.
        level-0-carrier)   echo "partitioned_pde_assert chain_carriers_fixture fractal_graph graph_set_map" ;;
        # ---- L1: + the relation nucleus
        level-1-relation)  echo "partitioned_pde_assert chain_carriers_fixture chain_relations_fixture fractal_graph graph_set_map graph_relation graph_binary_relation" ;;
        # ---- L2: + relation algebra
        level-2-relation-algebra) echo "partitioned_pde_assert chain_carriers_fixture chain_relations_fixture chain_algebra_fixture fractal_graph graph_set_map graph_relation graph_relation_algebra graph_binary_relation" ;;
        # ---- L3: + the ordinary graph realization. No partitioner.
        #          A level that REALIZES a stored_graph must describe
        #          that graph's own carriers - they are its declared
        #          domains, distinct from the oracle's - and describing
        #          means constructing a representation. So every level
        #          from here on earns graph_set_representation by
        #          building one, not by inheriting a permission.
        level-3-graph)     echo "partitioned_pde_assert chain_carriers_fixture chain_relations_fixture fractal_graph graph_set_representation graph_set_map graph_relation graph_binary_relation class_graph" ;;
        # ---- L4: + the partitioner. Graph-to-graph only; no field.
        #          Calls owned_vertices, which CARVES - hence both the
        #          inclusion map and a label map it never reads.
        level-4-graph-calculus) echo "partitioned_pde_assert chain_carriers_fixture chain_relations_fixture chain_algebra_fixture fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_ordinary_view graph_relation graph_relation_algebra graph_binary_relation class_graph class_graph_partitioner" ;;
        # ---- L5: + fields and transport. NO differential operator.
        #          The one level that both CARVES a probe subset of its
        #          own and ASSERTS where a transported domain came from.
        level-5-field-calculus) echo "partitioned_pde_assert fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_ordinary_view graph_field_calculus class_graph class_graph_field class_graph_partitioner class_graph_assembler" ;;
        # ---- L6: + the differential operator and its fixture.
        #          Still no solver. Transport carves; the refusal
        #          builds a foreign carrier, and is keyed below.
        level-6-discretization) echo "partitioned_pde_assert shifted_laplacian_fixture fractal_graph graph_set_representation graph_set_map graph_label_map graph_inclusion_map graph_ordinary_view graph_field_calculus class_graph class_graph_field class_graph_partitioner class_graph_assembler class_graph_differential_operator" ;;
        level-6-discretization/refusal.f90) echo "partitioned_pde_assert shifted_laplacian_fixture fractal_graph graph_set_representation graph_set_map graph_field_calculus class_graph class_graph_field class_graph_differential_operator" ;;
        # ---- L7: + minimization. Consumes Level 6, not Level 8.
        level-7-minimization) echo "partitioned_pde_assert shifted_laplacian_fixture fractal_graph graph_set_representation graph_set_map graph_ordinary_view graph_field_calculus class_graph class_graph_field class_graph_gmres" ;;
        # ---- L8: + the partitioned constitution, which consumes L6.
        #          It asks its domains WHICH and nothing else, so it is
        #          granted identity alone - no map of any kind.
        level-8-constitution) echo "partitioned_pde_assert shifted_laplacian_fixture partitioned_shifted_laplacian_fixture fractal_graph graph_ordinary_view graph_field_calculus class_graph class_graph_field class_graph_gmres" ;;
        # ---- L9: the statement, consuming Level 8 and the solver.
        level-9-statement) echo "partitioned_pde_assert shifted_laplacian_fixture partitioned_shifted_laplacian_fixture fractal_graph graph_set_representation graph_set_map graph_ordinary_view graph_field_calculus class_graph class_graph_field class_graph_gmres" ;;

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
    permits level-0-carrier fractal_graph
    permits level-0-carrier graph_set_map
    refuses level-0-carrier graph_carrier            # the retired source
    # A level takes its domains from the fixture, so it builds no
    # representation and is not granted the module that makes one.
    refuses level-0-carrier graph_set_representation
    # It neither carves nor asserts an embedding.
    refuses level-0-carrier graph_inclusion_map
    refuses level-0-carrier graph_label_map
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
    permits common/chain_carriers_fixture.f90 fractal_graph
    permits common/chain_carriers_fixture.f90 graph_set_representation
    permits common/chain_carriers_fixture.f90 graph_set_map
    refuses common/chain_carriers_fixture.f90 graph_inclusion_map
    refuses common/chain_carriers_fixture.f90 graph_label_map
    refuses common/chain_carriers_fixture.f90 graph_binary_relation
    refuses common/partitioned_pde_assert.f90 fractal_graph

    # THE OPERATION HOLDS IDENTITY AND A COUNT, AND NOTHING ELSE.
    refuses common/shifted_laplacian_fixture.f90 graph_set_map
    refuses common/shifted_laplacian_fixture.f90 graph_set_representation
    refuses common/shifted_laplacian_fixture.f90 graph_inclusion_map
    refuses common/shifted_laplacian_fixture.f90 graph_label_map

    # Carving compels holding both provenance maps, so the sources
    # that carve get them and the sources that do not are refused.
    permits common/partitioned_shifted_laplacian_fixture.f90 graph_set_representation
    permits common/partitioned_shifted_laplacian_fixture.f90 graph_inclusion_map
    permits common/partitioned_shifted_laplacian_fixture.f90 graph_label_map
    permits level-4-graph-calculus graph_inclusion_map
    permits level-5-field-calculus graph_set_representation
    permits level-3-graph graph_set_representation
    # L0-L2 take their domains from the fixture and build none.
    refuses level-1-relation graph_set_representation
    refuses level-2-relation-algebra graph_set_representation
    refuses level-7-minimization graph_inclusion_map
    refuses level-7-minimization graph_label_map
    refuses level-9-statement graph_inclusion_map
    # Level 8 asks only WHICH, so it is granted identity alone.
    refuses level-8-constitution graph_set_map

    # An unclassified source fails closed rather than silently open.
    allows common/not_a_real_fixture.f90 fractal_graph
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
