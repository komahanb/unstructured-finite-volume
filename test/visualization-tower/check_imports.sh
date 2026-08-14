#!/bin/bash
# The visualization tower's import gate, keyed PER LEVEL.
#
# Levels are the implementation architecture, so the dependency
# ceiling rises level by level - never gate by gate. A level source
# may `use` only the framework modules its own rung has been granted,
# and a directory with sources but no allowlist fails closed.
#
# common/ is NOT a hole in the stratification. Each shared fixture is
# keyed by FILE and classified by the first level that earns it:
#
#     visualization_assert             below everything (nothing)
#     visualization_carriers_fixture   earned at Level 0
#     visualization_relations_fixture  earned at Level 1
#     visualization_algebra_fixture    earned at Level 2
#     structural_renderer_fixture      earned at Level 4
#
# THE RENDERER IS THE ONE THAT MATTERS. It is earned at Level 4 and
# refused at 0-3, so no level below can quietly become a picture. A
# tower that could draw at Level 1 would never find out whether the
# relation algebra was necessary to draw at all.
#
#                  THE FRONTIER, HELD SHUT AT GATE A
#
# Two kinds of refusal live here, and the difference matters.
#
# STAGED refusals are ceilings that RISE, one rung at a time:
#
#     graph_relation, graph_binary_relation  earned at Level 1
#     graph_relation_algebra                 earned at Level 2
#     graph_structure                        earned at Level 3
#     structural_renderer_fixture            earned at Level 4
#
# UNIVERSAL refusals never lift at ANY level of this tower, and there
# are three families of them.
#
# 1. VALUES. graph_field_calculus and class_graph_field are Level 5's
#    business and Level 5 is unbuilt. Gate A's central claim is that
#    the structural picture precedes the numbers; a Gate-A level that
#    could reach a field could have quietly leaned on one.
#
# 2. OPERATORS. graph_grammar, class_graph_stencil, class_graph_step
#    and everything that solves, steps, linearizes or marches. The
#    types the brief names explicitly - discretization_operator and
#    stencil_operator - live inside class_graph_stencil,
#    class_graph_step, graph_calculus and graph_fitting, so refusing
#    those modules refuses the types; a second scan below refuses the
#    bare names too, in case a later level finds another road to them.
#    dependencies() is Level 6's question and Level 6 is unbuilt.
#
# 3. THE ORDINARY GRAPH. graph_profile and graph_algorithms, refused
#    at every level.
#
# THAT THIRD REFUSAL IS THE SHARPEST ASSERTION IN THIS GATE, and it
# is load-bearing evidence rather than hygiene.
#
# Gate A concludes that the relational nucleus already renders a
# rectangular typed dependency X -> Y, and that the ordinary-graph
# profile is intentionally too specialized to express one - both its
# readings demand a SINGLE vertex carrier. Because no level of this
# tower may name graph_profile, the pictures Level 4 produced cannot
# have leaned on it. The gate is what makes "no ordinary graph was
# required" checkable instead of promised.
#
# The derivative and adjoint fixtures stay refused for the older
# reason: this is not a derivative-family client, and its
# independence as evidence depends on that being mechanical.

here="$(cd "$(dirname "$0")" && pwd)"

intrinsics="iso_fortran_env iso_c_binding ieee_arithmetic ieee_exceptions ieee_features"

allowed_for() {
    case "$1" in
        # ---- shared fixtures, keyed by the level that earns them
        common/visualization_assert.f90) echo "" ;;
        common/visualization_carriers_fixture.f90) echo "graph_carrier" ;;
        common/visualization_relations_fixture.f90) echo "visualization_assert graph_carrier graph_relation graph_binary_relation" ;;
        common/visualization_algebra_fixture.f90) echo "graph_carrier graph_relation graph_binary_relation graph_relation_algebra" ;;
        common/structural_renderer_fixture.f90) echo "visualization_carriers_fixture graph_carrier graph_relation graph_binary_relation" ;;
        common)            echo "__no_allowlist__" ;;

        # ---- L0: carriers only. NOTHING relational - not the relation
        #          nucleus, and not the Level-1 fixture.
        level-0-carrier)   echo "visualization_assert visualization_carriers_fixture graph_carrier" ;;
        # ---- L1: + the relation nucleus and its binary specialization
        level-1-relation)  echo "visualization_assert visualization_carriers_fixture visualization_relations_fixture graph_carrier graph_relation graph_binary_relation" ;;
        # ---- L2: + relation algebra. Where D1, D2, D3 are DERIVED.
        level-2-relation-algebra) echo "visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture graph_carrier graph_relation graph_binary_relation graph_relation_algebra" ;;
        # ---- L3: + the relational graph container
        level-3-graph)     echo "visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture graph_carrier graph_relation graph_binary_relation graph_relation_algebra graph_structure" ;;
        # ---- L4: + the renderer, and NOTHING ELSE. No graph_profile:
        #          the level's conclusion is that the ordinary graph
        #          was not required, and the gate is what makes that
        #          claim mechanical.
        level-4-graph-calculus) echo "visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture structural_renderer_fixture graph_carrier graph_relation graph_binary_relation graph_relation_algebra graph_structure" ;;

        # ===== REVIEW GATE A ===== levels 5-9 are UNBUILT.

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

#---------------------------------------------------------------------
# THE NUMBERLESS LAW.
#
# Gate A's central claim is that the structural picture precedes the
# numerical coefficients entirely. Refusing class_graph_field and
# graph_field_calculus does not establish that: a level could simply
# have written
#
#     real(dp) :: w = 2.0_dp
#
# and helped itself. So the claim is checked directly - no real or
# complex declaration, and no literal carrying a decimal point or an
# exponent, in any Gate-A source.
#
# Comments are stripped first, because the tower's prose says the word
# "real" often and means the English one. The stripping cuts at the
# first '!', which is exact for these sources - none of them puts a
# bang inside a string.
#
# Integers are untouched. A carrier's size is 4, and that is not a
# coefficient.
#---------------------------------------------------------------------

holds_no_number() {
    local code
    code=$(sed 's/!.*//' "$1")
    echo "$code" | grep -qiE '\b(real|complex|double[[:space:]]+precision)\b' && return 1
    echo "$code" | grep -qE '[0-9]\.[0-9]|\.[0-9]+[eEdD]?|[0-9][eEdDqQ][+-]?[0-9]|_dp\b' && return 1
    return 0
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
# The gate's own test: the decision function, exercised on the three
# questions this tower must never get wrong - the fixture ladder, the
# rising ceiling, and the frontier. This proves the ALLOWLISTS say
# what they must; the bare scan below proves the scanner acts on them.
#---------------------------------------------------------------------

if [ "$1" = "--selftest" ]; then
    fail=0
    levels="level-0-carrier level-1-relation level-2-relation-algebra level-3-graph level-4-graph-calculus"
    before_one="level-0-carrier"
    before_two="$before_one level-1-relation"
    before_three="$before_two level-2-relation-algebra"
    before_four="$before_three level-3-graph"

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
    permits level-0-carrier visualization_carriers_fixture
    permits level-0-carrier visualization_assert
    permits level-0-carrier graph_carrier
    permits level-0-carrier iso_fortran_env
    refuses level-0-carrier visualization_relations_fixture   # the fixture ladder
    refuses level-0-carrier visualization_algebra_fixture
    refuses level-0-carrier graph_relation
    refuses level-0-carrier graph_binary_relation
    refuses level-0-carrier graph_structure

    # L1 stands on L0's fixture and adds its own; L2's is still above it.
    permits level-1-relation visualization_relations_fixture
    permits level-1-relation graph_binary_relation
    refuses level-1-relation visualization_algebra_fixture
    refuses level-1-relation graph_relation_algebra
    refuses level-1-relation graph_structure

    # L2 earns the algebra; the container is still one rung up.
    permits level-2-relation-algebra visualization_algebra_fixture
    permits level-2-relation-algebra graph_relation_algebra
    refuses level-2-relation-algebra graph_structure

    # L3 earns the container.
    permits level-3-graph graph_structure

    # ---- STAGED: the renderer, earned at Level 4 and nowhere below.
    #      A tower that could draw before it could derive would never
    #      find out whether the algebra was needed to draw at all.
    permits level-4-graph-calculus structural_renderer_fixture
    for lvl in $before_four; do
        refuses "$lvl" structural_renderer_fixture
    done
    for lvl in $before_three; do
        refuses "$lvl" graph_structure
    done
    for lvl in $before_two; do
        refuses "$lvl" graph_relation_algebra
        refuses "$lvl" visualization_algebra_fixture
    done
    for lvl in $before_one; do
        refuses "$lvl" graph_relation
        refuses "$lvl" graph_binary_relation
        refuses "$lvl" visualization_relations_fixture
    done

    # ---- UNIVERSAL 1: VALUES. Gate A holds no number at all, and the
    #      levels that would supply one are unbuilt.
    for lvl in $levels; do
        refuses "$lvl" graph_field_calculus
        refuses "$lvl" class_graph_field
    done

    # ---- UNIVERSAL 2: OPERATORS. Nothing that steps, solves,
    #      linearizes, marches, or owns dependencies().
    for lvl in $levels; do
        refuses "$lvl" graph_grammar
        refuses "$lvl" class_graph_stencil
        refuses "$lvl" class_graph_step
        refuses "$lvl" graph_calculus
        refuses "$lvl" graph_fitting
        refuses "$lvl" class_graph_linearization
        refuses "$lvl" graph_minimization
        refuses "$lvl" class_graph_gmres
        refuses "$lvl" class_graph_newton
        refuses "$lvl" class_graph_marcher
        refuses "$lvl" derivative_fixture
        refuses "$lvl" adjoint_fixture
    done

    # ---- UNIVERSAL 3: THE ORDINARY GRAPH. The load-bearing one.
    #      Level 4 concludes that a rectangular typed dependency needs
    #      no ordinary-graph reading; because no level may name the
    #      profile, that conclusion cannot have been reached with one.
    for lvl in $levels; do
        refuses "$lvl" graph_profile
        refuses "$lvl" graph_algorithms
        refuses "$lvl" class_graph
        refuses "$lvl" class_stored_graph
    done

    # The illegal upward import the brief names, spelled out.
    refuses level-2-relation-algebra graph_field_calculus

    # The fixtures themselves are keyed per file.
    permits common/visualization_carriers_fixture.f90 graph_carrier
    refuses common/visualization_carriers_fixture.f90 graph_binary_relation
    refuses common/visualization_carriers_fixture.f90 visualization_assert
    refuses common/visualization_assert.f90 graph_carrier
    permits common/structural_renderer_fixture.f90 visualization_carriers_fixture
    refuses common/structural_renderer_fixture.f90 graph_structure
    refuses common/structural_renderer_fixture.f90 graph_profile
    refuses common/visualization_relations_fixture.f90 graph_relation_algebra
    refuses common/visualization_algebra_fixture.f90 graph_structure

    # ---- THE NUMBERLESS LAW, tested on the checker itself.
    numbered=$(mktemp);   printf 'module m\n real(dp) :: w = 2.0_dp\nend module\n' > "$numbered"
    numberless=$(mktemp); printf '! a real picture, honestly\nmodule m\n integer, parameter :: N = 4\n x = reshape([1, 2], [2, 1])\n if (k .eq. 3) call go(a)\nend module\n' > "$numberless"
    exponent=$(mktemp);   printf 'module m\n t = 1.0e-3\nend module\n' > "$exponent"

    if holds_no_number "$numbered"; then
        echo " FAIL : the numberless law accepted a real coefficient"
        fail=1
    fi
    if holds_no_number "$exponent"; then
        echo " FAIL : the numberless law accepted an exponent literal"
        fail=1
    fi
    if holds_no_number "$numberless"; then :; else
        echo " FAIL : the numberless law refused integers and the word 'real' in a comment"
        fail=1
    fi
    rm -f "$numbered" "$numberless" "$exponent"

    # An unclassified source still fails closed rather than silently
    # open - the five built levels are named, and nothing else is.
    allows level-5-field-calculus graph_carrier
    if [ "$?" -ne 2 ]; then
        echo " FAIL : an unbuilt level did not fail closed"
        fail=1
    fi
    allows level-9-statement graph_carrier
    if [ "$?" -ne 2 ]; then
        echo " FAIL : an unbuilt level did not fail closed"
        fail=1
    fi

    if [ "$fail" -ne 0 ]; then
        echo "IMPORT GATE: the layering decision is wrong"
        exit 1
    fi
    echo "import gate: the ladder rises one rung at a time, and the ordinary graph is refused at every level"
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

        # Gate A holds no number, and this is where that is enforced
        # rather than asserted.
        if ! holds_no_number "$src"; then
            echo "IMPORT GATE: $(basename "$src") in $name carries a real declaration or a non-integer literal - Gate A holds no coefficients"
            violation=1
        fi

        # The two production TYPE names the brief forbids by name.
        # Refusing their modules already refuses them; this catches a
        # road to them the module scan would not see.
        for banned in discretization_operator stencil_operator; do
            if grep -qiE "(^|[^a-zA-Z0-9_])$banned([^a-zA-Z0-9_]|$)" "$src"; then
                echo "IMPORT GATE: $(basename "$src") in $name names '$banned' - that is Level 6's question"
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
