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
#     visualization_values_fixture     earned at Level 5
#     valued_renderer_fixture          earned at Level 5
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
#     relation_finitary, relation_binary  earned at Level 1
#     relation_algebra                 earned at Level 2
#     graph_fractal, view_relational   earned at Level 3
#     structural_renderer_fixture            earned at Level 4
#     field_stored, field_calculus, and the two
#     coefficient fixtures                   earned at Level 5
#
# UNIVERSAL refusals never lift at any level BUILT SO FAR, and there
# are two families of them - values used to be a third, and Level 5 is
# where it stopped being one.
#
# 0. VALUES, refused at L0-L4 and EARNED AT L5. field_calculus
#    and field_stored are Level 5's business. Gate A's claim is
#    that the structural picture precedes the numbers, so the four
#    levels below Level 5 must remain unable to reach a coefficient -
#    which is a ceiling that lifted, not a refusal that was deleted.
#
# 2. OPERATORS. view_directed, operation_stencil, operation_step
#    and everything that solves, steps, linearizes or marches. The
#    types the brief names explicitly - discretization_operator and
#    stencil_operator - live inside operation_stencil,
#    operation_step, graph_calculus and operation_fitting, so refusing
#    those modules refuses the types; a second scan below refuses the
#    bare names too, in case a later level finds another road to them.
#    dependencies() is Level 6's question and Level 6 is unbuilt.
#
# 3. THE ORDINARY GRAPH. and relation_algorithms, refused
#    at every level, Level 5 included.
#
# THAT THIRD REFUSAL IS THE SHARPEST ASSERTION IN THIS GATE, and it
# is load-bearing evidence rather than hygiene.
#
# Gate A concludes that the relational nucleus already renders a
# rectangular typed dependency X -> Y, and that the ordinary-graph
# profile is intentionally too specialized to express one - both its
# readings demand a SINGLE vertex carrier. Because no level of this
# tower may name, the pictures Level 4 produced cannot
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
        # 2026-08-16: graph_carrier is retired as a domain source. The
        # one grant splits into five capabilities, and each source gets
        # only what it uses:
        #
        #   graph_fractal            identity - broadly
        #   map_set_representation only where one is CONSTRUCTED
        #   map_set            only where one is BOUND or QUERIED
        #   map_inclusion      provenance ASSERTED, or a CARVE
        #   map_label          a label QUERIED, or a CARVE
        #
        # Carving implies all three maps because it binds extension,
        # label and embedding atomically, so no half-described set can
        # escape. A caller of a carving operation therefore HOLDS the
        # label map even where it reads no label.
        # 2026-08-16: graph_carrier is retired as a domain source. The
        # one grant splits into five capabilities, and each source gets
        # only what it uses:
        #
        #   graph_fractal            identity - broadly
        #   map_set_representation only where one is CONSTRUCTED
        #   map_set            only where one is BOUND or QUERIED
        #   map_inclusion      provenance ASSERTED, or a CARVE
        #   map_label          a label QUERIED, or a CARVE
        #
        # Carving implies all three maps because it binds extension,
        # label and embedding atomically, so no half-described set can
        # escape. A caller of a carving operation therefore HOLDS the
        # label map even where it reads no label.
        common/visualization_assert.f90) echo "" ;;
        common/visualization_carriers_fixture.f90) echo "graph_fractal map_set_representation map_set map_label" ;;
        common/visualization_relations_fixture.f90) echo "graph_fractal map_set visualization_assert relation_finitary relation_binary" ;;
        common/visualization_algebra_fixture.f90) echo "graph_fractal map_set  relation_finitary relation_binary relation_algebra" ;;
        common/structural_renderer_fixture.f90) echo "graph_fractal map_set map_label visualization_carriers_fixture relation_finitary relation_binary" ;;
        common/visualization_values_fixture.f90) echo "graph_fractal map_set  field_stored" ;;
        common/valued_renderer_fixture.f90) echo "graph_fractal map_set map_label visualization_carriers_fixture structural_renderer_fixture relation_finitary relation_binary field_calculus field_stored" ;;
        common/production_discretization_fixture.f90) echo "graph_fractal operation_action operation_stencil operation_step" ;;
        common/production_pattern_renderer_fixture.f90) echo "graph_fractal map_set map_label visualization_carriers_fixture structural_renderer_fixture relation_finitary view_directed" ;;
        common)            echo "__no_allowlist__" ;;

        # ---- L0: carriers only. NOTHING relational - not the relation
        #          nucleus, and not the Level-1 fixture.
        level-0-carrier)   echo "graph_fractal map_set_representation map_set map_label visualization_assert visualization_carriers_fixture" ;;
        # ---- L1: + the relation nucleus and its binary specialization
        level-1-relation)  echo "graph_fractal map_set map_label visualization_assert visualization_carriers_fixture visualization_relations_fixture relation_finitary relation_binary" ;;
        # ---- L2: + relation algebra. Where D1, D2, D3 are DERIVED.
        level-2-relation-algebra) echo "graph_fractal map_set map_label visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture relation_finitary relation_binary relation_algebra" ;;
        # ---- L3: + the relational graph container
        level-3-graph)     echo "graph_fractal map_set_representation map_set map_label map_inclusion visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture relation_finitary relation_binary relation_algebra view_relational" ;;
        # ---- L4: + the renderer, and NOTHING ELSE. No:
        #          the level's conclusion is that the ordinary graph
        #          was not required, and the gate is what makes that
        #          claim mechanical.
        level-4-graph-calculus) echo "graph_fractal map_set_representation map_set map_label map_inclusion visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture structural_renderer_fixture relation_finitary relation_binary relation_algebra view_relational" ;;

        # ===== REVIEW GATE A =====

        # ---- L5: + the field nucleus, and the two fixtures that carry
        #          coefficients. The ordinary graph and every operator
        #          module stay refused: values arrived, machinery did
        #          not.
        level-5-field-calculus) echo "graph_fractal map_set map_label visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture structural_renderer_fixture visualization_values_fixture valued_renderer_fixture relation_finitary relation_binary relation_algebra field_calculus field_stored" ;;

        # ---- L6: + PRODUCTION DISCRETIZATION, and only the three
        #          modules the level actually names. graph_calculus and
        #          view_directed_stored are NOT permitted: the test reaches the
        #          concrete citizens, never the abstract type or the
        #          stored graph directly, and a ceiling permits what is
        #          used rather than what is nearby. operation_fitting is
        #          refused outright, and so is everything that solves.
        level-6-discretization) echo "graph_fractal map_set_representation map_set map_label visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture structural_renderer_fixture production_discretization_fixture production_pattern_renderer_fixture relation_finitary relation_binary relation_algebra view_directed operation_stencil operation_step" ;;

        # ---- L7: + MINIMIZATION, and only the concrete the
        #          experiment uses. operation_minimization itself is not
        #          named (jacobi inherits attach/matvec/sweep_order/
        #          diagonal), and gauss_seidel is inspected in the
        #          census rather than run. Nothing that solves beyond
        #          jacobi: no gmres, no newton, no multigrid, no
        #          marcher, no linearization.
        level-7-minimization) echo "graph_fractal map_set_representation map_set map_label visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture structural_renderer_fixture production_pattern_renderer_fixture relation_finitary relation_binary relation_algebra view_directed view_directed_stored operation_stencil operation_jacobi" ;;

        # ---- L8: the constitution. Both axes and the context at once,
        #          on the vocabulary Levels 6 and 7 already earned.
        level-8-constitution) echo "graph_fractal map_set_representation map_set map_label visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture structural_renderer_fixture production_pattern_renderer_fixture relation_finitary relation_binary relation_algebra view_directed view_directed_stored operation_stencil operation_step operation_jacobi" ;;
        # ---- L9: the statement, on the same constitution plus the
        #          coefficient view Level 5 earned.
        level-9-statement) echo "graph_fractal map_set_representation map_set map_label visualization_assert visualization_carriers_fixture visualization_relations_fixture visualization_algebra_fixture structural_renderer_fixture valued_renderer_fixture visualization_values_fixture production_pattern_renderer_fixture relation_finitary relation_binary relation_algebra view_directed field_calculus view_directed_stored field_stored operation_stencil operation_step operation_jacobi" ;;

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
# THE NUMBERLESS LAW, AND WHERE IT STOPS.
#
# Gate A's central claim is that the structural picture precedes the
# numerical coefficients entirely. Refusing field_stored and
# field_calculus does not establish that: a level could simply
# have written
#
#     real(dp) :: w = 2.0_dp
#
# and helped itself. So the claim is checked directly - no real or
# complex declaration, and no literal carrying a decimal point or an
# exponent, in any source below Level 5.
#
# LEVEL 5 IS WHERE THE NUMBERS ARRIVE, AND THEY ARRIVE DELIBERATELY.
# The law is therefore a CEILING that lifts at exactly one rung, not a
# rule that was deleted when it became inconvenient:
#
#     L0-L4 and their fixtures     no coefficient may be written
#     L5 and its two fixtures      coefficients are the subject
#
# So the Gate-A claim stays mechanical after Gate A, which is the only
# way it stays worth anything: Level 4's renderer must still be unable
# to hold a number, or "Level 4 renders mathematics that has none"
# would become a promise again.
#
# Comments are stripped first, because the tower's prose says the word
# "real" often and means the English one. The stripping cuts at the
# first '!', which is exact for these sources - none of them puts a
# bang inside a string.
#
# Integers are untouched everywhere. A carrier's size is 4, and that
# is not a coefficient.
#---------------------------------------------------------------------

# numbers_allowed <level[/file]>
#     0  this source is entitled to carry coefficients
#     1  it is not
numbers_allowed() {
    case "$1" in
        level-5-field-calculus|level-5-field-calculus/*) return 0 ;;
        common/visualization_values_fixture.f90)         return 0 ;;
        common/valued_renderer_fixture.f90)              return 0 ;;
        level-6-discretization|level-6-discretization/*) return 0 ;;
        level-7-minimization|level-7-minimization/*)     return 0 ;;
        level-8-constitution|level-8-constitution/*)     return 0 ;;
        level-9-statement|level-9-statement/*)           return 0 ;;
        common/production_discretization_fixture.f90)    return 0 ;;
        *)                                               return 1 ;;
    esac
}

# discretization_allowed <level[/file]>
#     0  this source may name production discretization vocabulary
#     1  it may not
discretization_allowed() {
    case "$1" in
        level-6-discretization|level-6-discretization/*)  return 0 ;;
        level-7-minimization|level-7-minimization/*)      return 0 ;;
        level-8-constitution|level-8-constitution/*)      return 0 ;;
        level-9-statement|level-9-statement/*)            return 0 ;;
        common/production_discretization_fixture.f90)     return 0 ;;
        common/production_pattern_renderer_fixture.f90)   return 0 ;;
        *)                                                return 1 ;;
    esac
}

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
    with_five="$levels level-5-field-calculus"
    before_one="level-0-carrier"
    before_two="$before_one level-1-relation"
    before_three="$before_two level-2-relation-algebra"
    before_four="$before_three level-3-graph"

    permits() {
        if allows "$1" "$2"; then :; else
            echo "FAIL : the import gate refused '$2' at $1"
            fail=1
        fi
    }
    refuses() {
        if allows "$1" "$2"; then
            echo "FAIL : the import gate permitted '$2' at $1"
            fail=1
        fi
    }

    # L0 earns carriers and the carrier fixture, and NOTHING relational.
    permits level-0-carrier visualization_carriers_fixture
    permits level-0-carrier visualization_assert
    permits level-0-carrier graph_fractal
    permits level-0-carrier iso_fortran_env
    refuses level-0-carrier visualization_relations_fixture   # the fixture ladder
    refuses level-0-carrier visualization_algebra_fixture
    refuses level-0-carrier relation_finitary
    refuses level-0-carrier relation_binary
    refuses level-0-carrier view_relational

    # L1 stands on L0's fixture and adds its own; L2's is still above it.
    permits level-1-relation visualization_relations_fixture
    permits level-1-relation relation_binary
    refuses level-1-relation visualization_algebra_fixture
    refuses level-1-relation relation_algebra
    refuses level-1-relation view_relational

    # L2 earns the algebra; the container is still one rung up.
    permits level-2-relation-algebra visualization_algebra_fixture
    permits level-2-relation-algebra relation_algebra
    refuses level-2-relation-algebra view_relational

    # L3 earns the container.
    permits level-3-graph graph_fractal
    permits level-3-graph view_relational

    # ---- STAGED: the renderer, earned at Level 4 and nowhere below.
    #      A tower that could draw before it could derive would never
    #      find out whether the algebra was needed to draw at all.
    permits level-4-graph-calculus structural_renderer_fixture
    for lvl in $before_four; do
        refuses "$lvl" structural_renderer_fixture
    done
    for lvl in $before_three; do
        refuses "$lvl" view_relational
    done
    for lvl in $before_two; do
        refuses "$lvl" relation_algebra
        refuses "$lvl" visualization_algebra_fixture
    done
    for lvl in $before_one; do
        refuses "$lvl" relation_finitary
        refuses "$lvl" relation_binary
        refuses "$lvl" visualization_relations_fixture
    done

    # ---- UNIVERSAL 1: VALUES. Gate A holds no number at all, and the
    #      levels that would supply one are unbuilt.
    for lvl in $levels; do
        refuses "$lvl" field_calculus
        refuses "$lvl" field_stored
    done

    # ---- UNIVERSAL 2: OPERATORS. Nothing that steps, solves,
    #      linearizes, marches, or owns dependencies().
    for lvl in $with_five; do
        refuses "$lvl" operation_action
        refuses "$lvl" view_directed
        refuses "$lvl" operation_stencil
        refuses "$lvl" operation_step
        refuses "$lvl" graph_calculus
        refuses "$lvl" operation_fitting
        refuses "$lvl" operation_linearization
        refuses "$lvl" operation_minimization
        refuses "$lvl" operation_gmres
        refuses "$lvl" operation_newton
        refuses "$lvl" operation_marching
        refuses "$lvl" derivative_fixture
        refuses "$lvl" adjoint_fixture
    done

    # ---- UNIVERSAL 3: THE ORDINARY GRAPH. The load-bearing one.
    #      Level 4 concludes that a rectangular typed dependency needs
    #      no ordinary-graph reading; because no level may name the
    #      profile, that conclusion cannot have been reached with one.
    for lvl in $with_five; do
        refuses "$lvl" graph_profile
        refuses "$lvl" relation_algorithms
        refuses "$lvl" view_directed_stored
        refuses "$lvl" class_stored_graph
    done

    # The illegal upward import the brief names, spelled out.
    refuses level-2-relation-algebra field_calculus

    # The fixtures themselves are keyed per file.
    permits common/visualization_carriers_fixture.f90 graph_fractal
    refuses common/visualization_carriers_fixture.f90 relation_binary
    refuses common/visualization_carriers_fixture.f90 visualization_assert
    refuses common/visualization_assert.f90 graph_fractal
    permits common/structural_renderer_fixture.f90 visualization_carriers_fixture
    refuses common/structural_renderer_fixture.f90 view_relational
    refuses common/structural_renderer_fixture.f90 graph_profile
    refuses common/visualization_relations_fixture.f90 relation_algebra
    refuses common/visualization_algebra_fixture.f90 view_relational

    # ---- L5 earns the field nucleus and its two coefficient
    #      fixtures. Nothing below it may reach any of them.
    permits level-5-field-calculus field_stored
    permits level-5-field-calculus field_calculus
    permits level-5-field-calculus visualization_values_fixture
    permits level-5-field-calculus valued_renderer_fixture
    permits level-5-field-calculus structural_renderer_fixture
    for lvl in $levels; do
        refuses "$lvl" visualization_values_fixture
        refuses "$lvl" valued_renderer_fixture
    done

    # L5 does NOT lift the frontier. Values arrived; machinery did not.
    refuses level-5-field-calculus graph_profile
    refuses level-5-field-calculus relation_algorithms
    refuses level-5-field-calculus operation_step
    refuses level-5-field-calculus operation_stencil
    refuses level-5-field-calculus operation_action
    refuses level-5-field-calculus view_directed
    refuses level-5-field-calculus operation_linearization

    # The valued renderer stands ON the structural one, never inside
    # it - and the structural one may not learn about fields.
    permits common/valued_renderer_fixture.f90 structural_renderer_fixture
    permits common/valued_renderer_fixture.f90 field_stored
    refuses common/structural_renderer_fixture.f90 field_stored
    refuses common/structural_renderer_fixture.f90 field_calculus
    refuses common/structural_renderer_fixture.f90 valued_renderer_fixture
    refuses common/visualization_values_fixture.f90 structural_renderer_fixture

    # ---- L6 earns production discretization, and only what it uses.
    permits level-6-discretization operation_stencil
    permits level-6-discretization operation_step
    permits level-6-discretization view_directed
    permits level-6-discretization production_discretization_fixture
    permits level-6-discretization production_pattern_renderer_fixture
    permits level-6-discretization structural_renderer_fixture
    for lvl in $with_five; do
        refuses "$lvl" operation_stencil
        refuses "$lvl" operation_step
        refuses "$lvl" operation_action
        refuses "$lvl" view_directed
        refuses "$lvl" graph_calculus
        refuses "$lvl" production_discretization_fixture
        refuses "$lvl" production_pattern_renderer_fixture
    done

    # L6 does NOT lift the rest of the frontier. Nothing that solves,
    # and nothing merely adjacent to discretization machinery.
    refuses level-6-discretization operation_fitting
    refuses level-6-discretization operation_minimization
    refuses level-6-discretization operation_gmres
    refuses level-6-discretization operation_newton
    refuses level-6-discretization operation_marching
    refuses level-6-discretization operation_linearization
    refuses level-6-discretization graph_profile
    refuses level-6-discretization relation_algorithms
    refuses level-6-discretization graph_calculus
    refuses level-6-discretization view_directed_stored
    refuses level-6-discretization field_stored
    refuses level-6-discretization valued_renderer_fixture

    # The production fixtures stand where they are earned, and the
    # Level-4 renderer still knows nothing of production.
    permits common/production_pattern_renderer_fixture.f90 view_directed
    permits common/production_pattern_renderer_fixture.f90 structural_renderer_fixture
    permits common/production_discretization_fixture.f90 operation_stencil
    refuses common/structural_renderer_fixture.f90 operation_action
    refuses common/structural_renderer_fixture.f90 view_directed
    refuses common/structural_renderer_fixture.f90 operation_stencil
    refuses common/valued_renderer_fixture.f90 operation_step
    refuses common/production_pattern_renderer_fixture.f90 operation_stencil

    # ---- L7 earns minimization, and only the concrete it uses.
    permits level-7-minimization operation_jacobi
    permits level-7-minimization view_directed_stored
    permits level-7-minimization operation_stencil
    permits level-7-minimization production_pattern_renderer_fixture
    for lvl in $with_five level-6-discretization; do
        refuses "$lvl" operation_jacobi
        refuses "$lvl" operation_gauss_seidel
        refuses "$lvl" operation_minimization
    done

    # L7 does NOT lift the rest of the frontier.
    refuses level-7-minimization operation_gmres
    refuses level-7-minimization operation_newton
    refuses level-7-minimization operation_multigrid
    refuses level-7-minimization operation_marching
    refuses level-7-minimization operation_linearization
    refuses level-7-minimization operation_fitting
    refuses level-7-minimization graph_profile
    refuses level-7-minimization operation_gauss_seidel
    refuses level-7-minimization operation_minimization

    # ---- THE DISCRETIZATION VOCABULARY: WHERE IT LIFTS.
    for below in level-0-carrier/test.f90 level-4-graph-calculus/test.f90 \
                 level-5-field-calculus/test.f90 \
                 common/structural_renderer_fixture.f90 \
                 common/valued_renderer_fixture.f90; do
        if discretization_allowed "$below"; then
            echo "FAIL : production discretization vocabulary allowed at $below"
            fail=1
        fi
    done
    for at_six in level-6-discretization/test.f90 \
                  level-7-minimization/test.f90 \
                  common/production_discretization_fixture.f90; do
        if discretization_allowed "$at_six"; then :; else
            echo "FAIL : Level 6 was refused its own production vocabulary at $at_six"
            fail=1
        fi
    done

    # ---- THE NUMBERLESS LAW: WHERE IT HOLDS AND WHERE IT LIFTS.
    if numbers_allowed level-4-graph-calculus/test.f90; then
        echo "FAIL : the numberless law lifted below Level 5"
        fail=1
    fi
    for below in level-0-carrier/test.f90 level-1-relation/test.f90 \
                 level-2-relation-algebra/test.f90 level-3-graph/test.f90 \
                 common/visualization_assert.f90 \
                 common/structural_renderer_fixture.f90; do
        if numbers_allowed "$below"; then
            echo "FAIL : the numberless law lifted at $below"
            fail=1
        fi
    done
    for at_five in level-5-field-calculus/test.f90 \
                   common/visualization_values_fixture.f90 \
                   common/valued_renderer_fixture.f90 \
                   level-6-discretization/test.f90 \
                   level-7-minimization/test.f90 \
                   common/production_discretization_fixture.f90; do
        if numbers_allowed "$at_five"; then :; else
            echo "FAIL : Level 5 was refused its own coefficients at $at_five"
            fail=1
        fi
    done

    # ---- THE NUMBERLESS LAW, tested on the checker itself.
    numbered=$(mktemp);   printf 'module m\n real(dp) :: w = 2.0_dp\nend module\n' > "$numbered"
    numberless=$(mktemp); printf '! a real picture, honestly\nmodule m\n integer, parameter :: N = 4\n x = reshape([1, 2], [2, 1])\n if (k .eq. 3) call go(a)\nend module\n' > "$numberless"
    exponent=$(mktemp);   printf 'module m\n t = 1.0e-3\nend module\n' > "$exponent"

    if holds_no_number "$numbered"; then
        echo "FAIL : the numberless law accepted a real coefficient"
        fail=1
    fi
    if holds_no_number "$exponent"; then
        echo "FAIL : the numberless law accepted an exponent literal"
        fail=1
    fi
    if holds_no_number "$numberless"; then :; else
        echo "FAIL : the numberless law refused integers and the word 'real' in a comment"
        fail=1
    fi
    rm -f "$numbered" "$numberless" "$exponent"

    # An unclassified source still fails closed rather than silently
    # open - the five built levels are named, and nothing else is.
    allows level-10-nowhere graph_fractal
    if [ "$?" -ne 2 ]; then
        echo "FAIL : an unbuilt level did not fail closed"
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

        # Below Level 5 there are no numbers, and this is where that
        # is enforced rather than asserted. At Level 5 there are, and
        # the ceiling lifts for exactly those three sources.
        if ! numbers_allowed "$key"; then
            if ! holds_no_number "$src"; then
                echo "IMPORT GATE: $(basename "$src") in $name carries a real declaration or a non-integer literal - only Level 5 holds coefficients"
                violation=1
            fi
        fi

        # THE PRODUCTION DISCRETIZATION VOCABULARY, level-sensitive.
        #
        # Refusing the modules already refuses the types; this catches
        # a road to them the module scan would not see. Like the
        # numberless law, it is a CEILING THAT LIFTS rather than a rule
        # deleted when Level 6 needed it: forbidden L0-L5, allowed L6.
        if ! discretization_allowed "$key"; then
            for banned in discretization_operator stencil_operator step_operator; do
                if grep -qiE "(^|[^a-zA-Z0-9_])$banned([^a-zA-Z0-9_]|$)" "$src"; then
                    echo "IMPORT GATE: $(basename "$src") in $name names '$banned' - production discretization first enters at Level 6"
                    violation=1
                fi
            done
        fi

        # NOTHING IN THIS TOWER IS EVER APPLIED DIRECTLY. Level 6
        # constructs production operators and interrogates them;
        # Level 7 probes numerically, but only through the minimizer's
        # own vocabulary - matvec, sweep_order, diagonal - so the
        # measurement is about the minimizer rather than about the
        # action. A level that called apply itself would have bypassed
        # the very interface it is measuring.
        #
        # Comments are stripped first, as the numberless law strips
        # them: this scan measures CODE, and the tower's prose says
        # "% apply()" whenever it explains why it does not call one.
        if sed 's/!.*//' "$src" | grep -qE "% *apply *\("; then
            echo "IMPORT GATE: $(basename "$src") in $name applies an operation - this tower interrogates structure, it never evaluates"
            violation=1
        fi
    done
done

if [ "$violation" -ne 0 ]; then
    echo "IMPORT GATE: the tower layering is violated"
    exit 1
fi
echo "import gate: every source imports only its level and below"
