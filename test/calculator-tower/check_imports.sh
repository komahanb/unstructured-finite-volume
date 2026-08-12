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
        common)                  echo "" ;;
        level-0-carrier)         echo "calculator_assert graph_carrier" ;;
        level-1-relation)        echo "calculator_assert graph_carrier graph_relation" ;;
        level-2-relation-algebra) echo "calculator_assert graph_carrier graph_relation graph_relation_algebra" ;;
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
