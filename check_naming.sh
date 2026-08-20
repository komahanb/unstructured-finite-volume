#!/bin/bash
# The naming law (doc/coding-standards.md), statically enforced.
#
#   usage: check_naming.sh [tree]
#
# tree is a directory holding src/ (default: this script's own
# directory). Two callers share this one implementation, so the law
# has exactly one set of teeth:
#
#   githooks/pre-commit           runs it on the STAGED tree
#   test/graph-ordinary/run.sh    runs it on the working tree
#
# Each check states what it refuses. A check that cannot be made
# mechanical (english word order, num_ for a quantity that IS a
# cardinality) lives in review, not here.
set -e

root="${1:-$(cd "$(dirname "$0")" && pwd)}"
srcdir="$root/src"
objects="$srcdir/OBJECTS"

[ -d "$srcdir" ] || { echo " FAIL : no src/ under $root"; exit 1; }

# Comments stripped and continuations merged, so a name split
# across & lines cannot hide from the token checks below.
joined() {
    sed 's/!.*//' "$1" | sed -e ':a' -e '/&[ ]*$/{N;s/&[ ]*\n[ ]*&\{0,1\}/ /;ta}'
}

#---------------------------------------------------------------------
# 1. Modules read namespace order: every src file begins with a prime
#    (or util_). Types read english order: no type wears a namespace
#    prefix, so a type can never collide with a module name.
#---------------------------------------------------------------------

primes='graph|token|relation|field|operation|transform|map|view|util'

bad_modules=$(ls "$srcdir"/*.f90 | xargs -n1 basename | sed 's/\.f90//' \
    | grep -vE "^($primes)_" || true)
if [ -n "$bad_modules" ]; then
    echo " FAIL : modules outside the prime namespace: $bad_modules"
    exit 1
fi

bad_types=$(grep -hE '^ *type(, *(abstract|extends\([a-z_]+\)|public|private))* *:: *[a-z_]+' \
        "$srcdir"/*.f90 \
    | grep -oE ':: *[a-z_]+' | sed 's/:: *//' | sort -u \
    | grep -E "^($primes)_" || true)
if [ -n "$bad_types" ]; then
    echo " FAIL : types wearing a module namespace: $bad_types"
    exit 1
fi
echo " PASS : modules read namespace order, types read english order"

#---------------------------------------------------------------------
# 2. OBJECTS is complete and dependency-ordered: every entry has a
#    source, every source is listed, and a file may use only modules
#    already built before it. The build has no dependency tracking -
#    the list order IS the dependency proof - and the incremental
#    build's stale .mod files can mask a clean-build order bug, which
#    is why this check exists.
#---------------------------------------------------------------------

listed=$(grep -oE '^[a-z_0-9]+\.o' "$objects" | sed 's/\.o$//')
for m in $listed; do
    [ -f "$srcdir/$m.f90" ] || { echo " FAIL : OBJECTS lists $m.o with no $m.f90"; exit 1; }
done
for f in "$srcdir"/*.f90; do
    b=$(basename "$f" .f90)
    echo "$listed" | grep -qx "$b" || { echo " FAIL : $b.f90 is not in OBJECTS"; exit 1; }
done
seen=" "
for m in $listed; do
    for u in $(sed 's/!.*//' "$srcdir/$m.f90" | grep -oE '^ *use +[a-z_0-9]+' | awk '{print $2}' | sort -u); do
        [ -f "$srcdir/$u.f90" ] || continue
        case "$seen" in
        *" $u "*) : ;;
        *) echo " FAIL : OBJECTS builds $m before its import $u"; exit 1 ;;
        esac
    done
    seen="$seen$m "
done
echo " PASS : OBJECTS is complete and dependency-ordered"

#---------------------------------------------------------------------
# 3. One type, one name: no import alias in src except the kind
#    constant dp => real64, the law's one named exception.
#---------------------------------------------------------------------

for f in "$srcdir"/*.f90; do
    if joined "$f" | grep -iE 'only *:.*=>' | grep -ivE 'dp *=> *real64' | grep -q .; then
        echo " FAIL : $(basename "$f") carries an import alias other than dp => real64"
        exit 1
    fi
done
echo " PASS : one type, one name - no import aliases in src"

#---------------------------------------------------------------------
# 4. Readers are bare nouns: no get_ declaration in src.
#---------------------------------------------------------------------

if grep -nE '(procedure[^!]*:: *get_|function +get_|subroutine +get_|public[^!]*\bget_)' \
        "$srcdir"/*.f90 | grep -q .; then
    grep -nE '(procedure[^!]*:: *get_|function +get_|subroutine +get_|public[^!]*\bget_)' \
        "$srcdir"/*.f90 | head -3
    echo " FAIL : a get_ declaration exists in src"
    exit 1
fi
echo " PASS : readers are bare nouns - no get_ declarations in src"

#---------------------------------------------------------------------
# 5. No module re-exports a name it does not define: a public name
#    that also appears in the module's own use-only imports is a
#    re-export.
#---------------------------------------------------------------------

for f in "$srcdir"/*.f90; do
    imports=$(joined "$f" | grep -iE '^ *use .*only *:' | sed 's/.*only *://I' \
        | tr ',' '\n' | sed 's/=>.*//; s/ //g' | grep -v '^$' | sort -u)
    publics=$(joined "$f" | grep -E '^ *public *::' | sed 's/.*:://' \
        | tr ',' '\n' | sed 's/ //g' | grep -v '^$' | sort -u)
    both=$(comm -12 <(echo "$imports") <(echo "$publics") | grep -v '^$' || true)
    if [ -n "$both" ]; then
        echo " FAIL : $(basename "$f") re-exports a name it does not define: $both"
        exit 1
    fi
done
echo " PASS : no module re-exports a name it does not define"

#---------------------------------------------------------------------
# 6. British spelling, oxford -ize: the closed list of refused
#    spellings (identifiers and comments alike).
#---------------------------------------------------------------------

if grep -inwE 'fiber|fibers|center|centers|centered|coloring|colored|neighbor|neighbors' \
        "$srcdir"/*.f90 | grep -q .; then
    grep -inwE 'fiber|fibers|center|centers|centered|coloring|colored|neighbor|neighbors' \
        "$srcdir"/*.f90 | head -3
    echo " FAIL : an american spelling exists in src (fibre, centre, colouring, neighbour)"
    exit 1
fi
echo " PASS : british spelling holds (fibre, centre, colouring, neighbour)"

echo " PASS : the naming law holds"
