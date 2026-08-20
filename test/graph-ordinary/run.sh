#!/bin/bash
# build the library and the ordinary profile suite, run the laws,
# then run every refusal and assert each dies for its stated reason.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

# The schema refusals retired with graph_profile: the stored graph
# makes a two-tailed or two-headed edge unrepresentable rather than
# refused. This guard keeps the profile deleted: the file must not
# return and nothing may import it.
if [ -e "$here/../../src/graph_profile.f90" ]; then
    echo " FAIL : the deleted graph_profile has returned"
    exit 1
fi
grep -q "use graph_profile" "$here"/../../src/*.f90 "$here"/../../test/*/*.f90 "$here"/../../test/*/*/*.f90 2>/dev/null \
    && { echo " FAIL : a reference to the deleted graph_profile exists"; exit 1; } \
    || echo " PASS : graph_profile stays deleted"

# The naming law (doc/coding-standards.md): every src module reads
# namespace order - its first word is a prime or util - and no type
# reads namespace order, so a type can never collide with a module.
bad_modules=$(ls "$here"/../../src/*.f90 | xargs -n1 basename | sed 's/\.f90//' \
    | grep -vE '^(graph|token|relation|field|operation|transform|map|view|util)_' || true)
if [ -n "$bad_modules" ]; then
    echo " FAIL : modules outside the prime namespace: $bad_modules"
    exit 1
fi
bad_types=$(grep -hE '^ *type(, *(abstract|extends\([a-z_]+\)|public|private))* *:: *[a-z_]+' \
        "$here"/../../src/*.f90 \
    | grep -oE ':: *[a-z_]+' | sed 's/:: *//' | sort -u \
    | grep -E '^(graph|token|relation|field|operation|transform|map|view|util)_' || true)
if [ -n "$bad_types" ]; then
    echo " FAIL : types wearing a module namespace: $bad_types"
    exit 1
fi
echo " PASS : modules read namespace order, types read english order"

# OBJECTS is complete and dependency-ordered: every entry has a
# source, every source is listed, and a file may use only modules
# already built before it (the build has no dependency tracking, so
# the list order IS the dependency proof).
objects="$here/../../src/OBJECTS"
srcdir="$here/../../src"
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

# Joined view of a source: comments stripped, continuations merged,
# so a name split across & lines cannot hide from the checks below.
joined() {
    sed 's/!.*//' "$1" | sed -e ':a' -e '/&[ ]*$/{N;s/&[ ]*\n[ ]*&\{0,1\}/ /;ta}'
}

# One type, one name: no import alias in src except dp => real64.
for f in "$srcdir"/*.f90; do
    if joined "$f" | grep -iE 'only *:.*=>' | grep -ivE 'dp *=> *real64' | grep -q .; then
        echo " FAIL : $(basename "$f") carries an import alias other than dp => real64"
        exit 1
    fi
done
echo " PASS : one type, one name - no import aliases in src"

# Readers are bare nouns: no get_ declaration in src.
if grep -nE '(procedure[^!]*:: *get_|function +get_|subroutine +get_|public[^!]*\bget_)' "$srcdir"/*.f90 | grep -q .; then
    grep -nE '(procedure[^!]*:: *get_|function +get_|subroutine +get_|public[^!]*\bget_)' "$srcdir"/*.f90 | head -3
    echo " FAIL : a get_ declaration exists in src"
    exit 1
fi
echo " PASS : readers are bare nouns - no get_ declarations in src"

# No module re-exports a name it does not define: a public name that
# also appears in the module's own use-only imports is a re-export.
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
