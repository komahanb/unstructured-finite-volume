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
