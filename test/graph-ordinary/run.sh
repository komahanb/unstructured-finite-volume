#!/bin/bash
# build the library and the directed reading suite, run the laws,
# then the deletion guard and the naming law.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$here/../.." && ./build.sh >/dev/null )
fi

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


# The naming law, one implementation shared with the pre-commit
# hook: check_naming.sh refuses modules outside the prime namespace,
# types with a module namespace, an incomplete or misordered OBJECTS,
# import aliases, get_ declarations, re-exports, american forms
# and colloquial vocabulary.
bash "$here/../../check_naming.sh" "$here/../.."
