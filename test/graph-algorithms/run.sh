#!/bin/bash
# build the library and the algorithms suite, run the laws, then run
# every refusal and assert each terminates for its stated reason.
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null )
fi

make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null

cd "$suite_dir" && ./run

declare -A reason=(
  [notbinary]="relation_algorithms: the adjacency is a binary relation"
  [notsquare]="relation_algorithms: the adjacency runs over one domain"
  [cycle]="a topological order needs an acyclic graph"
)

for case in notbinary notsquare cycle; do
    if ./refusal "$case" >refusal.out 2>&1; then
        echo " FAIL : '$case' was accepted"
        exit 1
    fi
    if grep -q "${reason[$case]}" refusal.out; then
        echo " PASS : '$case' is refused, loudly"
    else
        echo " FAIL : '$case' terminated for the wrong reason"
        cat refusal.out
        exit 1
    fi
done
rm -f refusal.out
