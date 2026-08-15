#!/bin/bash
# build the library and the structure suite, run the laws, then run
# every refusal and assert that each dies for its stated reason.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

declare -A reason=(
  [unheld-domain]="a relation must relate the graph's own member sets"
  [empty-relation-family]="a related graph declares at least one relation"
  [map-empty-family]="a related graph declares at least one relation"
  [map-unheld-domain]="a relation must relate the graph's own member sets"
  [map-view]="a view cannot be owned"
  [twice]="a graph never signs twice"
  [undeclared]="a graph holds declared domains only"
  [dupset]="a graph holds each domain once"
  [duprel]="a graph holds each relation once"
  [view]="a view cannot be owned"
)

for case in unheld-domain empty-relation-family map-empty-family map-unheld-domain map-view twice undeclared dupset duprel view; do
    if ./refusal "$case" >refusal.out 2>&1; then
        echo " FAIL : '$case' was accepted"
        exit 1
    fi
    if grep -q "${reason[$case]}" refusal.out; then
        echo " PASS : '$case' is refused, loudly"
    else
        echo " FAIL : '$case' died for the wrong reason"
        cat refusal.out
        exit 1
    fi
done
rm -f refusal.out
