#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [unbound]="gti_attached_value_change: value map is bound"
  [unboundrun]="gti_attached_value_change: value map is bound"
  [emptyvalue]="gti_attached_value_map: a known value has values"
  [undeclaredgraph]="gti_attached_value_map: a value map is keyed on assigned identity"
)
for case in unbound unboundrun emptyvalue undeclaredgraph; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

src="$here/../../src/gti_attached_value_change.f90"

# the client law, statically: the change rides the protocol and the
# map, defines no controller of its own, and names no time module,
# no solver, no mesh, no form, no evaluator, no driver.
grep -q "use gti_change_protocols" "$src" \
    && echo " PASS : the change imports gti_change_protocols" \
    || { echo " FAIL : the protocol import is missing"; exit 1; }
grep -q "use gti_attached_value_maps" "$src" \
    && echo " PASS : the change imports gti_attached_value_maps" \
    || { echo " FAIL : the map import is missing"; exit 1; }
grep -v '^ *!' "$src" | grep -q "gti_change_controller" \
    && { echo " FAIL : the change defines or imports a controller"; exit 1; } \
    || echo " PASS : no controller is defined or imported here"
grep -qE "gti_time_|use .*(solver|mesh|flux|form|evaluator|driver)" "$src" \
    && { echo " FAIL : a forbidden seat is imported"; exit 1; } \
    || echo " PASS : no time, solver, mesh, form, evaluator, or driver import"
