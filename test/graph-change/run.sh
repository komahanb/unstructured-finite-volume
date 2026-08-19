#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [attachtwice]="graph_value_map: a value seat is attached once"
  [updatefree]="graph_value_map: an update touches an attached seat"
  [readunknown]="graph_value_map: a known value is read"
  [detachfree]="graph_value_map: a detach removes an attached seat"
  [emptyknown]="graph_value_map: a known value has values"
  [undeclared]="graph_value_map: a value map is keyed on assigned identity"
  [liarapply]="change_controller: applied change reports applied"
  [liarrevert]="change_controller: reverted change reports reverted"
  [impossible]="change_result: terminal state is consistent"
  [unbound]="value_change: value map is bound"
)
for case in attachtwice updatefree readunknown detachfree emptyknown \
            undeclared liarapply liarrevert impossible unbound; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

# the book stays closed, statically: this suite speaks only the
# living tower - no GTI vocabulary anywhere in its sources.
grep -q "gti_" "$here"/*.f90 \
    && { echo " FAIL : a GTI name appears in the suite"; exit 1; } \
    || echo " PASS : the suite speaks only the living tower"
