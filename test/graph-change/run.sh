#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [attachtwice]="map_value: a value row is attached once"
  [updatefree]="map_value: an update touches an attached row"
  [readunknown]="map_value: a known value is read"
  [detachfree]="map_value: a detach removes an attached row"
  [emptyknown]="map_value: a known value has values"
  [undeclared]="map_value: a value map is keyed on assigned identity"
  [silentapply]="run_change: applied change reports applied"
  [silentrevert]="run_change: reverted change reports reverted"
  [impossible]="change_record: terminal state is consistent"
  [unbound]="value_change: value map is bound"
)
for case in attachtwice updatefree readunknown detachfree emptyknown \
            undeclared silentapply silentrevert impossible unbound; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

# Regression guard: the test sources must not reference the
# deleted gti_* modules.
grep -q "gti_" "$here"/*.f90 \
    && { echo " FAIL : a gti_ reference appears in the suite sources"; exit 1; } \
    || echo " PASS : no gti_ reference in the suite sources"
