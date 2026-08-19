#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [attachundeclared]="gti_attached_value_map: a value map is keyed on assigned identity"
  [attachtwice]="gti_attached_value_map: a value seat is attached once"
  [markknownundeclared]="gti_attached_value_map: a value map is keyed on assigned identity"
  [markknownunattached]="gti_attached_value_map: an update touches an attached seat"
  [markunknownunattached]="gti_attached_value_map: an update touches an attached seat"
  [detachundeclared]="gti_attached_value_map: a value map is keyed on assigned identity"
  [detachunattached]="gti_attached_value_map: a detach removes an attached seat"
  [readunattached]="gti_attached_value_map: a known value is read"
  [readunknown]="gti_attached_value_map: a known value is read"
  [knownempty]="gti_attached_value_map: a known value has values"
)
for case in attachundeclared attachtwice markknownundeclared markknownunattached \
            markunknownunattached detachundeclared detachunattached \
            readunattached readunknown knownempty; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

src="$here/../../src/gti_attached_value_map.f90"

# the identity law, statically: rows key on graph_identity tokens,
# and no integer vertex position exists anywhere in this seat.
grep -q "use graph_identity" "$src" \
    && echo " PASS : the map keys on graph_identity tokens" \
    || { echo " FAIL : the identity import is missing"; exit 1; }
grep -qi "vertex" "$src" \
    && { echo " FAIL : the map names a vertex position"; exit 1; } \
    || echo " PASS : no vertex position exists in the map's vocabulary"

# not a time-graph refactor: the map never names a time module.
grep -q "gti_time_" "$src" \
    && { echo " FAIL : the map names a time module"; exit 1; } \
    || echo " PASS : the map names no GTI time module"
