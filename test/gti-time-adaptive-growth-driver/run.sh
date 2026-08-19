#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [novertex]="graph has an initial vertex"
  [wrongvertex]="appended vertex is last"
  [wrongrelation]="appended relation is last"
  [nosentinel]="candidate relation references appended vertex"
  [wrongunknown]="candidate unknown is appended vertex"
  [badmotif]="gti_time_graph: motif rule weight count matches relation arity"
  [unsolvedhistory]="gti_time_forward_driver: history vertex has solution"
)
for case in novertex wrongvertex wrongrelation nosentinel wrongunknown badmotif unsolvedhistory; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
