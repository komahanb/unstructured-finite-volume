#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [norelations]="at least one relation is required"
  [idx0]="forward_driver: relation index is in range"
  [idxhigh]="forward_driver: relation index is in range"
  [unsolvedhistory]="history vertex has solution"
  [noq]="unknown initial q has values"
  [ncomp2]="unknown initial q is a vector"
)
for case in norelations idx0 idxhigh unsolvedhistory noq ncomp2; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
