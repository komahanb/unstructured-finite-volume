#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [nosamples]="at least one sample is required"
  [idx0]="unknown sample index is in range"
  [idxhigh]="unknown sample index is in range"
  [missingq]="unknown sample provides q"
  [novalues]="trial q has values"
  [entries]="trial q shape matches unknown q"
  [comp]="trial q shape matches unknown q"
  [storage]="trial q shape matches unknown q"
  [badfield]="q field dynamic type is supported"
)
for case in nosamples idx0 idxhigh missingq novalues entries comp storage badfield; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
