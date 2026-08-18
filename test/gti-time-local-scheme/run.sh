#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [nosamples]="at least one sample is required"
  [unknowncomp]="unknown state component is refused"
  [dupq]="duplicate state component rule is refused"
  [noweights]="a rule has weights"
  [wrongcount]="rule weight count matches sample count"
  [missingq]="every sample provides q"
  [qshape]="sample q shapes match"
)
for case in nosamples unknowncomp dupq noweights wrongcount missingq qshape; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
