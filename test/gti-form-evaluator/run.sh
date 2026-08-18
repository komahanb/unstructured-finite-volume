#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [siglen]="output_signature has shape"
  [negent]="output nentries is nonnegative"
  [zerocomp]="output ncomp is at least one"
  [valueshape]="form output does not match output_signature"
  [partialshape]="form output does not match output_signature"
  [storage]="output storage does not match output shape"
)
for case in siglen negent zerocomp valueshape partialshape storage; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
