#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [negiter]="max_iterations is positive"
  [negrtol]="residual tolerance is nonnegative"
  [negstep]="step tolerance is nonnegative"
  [negsing]="singular tolerance is positive"
  [novalues]="newton_driver: trial q has values"
  [dqshape]="direction q shape matches trial q"
  [wide]="residual is a vector"
  [short]="residual size matches unknown size"
  [singular]="dense Jacobian pivot is nonsingular"
)
for case in negiter negrtol negstep negsing novalues dqshape wide short singular; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
