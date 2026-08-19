#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [negsing]="singular tolerance is positive"
  [noqstar]="q star has values"
  [noeta]="design direction has values"
  [wide]="tangent rhs is a vector"
  [rhssize]="rhs size matches unknown size"
  [singular]="tangent_driver: dense Jacobian pivot is nonsingular"
)
for case in negsing noqstar noeta wide rhssize singular; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
