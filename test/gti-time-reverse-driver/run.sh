#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [negsing]="reverse_driver: singular tolerance is positive"
  [norelations]="reverse_driver: at least one relation is required"
  [idx0]="reverse_driver: relation index is in range"
  [idxhigh]="reverse_driver: relation index is in range"
  [seedsize]="vertex seed size matches graph vertices"
  [noeta]="reverse_driver: design direction has values"
  [unsolved]="reverse_driver: unknown vertex has solution"
  [noseed]="unknown seed has values"
  [seedshape]="unknown seed size matches unknown size"
  [nolambda]="reverse_driver: lambda has values"
  [propshape]="propagated seed shape matches target"
  [rdsize]="reverse_driver: residual design action size matches unknown size"
  [singular]="reverse_driver: dense Jacobian pivot is nonsingular"
)
for case in negsing norelations idx0 idxhigh seedsize noeta unsolved noseed seedshape nolambda propshape rdsize singular; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
