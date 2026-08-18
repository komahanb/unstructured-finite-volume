#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [bdf0]="BDF order is supported"
  [bdf3]="BDF order is supported"
  [bdfh0]="step size is positive"
  [bdfhneg]="step size is positive"
  [dirkh0]="step size is positive"
  [dirkg0]="DIRK diagonal is nonzero"
  [abm0]="ABM order is supported"
  [abm3]="ABM order is supported"
  [abmh0]="step size is positive"
)
for case in bdf0 bdf3 bdfh0 bdfhneg dirkh0 dirkg0 abm0 abm3 abmh0; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
