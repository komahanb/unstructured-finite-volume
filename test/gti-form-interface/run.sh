#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [rdegree]="beyond the declared max_degree"
  [fdegree]="beyond the declared max_degree"
  [dircount]="one direction per order"
  [mismatch]="must tell one story"
  [badkind]="unknown argument kind"
  [badcomp]="unknown state component"
  [bshape]="must fill entries times components exactly"
)
for case in rdegree fdegree dircount mismatch badkind badcomp bshape; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
