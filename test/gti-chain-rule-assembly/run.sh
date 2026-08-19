#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [negdeg]="degree is supported"
  [highdeg]="degree is supported"
  [inconsistent]="inconsistent argument channel"
  [dupq]="duplicate argument channel"
  [dupxi]="duplicate argument channel"
  [negent]="output nentries is nonnegative"
  [zerocomp]="output ncomp is at least one"
  [badseedmeta]="derivative seat names its argument"
  [inconsistent4]="inconsistent argument channel"
  [shape4]="form output does not match output_signature"
  [shifty]="output accumulation shape mismatch"
)
for case in negdeg highdeg inconsistent dupq dupxi negent zerocomp badseedmeta inconsistent4 shape4 shifty; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
