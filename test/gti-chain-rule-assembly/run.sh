#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [negdeg]="degree is supported"
  [beyondform]="a partial request beyond the declared max_degree is refused"
  [hugedeg]="pattern coefficient is representable"
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
for case in negdeg beyondform hugedeg inconsistent dupq dupxi negent zerocomp badseedmeta inconsistent4 shape4 shifty; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

src="$here/../../src/gti_chain_rule_assembly.f90"

# the retirement, statically: no hand table, no artificial cap.
grep -qE "degree > 4|degree <= 4|max_degree > 4" "$src" \
    && { echo " FAIL : a degree-4 cap survives in the source"; exit 1; } \
    || echo " PASS : no degree-4 cap exists in the source"
grep -q "select case (degree)" "$src" \
    && { echo " FAIL : a tabulated select case over the degree survives"; exit 1; } \
    || echo " PASS : the pattern table is generated, not tabulated"
