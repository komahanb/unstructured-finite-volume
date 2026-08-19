#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [deg0]="degree_r_directional_driver: degree is supported"
  [deg5]="degree_r_directional_driver: degree is supported"
  [negsing]="degree_r_directional_driver: singular tolerance is positive"
  [norelations]="degree_r_directional_driver: at least one relation is required"
  [idx0]="degree_r_directional_driver: relation index is in range"
  [idxhigh]="degree_r_directional_driver: relation index is in range"
  [degdim]="degree_r_directional_driver: derivative array degree matches options"
  [vertdim]="degree_r_directional_driver: derivative array vertex count matches graph"
  [noeta]="degree_r_directional_driver: design direction has values"
  [unsolved]="degree_r_directional_driver: unknown vertex has solution"
  [unsolvedhistory]="degree_r_directional_driver: history vertex has solution"
  [noq]="degree_r_directional_driver: unknown q has values"
  [badseed]="degree_r_directional_driver: derivative shape matches vertex q"
  [shortres]="degree_r_directional_driver: degree residual size matches unknown size"
  [singular]="degree_r_directional_driver: dense Jacobian pivot is nonsingular"
)
for case in deg0 deg5 negsing norelations idx0 idxhigh degdim vertdim noeta \
            unsolved unsolvedhistory noq badseed shortres singular; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

src="$here/../../src/gti_time_degree_r_directional_driver.f90"

# check 23: the higher-order chain-rule seat is genuinely the
# assembler of every B^(s) - the import is in the production source.
grep -q "use gti_chain_rule_assemblies" "$src" \
    && echo " PASS : the production source imports gti_chain_rule_assemblies" \
    || { echo " FAIL : chain-rule assembler import missing"; exit 1; }

# check 24: the general driver supersedes the degree-2 seat without
# leaning on it - the old driver's name is absent from production.
grep -q "gti_time_degree2_directional_drivers" "$src" \
    && { echo " FAIL : the production source names the degree-2 driver"; exit 1; } \
    || echo " PASS : the production source does not import the degree-2 driver"

# check 25: no primal residual evaluation exists in the source, so
# nothing can be differenced - every number is a partial action.
grep -qF "% value(" "$src" \
    && { echo " FAIL : the production source evaluates a primal value"; exit 1; } \
    || echo " PASS : no primal evaluation, hence no finite differencing, exists"

# check 26: J_u is built once, before the degree loop - textually.
jline=$(grep -n "call build_unknown_jacobian" "$src" | cut -d: -f1)
sline=$(grep -n "do s = 1, options % max_degree" "$src" | cut -d: -f1)
[ "$(echo "$jline" | wc -l)" -eq 1 ] && [ -n "$jline" ] && [ -n "$sline" ] && [ "$jline" -lt "$sline" ] \
    && echo " PASS : J_u is built once per relation, before the degree loop" \
    || { echo " FAIL : the Jacobian build does not precede the degree loop"; exit 1; }
