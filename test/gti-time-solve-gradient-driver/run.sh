#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [fwdopt]="max_iterations is positive"
  [revopt]="reverse_driver: singular tolerance is positive"
  [badgraph]="gti_time_graph: relation unknown sample is in range"
  [zeroterm]="gti_time_functional: at least one term is required"
  [noeta]="design direction has values"
)
for case in fwdopt revopt badgraph zeroterm noeta; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
# The composition's own two scalar guards protect against future
# lower-driver changes; every lower driver produces those actions
# as scalars by construction, so no public path can trip the
# guards today. Their presence in source is proven instead.
grep -q "gti_time_solve_gradient_driver: explicit design action is scalar" \
    "$here/../../src/gti_time_solve_gradient_driver.f90" \
    && echo " PASS : the explicit-action scalar guard exists in source" \
    || { echo " FAIL : explicit-action guard missing"; exit 1; }
grep -q "gti_time_solve_gradient_driver: residual design action is scalar" \
    "$here/../../src/gti_time_solve_gradient_driver.f90" \
    && echo " PASS : the residual-action scalar guard exists in source" \
    || { echo " FAIL : residual-action guard missing"; exit 1; }
