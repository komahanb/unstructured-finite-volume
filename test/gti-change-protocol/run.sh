#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [noapply]="gti_change_controller: applied change reports applied"
  [nocheck]="gti_change_controller: checked change reports checked"
  [nokeep]="gti_change_controller: kept change reports kept"
  [norevert]="gti_change_controller: reverted change reports reverted"
  [failapply_norevert]="gti_change_controller: reverted change reports reverted"
  [failcheck_norevert]="gti_change_controller: reverted change reports reverted"
  [vetocheck_norevert]="gti_change_controller: reverted change reports reverted"
  [badterminal]="gti_change_result: terminal state is consistent"
)
for case in noapply nocheck nokeep norevert failapply_norevert \
            failcheck_norevert vetocheck_norevert badterminal; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out

src="$here/../../src/gti_change_protocol.f90"
growth="$here/../../src/gti_time_adaptive_growth_driver.f90"

# check 11: the controller owns only the lifecycle - the run verb
# never inspects what the change declares it touches.
sed -n '/subroutine run(/,/end subroutine run/p' "$src" | grep -q 'touches_' \
    && { echo " FAIL : the controller inspects the touch flags"; exit 1; } \
    || echo " PASS : the controller never inspects touches_structure/touches_value"

# check 12: lifecycle only - the protocol imports nothing at all,
# and in particular no GTI time module and no graph/mesh/solver seat.
grep -qE '^ *use ' "$src" \
    && { echo " FAIL : the protocol imports something"; exit 1; } \
    || echo " PASS : the protocol source has no imports at all"
grep -q 'gti_time_' "$src" \
    && { echo " FAIL : the protocol names a time module"; exit 1; } \
    || echo " PASS : the protocol source names no GTI time module"

# check 13: adaptive growth is the first concrete client.
grep -q 'use gti_change_protocols' "$growth" \
    && echo " PASS : adaptive growth imports gti_change_protocols" \
    || { echo " FAIL : adaptive growth does not import the protocol"; exit 1; }
