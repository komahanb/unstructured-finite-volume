#!/bin/bash
# Gate C's refusals: the composite belongs to ONE decomposition.
# A same-sized graph of a different shape must die, and die for its
# stated reason - the message is checked, not merely the exit code.
set -e

here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"

declare -A reason=(
  [unequal-host-apply]="this decomposition belongs to another graph"
  [unequal-host-attach]="this decomposition belongs to another graph"
)

for case in unequal-host-apply unequal-host-attach; do
    if ./refusal "$case" >refusal.out 2>&1; then
        echo " FAIL : '$case' was accepted"
        exit 1
    fi
    if grep -q "${reason[$case]}" refusal.out; then
        echo " PASS : '$case' is refused, loudly"
    else
        echo " FAIL : '$case' died for the wrong reason"
        cat refusal.out
        exit 1
    fi
done
rm -f refusal.out
