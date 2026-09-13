#!/bin/bash
# Level 1 refusals: each case must die, and die for its stated
# reason - the message is checked, not merely the nonzero exit.
set -e

here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"

declare -A reason=(
  [arity]="requires one row per domain"
  [member]="requires every tuple entry to belong to its domain"
  [undeclared]="requires a signature of declared domains only"
)

for case in arity member undeclared; do
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
