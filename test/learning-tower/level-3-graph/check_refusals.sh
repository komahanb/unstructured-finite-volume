#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"
declare -A reason=(
  [foreign]="a relation must relate the graph's own member sets"
  [dupset]="a graph holds each domain once"
  [duprel]="a graph holds each relation once"
  [view]="a view cannot be owned"
)
for case in foreign dupset duprel view; do
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
