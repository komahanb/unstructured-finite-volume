#!/bin/bash
# Level 8 refusal: a doubly-located residual row must die by name.
set -e
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"
if ./refusal >refusal.out 2>&1; then
    echo " FAIL : a doubly-located row was accepted"
    exit 1
fi
grep -q "one location per residual row" refusal.out \
    && echo " PASS : a doubly-located row is refused, loudly" \
    || { cat refusal.out; exit 1; }
rm -f refusal.out
