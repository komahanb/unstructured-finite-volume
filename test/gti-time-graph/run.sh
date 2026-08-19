#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [novertices]="at least one vertex is required"
  [buildidx0]="relation index is in range"
  [buildidxhigh]="relation index is in range"
  [novertextuple]="relation has sample vertices"
  [vertexidx0]="relation sample vertex is in range"
  [vertexidxhigh]="relation sample vertex is in range"
  [unknown0]="relation unknown sample is in range"
  [unknownhigh]="relation unknown sample is in range"
  [norules]="relation motif has rules"
  [noweights]="motif rule has weights"
  [weightcount]="motif rule weight count matches relation arity"
  [missingq]="referenced vertex provides q"
)
for case in novertices buildidx0 buildidxhigh novertextuple vertexidx0 vertexidxhigh unknown0 unknownhigh norules noweights weightcount missingq; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
