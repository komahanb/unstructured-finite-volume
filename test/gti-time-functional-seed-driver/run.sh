#!/bin/bash
set -e
here="$(cd "$(dirname "$0")" && pwd)"
( cd "$here/../.." && ./build.sh >/dev/null 2>&1 )
make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null
cd "$here" && ./run
declare -A reason=(
  [zeroterm]="at least one term is required"
  [vertex0]="term vertex is in range"
  [vertexhigh]="term vertex is in range"
  [unsolved]="term vertex has solution"
  [novertexq]="term vertex provides q"
  [seedsize]="seed_driver: vertex seed size matches graph vertices"
  [noeta]="seed_driver: design direction has values"
  [qnovalues]="seed_driver: vertex q has values"
  [valshort]="functional value is scalar"
  [valwide]="functional value is scalar"
  [gradshort]="seed_driver: functional action is scalar"
  [gradwide]="seed_driver: functional action is scalar"
  [dashort]="seed_driver: functional action is scalar"
  [dawide]="seed_driver: functional action is scalar"
  [propshape]="seed_driver: propagated seed shape matches target"
)
for case in zeroterm vertex0 vertexhigh unsolved novertexq seedsize noeta qnovalues valshort valwide gradshort gradwide dashort dawide propshape; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
