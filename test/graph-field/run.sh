#!/bin/bash
set -e
suite_dir="$(cd "$(dirname "$0")" && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null 2>&1 )
fi
make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null
cd "$suite_dir" && ./run
declare -A reason=(
  [ishape]="a value vector must fill its domain exactly"
  [rshape]="a value vector must fill its domain exactly"
  [cshape]="a value vector must fill its domain exactly"
  [lshape]="a value vector must fill its domain exactly"
  [sshape]="a value vector must fill its domain exactly"
  [unsigned]="a field requires a declared domain"
)
for case in ishape rshape cshape lshape sshape unsigned; do
    if ./refusal "$case" >refusal.out 2>&1; then echo " FAIL : '$case' accepted"; exit 1; fi
    grep -q "${reason[$case]}" refusal.out && echo " PASS : '$case' is refused, loudly" || { cat refusal.out; exit 1; }
done
rm -f refusal.out
