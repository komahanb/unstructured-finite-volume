#!/bin/bash
# build the library and the fractal spike, run the laws and the four
# inhabitations, then run every refusal and assert that each dies
# for its stated reason.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here" && ./run

declare -A reason=(
  [unknownanswer]="an unknown branch answers no realization - bottom is not a value"
  [nullanswer]="a null branch answers no realization - absence is not a value"
  [realizenull]="absence is not ignorance - a null branch never realizes"
  [realizetwice]="knowledge grows once - a known branch never realizes twice"
  [foreign]="a branch realizes within its own universe"
  [unhatched]="an unhatched graph answers nothing"
  [novalue]="this graph carries no value"
  [fourth]="there is no fourth status"
  [cyclelength]="this spine never ends - a cycle has no extent"
  [unknowntail]="an unknown tail has no extent yet"
  [boundaryhead]="a boundary occurrence heads nowhere"
  [unbound]="this environment binds no value to that member"
  [lawless]="this reading binds no law to that operation"
  [cyclicexpr]="this expression never bottoms out - a cycle has no value"
  [unlearnedfar]="an unlearned far side has no answer yet"
  [numeralmember]="a compressed numeral is never a member"
  [blockzero]="a block counts from one"
  [staleuniverse]="this universe is not the one that minted you"
)

for case in unknownanswer nullanswer realizenull realizetwice foreign \
            unhatched novalue fourth cyclelength unknowntail \
            boundaryhead unbound lawless cyclicexpr unlearnedfar \
            numeralmember blockzero staleuniverse; do
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
