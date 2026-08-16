#!/bin/bash
# Build the library and the identity-map suites, run the laws, then run
# each under valgrind and assert that a map outliving its binders reads
# nothing it does not own. Refusals last.
set -e

here="$(cd "$(dirname "$0")" && pwd)"

( cd "$here/../.." && ./build.sh >/dev/null )

make -C "$here" clean >/dev/null 2>&1 || true
make -C "$here" >/dev/null

cd "$here"
./lifetime
echo ''
./labels

echo ''
if command -v valgrind >/dev/null 2>&1; then
    # The law is about STORAGE, so the allocator must agree with the
    # native run. Before the gate the lifetime probe reported 170
    # invalid reads from 48 contexts and still printed ALL PROPOSITIONS
    # HOLD.
    for probe in lifetime labels; do
        if valgrind --error-exitcode=9 --leak-check=full \
                    --errors-for-leak-kinds=definite \
                    ./$probe >valgrind.out 2>&1; then
            echo " PASS : $probe reads no storage it does not own"
        else
            echo " FAIL : $probe still borrows its keys"
            grep -E "ERROR SUMMARY|Invalid read|definitely lost" valgrind.out || true
            exit 1
        fi
    done
    rm -f valgrind.out
else
    echo " SKIP : valgrind is not installed; the native runs stand alone"
fi

declare -A reason=(
  [unsigned]="keyed on assigned identity"
  [twice]="a set is named once"
)

echo ''
for case in unsigned twice; do
    if ./refusal "$case" >refusal.out 2>&1; then
        echo " FAIL : '$case' was admitted"
        exit 1
    fi
    if grep -q "${reason[$case]}" refusal.out; then
        echo " PASS : '$case' is refused"
    else
        echo " FAIL : '$case' is refused for a different reason"
        cat refusal.out
        exit 1
    fi
done
rm -f refusal.out
