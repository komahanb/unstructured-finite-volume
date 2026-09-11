#!/bin/bash
# Build the library and the map prototypes, then run them. Analysis
# evidence for the member-set and relation fractal maps - not
# production laws.
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null )
fi

make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null

cd "$suite_dir"
./set
echo ''
./relation
echo ''
./scale
