#!/bin/bash
set -e
suite_dir="$(cd "$(dirname "$0")" && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null 2>&1 )
fi
make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null
cd "$suite_dir" && ./run
