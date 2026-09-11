#!/bin/bash
# Build the library and the partition-law suite, then run the laws.
set -e

suite_dir="$(cd "$(dirname "$0")" && pwd)"

if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    ( cd "$suite_dir/../.." && ./build.sh >/dev/null )
fi

make -C "$suite_dir" clean >/dev/null 2>&1 || true
make -C "$suite_dir" >/dev/null

cd "$suite_dir" && ./run

declare -A expected=(
  [wrong-whole]="assemble: this relation was not written for this whole"
  [smaller-whole]="assemble: this relation was not written for this whole"
  [whole-vertex-count]="assemble: this relation was not written for this whole"
  [whole-edge-count]="assemble: this relation was not written for this whole"
  [part-count]="assemble: this relation was not written for this part"
  [out-of-range-full]="assemble: a relation must map into the whole carrier"
  [out-of-range-subset]="assemble: a relation must map into the whole carrier"
  [field-count]="assemble: a full field must fill the part carrier"
)
for mode in wrong-whole smaller-whole whole-vertex-count whole-edge-count \
            part-count out-of-range-full out-of-range-subset field-count; do
  if output=$(./refusal "$mode" 2>&1); then
    echo "FAIL: $mode was accepted"
    exit 1
  fi
  case "$output" in
    *"${expected[$mode]}"*) echo "PASS: $mode refused" ;;
    *) echo "FAIL: $mode failed for another reason"; echo "$output"; exit 1 ;;
  esac
done
