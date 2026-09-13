#!/bin/bash
# Verify immutable topology, incidence laws and allocation-free readers.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
if [ "${UFVM_SKIP_LIBRARY_BUILD:-0}" != 1 ]; then
    (cd "$here/../.." && ./build.sh >/dev/null)
fi
make -C "$here" clean >/dev/null
make -C "$here" >/dev/null
(cd "$here" && ./run && ./traversal 2000 10)
make -C "$here" refusals
declare -A reason=(
    [negative_count]='the vertex count must be nonnegative'
    [endpoint_extents]='head and tail arrays must have equal extent'
    [zero_tail]='every tail must belong to the vertex set'
    [outside_tail]='every tail must belong to the vertex set'
    [vertex_tags]='one tag is required per vertex'
    [edge_tags]='one tag is required per edge'
    [global_vertex_extent]='one global index is required per vertex'
    [global_vertex_range]='every global vertex index must belong to the whole carrier'
    [global_vertex_duplicates]='global vertex indices must be injective'
    [global_edge_extent]='one global index is required per edge'
    [global_edge_range]='every global edge index must belong to the whole carrier'
    [global_edge_duplicates]='global edge indices must be injective'
    [vertex_owner_extent]='one owner is required per vertex'
    [vertex_owner_range]='every vertex owner must belong to the set of parts'
    [edge_owner_extent]='one owner is required per edge'
    [edge_owner_range]='every edge owner must belong to the set of parts'
    [missing_partition_vertices]='was passed without vglobal'
    [partition_number]='the part number must lie in 1..np'
    [undeclared_whole_vertices]='whole_vertices was passed but its identity is undeclared'
    [undeclared_whole_edges]='whole_edges was passed but its identity is undeclared'
    [inconsistent_vertex_extent]="the vertex carrier's extent must equal nv"
    [inconsistent_edge_extent]="the edge carrier's extent must equal size(tails)"
    [boundary_transpose]='this graph cannot be transposed'
    [boundary_reverse]='this graph cannot be transposed'
    [zero_fibre_index]="the index argument lies outside the fibre's extent"
    [outside_fibre_index]="the index argument lies outside the fibre's extent"
    [empty_fibre_index]='fibre_member was called on a fibre with no entries associated'
    [hierarchy_shared_extension]='requires a hierarchy with a sole owner'
    [hierarchy_released_twin]="this storage's hierarchy has been released"
    [hierarchy_empty_index]='requires an allocated hierarchy, but this storage has none'
)
refusal_output="$(mktemp)"
trap 'rm -f "$refusal_output"' EXIT
for case_name in "${!reason[@]}"; do
    if "$here/refusal" "$case_name" >"$refusal_output" 2>&1; then
        echo "FAIL : invalid topology accepted: $case_name"
        exit 1
    fi
    if ! grep -Fq "${reason[$case_name]}" "$refusal_output"; then
        cat "$refusal_output"
        echo "FAIL : unrelated refusal: $case_name"
        exit 1
    fi
    echo "PASS : invalid topology refused: $case_name"
done
