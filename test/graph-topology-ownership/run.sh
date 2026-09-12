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
    [negative_count]='the vertex count is nonnegative'
    [endpoint_extents]='head and tail arrays have equal extent'
    [zero_tail]='every tail belongs to the vertex set'
    [outside_tail]='every tail belongs to the vertex set'
    [vertex_tags]='one tag per vertex'
    [edge_tags]='one tag per edge'
    [global_vertex_extent]='one global index per vertex'
    [global_vertex_range]='global vertices belong to the whole carrier'
    [global_vertex_duplicates]='global vertex indices are injective'
    [global_edge_extent]='one global index per edge'
    [global_edge_range]='global edges belong to the whole carrier'
    [global_edge_duplicates]='global edge indices are injective'
    [vertex_owner_extent]='one owner per vertex'
    [vertex_owner_range]='vertex owners belong to the set of parts'
    [edge_owner_extent]='one owner per edge'
    [edge_owner_range]='edge owners belong to the set of parts'
    [missing_partition_vertices]='a partition relation requires global vertex indices'
    [partition_number]='the part number belongs to the set of parts'
    [undeclared_whole_vertices]='the whole vertex carrier is declared'
    [undeclared_whole_edges]='the whole edge carrier is declared'
    [inconsistent_vertex_extent]='one vertex carrier has one extent'
    [inconsistent_edge_extent]='one edge carrier has one extent'
    [boundary_transpose]='a graph with an edge without a head has no transpose'
    [boundary_reverse]='a graph with an edge without a head has no transpose'
    [zero_fibre_index]='fibre index is outside its extent'
    [outside_fibre_index]='fibre index is outside its extent'
    [empty_fibre_index]='fibre index is outside its extent'
    [hierarchy_shared_extension]='a hierarchy is extended by its sole owner'
    [hierarchy_released_twin]="this storage's hierarchy has been released"
    [hierarchy_empty_index]='the index names a node this storage owns'
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
