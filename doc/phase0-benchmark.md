# Phase 0 hot-traversal baseline (AGENTS.md §59, §66)

Recorded 2026-08-12, before any relational-refactor change. The
program is `test/graph-benchmark/bench.f90`; rerun it after each
migration step and compare against this table. AGENTS.md §66: the new
abstraction must not impose material regression on these paths
without a documented reason and an optimization plan.

## Setup

- Graph: structured 700 x 700 grid as a stored_graph — 490,000
  vertices, 978,600 edges (right + down neighbours).
- Machine: Intel Core i5-7500 @ 3.40 GHz, Linux 5.15.
- Compiler: GNU Fortran 15.2.0, `-std=f2023 -fcoarray=single -cpp -fPIC`.
- Two builds: the repo default (no `-O` flag) and `OPTIMIZE=yes` (`-O3`).
- Traversal sweeps visit every vertex three times; "ns/item" is per
  query for sweeps, per edge for operator applies, per vertex for
  construction/partition rows.

## Numbers

| act                              | default | -O3     | ns/item (default) | ns/item (-O3) |
|----------------------------------|---------|---------|-------------------|---------------|
| stored_graph construction        | 0.137 s | 0.046 s | 278.8             | 94.3          |
| incident_edges sweep (x3)        | 0.055 s | 0.034 s | 37.6              | 22.9          |
| adjacent_vertices sweep (x3)     | 0.050 s | 0.033 s | 34.2              | 22.7          |
| outgoing_edges sweep (x3)        | 0.048 s | 0.033 s | 32.3              | 22.3          |
| incoming_edges sweep (x3)        | 0.048 s | 0.033 s | 32.4              | 22.4          |
| divergence apply (x3)            | 0.163 s | 0.099 s | 55.4              | 33.7          |
| laplacian apply (x3)             | 0.195 s | 0.110 s | 66.4              | 37.5          |
| partition_graph, 4 parts         | 0.269 s | 0.127 s | 549.4             | 259.4         |
| carry + assemble field, 4 parts  | 0.294 s | 0.137 s | 600.4             | 279.1         |

## Phase 3 comparison (added 2026-08-12)

The CSR binary relation (`src/graph_binary_relation.f90`) on the same
topology — the tail relation V x E, image(v) = outgoing edges,
preimage(e) = tail vertex — against the old stored_graph queries:

| act                              | default | -O3     | ns/item (default) | ns/item (-O3) |
|----------------------------------|---------|---------|-------------------|---------------|
| csr_relation construction        | 0.082 s | 0.038 s | 167.2             | 78.1          |
| csr image sweep (x3)             | 0.068 s | 0.036 s | 46.3              | 24.4          |
| csr preimage sweep (x3)          | 0.138 s | 0.074 s | 47.0              | 25.4          |
| csr image_view sweep (x3)        | 0.033 s | 0.013 s | 22.1              | 9.0           |
| csr preimage_view sweep (x3)     | 0.067 s | 0.026 s | 22.7              | 8.9           |

Two tiers, and the split the Phase 0 notes could only conjecture is
now measured: the same fibre costs 24.4 ns owned (allocating `image`)
and 9.0 ns borrowed (`image_view`, a pointer slice) at -O3 — the
difference IS the allocation. The borrow beats the old
`outgoing_edges` (22.5 ns) by 2.5x; the allocating convenience sits
at parity with the old API; construction is faster than
`stored_graph`'s (78 vs 100 ns per vertex). The §66 gate holds with
room to spare, and hot Phase-4 topology should ride the views.

Complexity, stated honestly: every fibre first pays the carrier's
`local_index(member)`, so T_image(a) = T_local_index(a) + O(deg a).
Counted carriers answer local_index in O(1) — the V/E mesh path —
so the promise collapses to O(deg) exactly where it matters.

## Reading

- Every traversal query today returns a freshly allocated index
  array, so the ~22 ns/query at -O3 plausibly carries substantial
  allocation overhead — this benchmark does not isolate that split.
  A future CSR-slice view (AGENTS.md §33) is expected to improve on
  it; the requirement is that hot queries must not materially
  regress above these times.
- Divergence and laplacian walk the edge list directly; their per-edge
  cost is the honest hot-loop reference for the Level-6 migration.
- The partition rows time the whole four-part construction including
  part-graph assembly; they bound the Phase 7 transport rework.
