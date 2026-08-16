# First production migration
Four commits on `tower-graph-as-sets-relations`, from `a1478ff`. Each builds independently
(verified by checking out each and running `build.sh`).

```
b77e9ee  A  promote fractal_graph to src; tests import it
687678e  B  one ontology; obsolete two-graph doctrine removed
8e4f4fe  C  retire computational_graph; epistemic view; graph-state retargeted
5dfeb69  D  relational_graph migration map only
```

## 1. Production fractal_graph LOC
`src/fractal_graph.f90` — **67 code lines**, 49 comment, 54 blank. Byte-identical to the tested
kernel at `a1478ff` (verified by `diff` against `a1478ff:test/fractal-graph/fractal_graph.f90`).
Imports `graph_identity` and nothing else.

## 2. graph_state LOC before/after
| | code lines |
|---|---:|
| `src/graph_state.f90` before | 161 |
| after | deleted |
| `src/graph_epistemic_view.f90` (replacement) | 63 |

Net −98 production code lines, and one fewer type in the tree.

## 3. Public symbols deleted
From `graph_state`: `computational_graph`, `GRAPH_STATE_VOID`, `GRAPH_STATE_DATA`,
`GRAPH_STATE_OPERATOR`, `GRAPH_STATE_REALIZED`, `state_name`, `void_graph`, `data_graph`,
`operator_graph`, `realized_graph`, and the bindings `state`, `structure`, `declare`, `id`,
`same_as`, `name`. From `test/fractal-kernel`: a second `module fractal_graph` exporting
`graph_arena`, `graph`, `branch_spec`.

**One capability is deleted, not refactored.** The old seats were `class(*)` and could hold any
type. The new branches hold graphs. Anything that wanted an arbitrary payload in a seat must now
bind it externally through a map keyed on identity — §5 of the migration, stated so no later
reader assumes it survived.

Also deleted: the `host => relational_graph` member. A graph does not ride on a graph of another
kind (§6).

## 4. Public symbols retained as view vocabulary
In `graph_epistemic_view` — six free functions over `type(graph)`, no type, no state:
`has_data`, `has_operator`, `epistemic_defined`, `epistemic_name`, `data_of`, `residual_of`.

`void`, `data`, `operator`, `realized` survive **only as strings returned by `epistemic_name`**,
for the four combinations where both branches are UNKNOWN or KNOWN. They are not constructors and
not constants. The five combinations involving NULL are outside the reading's domain:
`epistemic_defined` answers so, and `epistemic_name` refuses rather than inventing a fifth name
(§4 — the 3×3 space is not forced back into the old 2×2).

## 5. Test results
Baseline at `a1478ff`: 27 of 28 suites green. `graph-benchmark` was **already red** — a missing
`class_graph_support.mod` that predates this work.

| after | suites green | red |
|---|---:|---|
| A | 26 | graph-benchmark (pre-existing) |
| B | 26 | same |
| C | 26 | same |
| D | 26 | same |

26 rather than 27 because `test/fractal-kernel` is deleted. Every other suite's result is
byte-identical to its baseline at each step (`diff` of the pass/fail lists).

- `test/fractal-graph`: 86 PASS, 0 FAIL — 9 language fixtures, 6 immutability candidates,
  45 assertions, 18 mutation experiments, 8 refusals. All now compile against `src/`.
- `test/graph-state` retargeted: 23 assertions, 3 refusals.

## 6. Dependency direction
```
graph_identity
     ^
fractal_graph            src/, 67 lines, imports graph_identity only
     ^
graph_epistemic_view     src/, first production view
     ^
graph_carrier, graph_relation, graph_structure, ...   unmigrated
```

`fractal_graph` imports `graph_identity` and nothing else — verified: it names no
`graph_relation`, `graph_carrier`, `graph_structure`, `graph_state`, field, operator or minimizer.

One asymmetry, deliberate: the spike's other views (`graph_views.f90` — `attribute_map`,
`relation_view`, `residual`, `compile_csr`, `csr`) stay in `test/fractal-graph/`. This pass
promotes only the kernel, so `graph_epistemic_view` is the first production view while the rest
are still test-local. That is scope, not an unfinished edge.

## 7. relational_graph migration map
`doc/relational-graph-migration-map.md` (commit D). Result: **nothing in `relational_graph` is
ontology.** `held_set`/`held_relation` are container storage forced by Fortran's one-dynamic-type
arrays; `sets(:)`/`relations(:)` are sequences the kernel already encodes with two branches;
`num_member_sets`, `member_set_at`, `num_relations`, `relation_at`, `holds_set` are traversals
over that sequence or indexing into its compiled representation; `label` is an attribute for an
external map; the signature-validity law is a property of the view, checked when the view is
built. It collapses to a view plus a compiled representation, not a second graph type.

## 8. Next safest root to delete
`relational_graph` — but not yet, and the measurement says why: `graph_structure` is named by 36
files, of which only **`graph_profile.f90`** consumes it in production, through `relation_at`
alone (`:162`, `:170`, `:478`). The other 34 are test call sites.

The next commit-sized unit is a sequence view over `type(graph)` giving `num_relations` /
`relation_at` semantics, in the `graph_views` line rather than in the kernel. Then retarget
`graph_profile.f90` (three call sites), then the suites, then delete.

## Acceptance
- `grep -R "type.*computational_graph" src test` → **empty**.
- Exactly one fundamental `type :: graph`, in `src/fractal_graph.f90`.
- Two abstract types still named `graph` (`graph_grammar.f90:228`, `interface_graph.f90:95`) each
  carry an explicit **MIGRATION DEBT** banner and are not presented as ontology.
- No performance change: no hot path was made to pointer-walk a graph; CSR and contiguous arrays
  remain compiled representations (§9).
- Level numbering untouched (§3). The releveling proposal and naming audit are marked superseded
  rather than acted on.
