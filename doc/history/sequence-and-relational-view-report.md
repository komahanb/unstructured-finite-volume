# Sequence and relational views
Commits A and B of the relational_graph retirement, from `f04e980`. Both **additive**: no consumer
touched, `graph_structure` untouched, `src/fractal_graph.f90` unchanged.
```
8a10381  A  graph_sequence_view + suite
f4791eb  B  graph_relational_view + suite, alongside graph_structure
```
**C onward is not in these commits.** See item 10 — the measurement changed the plan.

## 1. Sequence representation law
A sequence is represented by a **branch**, not a graph:
```
branch    NULL            the empty sequence
          KNOWN -> cell   a nonempty sequence
          UNKNOWN         the sequence is not known
cell      branch(1) = KNOWN -> element      always KNOWN
          branch(2) = the rest, again a sequence branch
```
`[a,b,c]` is `(a,(b,(c,NULL)))`. The empty sequence has no cell and needs no graph.

Two failures kept apart: a cell whose `branch(1)` is not KNOWN is **malformed**, always refused; a
spine reaching UNKNOWN is well formed and **undetermined**, refused only where the answer depends
on the unknown part. Both API forms are built and compared on every run by
`test/graph-sequence/run.sh`:

| form | module code lines | lines per call site | empty sequence |
|---|---:|---:|---|
| **branch** (shipped) | 19 | 1 | natural: NULL |
| graph (rejected) | 22 | 5 | needs a guard, and still needs a branch to decide |

## 2. Sequence API
`src/graph_sequence_view.f90` — four public free functions, no type:
```fortran
sequence_defined(b)        ! extent known: the spine reaches NULL
sequence_size(b)           ! refuses on an unknown extent
sequence_element(b, k)     ! type(graph), pointer; only k cells traversed
sequence_contains(b, g)    ! by identity, one traversal
```
`sequence_contains` is in the view rather than derived outside because a caller looping
`sequence_element` would be O(n²).

## 3. Relational-view API
`src/graph_relational_view.f90`. `branch(1)` = member-set sequence, `branch(2)` = relation
sequence. Counting and indexing delegate to the sequence view; nothing here traverses a spine.
```fortran
num_member_sets(g)              num_relations(g)
member_set_at(g, b, k)          relation_at(g, b, k)
holds_set(g, b, s)              relational_valid(g, b)
```
`b` is a `relational_binding`: element graph → the legacy `member_set` / `relation` it denotes. It
**owns** its objects — `graph_profile` stores the pointer `relation_at` returns for its view's
lifetime, so a borrowing binding would dangle. This is where `held_set` / `held_relation` belong:
storage, never ontology. `relational_valid` **answers**, it does not refuse: the old `create_graph`
refused because it was a constructor, and a view over a graph it did not construct reports. That
draws §7's distinction structurally — malformed sequence → refusal from the sequence view;
relationally invalid → `.false.` here.

## 4. graph_profile call-site migration
**Not done — inspected, not edited.** `create_view` (`:162`, `:170`) and `create_adjacency_view`
(`:478`) each do `r => g % relation_at(k)`, `select type` to `binary_relation`, then store the
pointer in the view. The replacement needs a `binding` argument, so `graph_profile`'s public
signatures change. That is the cascade in item 10.

## 5. Test groups migrated
**None.** Two new suites, neither replacing anything: `test/graph-sequence` (15 assertions,
9 refusals), `test/graph-relational` (18 assertions).

## 6. graph_structure LOC deleted
**0 — not yet deleted.** `src/graph_structure.f90` is 344 lines, untouched. Remaining surface:

| | count |
|---|---:|
| files naming `graph_structure` | 33 |
| `relational_graph` hits | 163 |
| `relational_graph(...)` construction sites | 65 |
| `type/class(relational_graph)` declarations | 59 |
| `held_set` / `held_relation` uses | 184 |

## 7. New view LOC added
| module | code lines |
|---|---:|
| `src/graph_sequence_view.f90` | 100 |
| `src/graph_relational_view.f90` | 147 |
| `src/fractal_graph.f90` | **67, unchanged** (verified by `git diff` at each commit) |

## 8. Performance and complexity
`sequence_size` O(n), `sequence_element(k)` O(k), `sequence_contains` O(n), `holds_set` O(m+n),
`relational_valid` O(|P|·arity·(m+n)). This is the semantic view. No hot loop was made to
pointer-walk a graph; CSR and contiguous arrays remain compiled representations. No production path
changed, so none regressed.

## 9. grep audit
Deletion has not happened, so §14's audit is expected non-empty:
```
relational_graph          163 hits / 33 files
use graph_structure        32 hits
held_set / held_relation  184 hits
src/fractal_graph.f90     unchanged
```
All suites match their baseline at every commit; the two new suites are additions.
`graph-benchmark` remains red for the pre-existing missing `class_graph_support.mod` — baseline
debt, neither caused nor repaired here.

## 10. Next root to map — and a correction
**Correction to the previous report.** It said `relational_graph` had *one* production consumer.
That is wrong: `src/graph_algorithms.f90` uses `graph_profile`, so giving `graph_profile` a binding
argument cascades to it and to **34** `create_view` / `create_adjacency_view` call sites across 16
test files. The one-consumer figure counted direct `relation_at` users only.

That cascade is why C was not started. A half-retargeted `graph_profile` leaves a tree that does not
build, and every commit must build independently.

Remaining order, unchanged but now correctly sized:

1. `graph_profile` — add the binding argument; retarget its 3 `relation_at` sites and both
   constructors. The same commit must carry `src/graph_algorithms.f90` and the 34 view call sites,
   or the tree does not build. Largest single unit in the migration.
2. Test groups by capability: construction/validity, access/identity, relation algebra, profile,
   algorithms, field calculus, towers.
3. Delete `graph_structure.f90`, `relational_graph`, `held_set`, `held_relation`, `create_graph`;
   run the §14 audit.

Then map `interface_graph.f90`'s abstract `graph` — an interface rather than a container, so it
migrates differently.
