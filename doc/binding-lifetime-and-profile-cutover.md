# Binding lifetime, then the graph_profile cutover
From `5627cac`. `src/fractal_graph.f90` unchanged throughout.
```
d88115b  gate     binding lifetime measured, then fixed
2a363b5  atomic   graph_profile + 34 constructor sites + 11 allowlist grants
```

## 1. Pointer-stability result
**The old binding was unsafe, and silently so.** Rows held the object as an `allocatable`
component; `bind_*` grew the row array by `move_alloc`, relocating it. Measured before the fix:
```
A before growth : same_as = T, name = R1
B after  growth : same_as = F            <- silently wrong
  name through the same pointer          <- SIGSEGV
```
`graph_profile` keeps exactly such a pointer (`this % tails => r`) for its view's lifetime. Had
the cutover gone first this would have been a corrupted profile, not a compile error.

## 2. Binding lifetime law
> A borrowed pointer stays valid for the life of the binding, across any number of later bindings.

A property of the **storage**, not of graph. A graph mutates freely under stable identity; an
object that lends pointers may impose stricter conditions than the graph does. Held on every run
by `test/graph-relational/lifetime.f90` — 12 assertions over all seven §2 cases: borrow before
growth, relation-storage growth, independent set-storage growth, intrinsic assignment, overwrite,
the `select type` narrow-then-grow pattern, and the call chain.

## 3. TARGET requirement
**None — and the fix removed the one that existed.** A row now holds a pointer to the object, so
the pointer returned does not point into the binding. `set_for` / `relation_for` lost their
`target` attribute, and no caller needs `TARGET` on a binding. Exercised by case G with the
binding declared without it.

## 4. Chosen storage fix
**Rows hold a `pointer` to an individually allocated object.** Growth copies rows; objects never
move. Obliges the binding to free its objects (`final`) and deep-copy on assignment
(`assignment(=)`), without which a copy would free objects the original still lends. Storage
machinery, not ontology.

| candidate | verdict |
|---|---|
| **individual allocation** | chosen — a mechanism, not a convention |
| build-then-borrow | rejected: convention only; the same silent corruption is one `bind_*` away |
| construct in one shot | rejected: same |
| copy into `graph_profile` | rejected: trades reference semantics away to avoid the question rather than answer it |

Diagnostics: none at `-Wall -pedantic -Wtarget-lifetime`. Runtime: valgrind 324 allocs, 324
frees, no leaks, 0 errors.

## 5. graph_profile signatures
```fortran
! before
ordinary_graph_view(g, tail_at, head_at)          class(relational_graph), target
directed_adjacency_view(g, selector)              class(relational_graph), target
! after
ordinary_graph_view(g, binding, tail_at, head_at) type(graph), type(relational_binding)
directed_adjacency_view(g, binding, selector)     type(graph), type(relational_binding)
```
Semantics unchanged: same three refusals, same fibre reading, same borrowed binary relations.
Neither view stores `g` — they never did — so nothing was added to keep (§8).

## 6. Call sites changed
34 constructions across 16 test files, plus 13 Makefiles and 5 tower import gates.
`test/support/relational_fixture.f90` reads a `relational_graph` a suite already builds and
answers the same structure as a graph plus its binding, so no site hand-builds one. It returns
**pointers** into storage it owns, so a site may reuse one variable and earlier views keep
working. Scaffolding: under `test/`, never in libufvm, dies with the capability passes.

The tower import gates refused the new imports. Eleven levels across five towers are granted
`fractal_graph`, `graph_relational_view`, `relational_fixture` — one level at a time, recorded in
the gate files, which is what those allowlists are for.

## 7. graph_algorithms: not changed
**Audit result, not an omission.** It declares `type(directed_adjacency_view), intent(in)` at four
sites and never constructs one, never names `relational_graph`, never uses `graph_structure`.
My previous report listed it in the cascade; that counted a *source* dependency as a
*constructor* dependency. §10 says do not edit it, so it is untouched.

## 8. Suites
Every suite matches the gate commit. `graph-benchmark` remains red on the pre-existing missing
`class_graph_support.mod` — baseline debt, neither caused nor repaired here. New this pass:
`test/graph-relational/lifetime` (12 assertions).

## 9. Fresh graph_structure surface
Measured after the cutover; do not reuse the old 33/163/65/59/184.

| | before | after |
|---|---:|---:|
| **`use graph_structure` in `src/`** | 1 | **0** |
| files naming `graph_structure` | 33 | 34 |
| `relational_graph` hits | 163 | 165 |
| construction sites | 65 | 65 |
| declarations | 59 | 58 |
| `held_set` / `held_relation` uses | 184 | 184 |

The counts rose by one file because the fixture reads `relational_graph` by design. The number
that moved is the first: **no production module names `graph_structure` any more.** The only
`src/` file mentioning it is itself.

## 10. Next deletion unit
All 34 remaining files are tests plus the fixture. Retarget by capability, in this order — the
first two are where `relational_graph` is only a container and the swap is mechanical:

1. **construction / validity** — `graph-structure`, `graph-ordinary` refusals: `relational_valid`
   replaces `create_graph`'s refusals.
2. **access / identity** — `num_member_sets`, `member_set_at`, `relation_at`, `holds_set` sites.
3. **relation algebra, field calculus, towers** — the levels that hold a graph only to pass it on.
4. Delete `test/support/relational_fixture.f90` with the last suite that needs it, then
   `graph_structure.f90`, `relational_graph`, `held_set`, `held_relation`, `create_graph`, and run
   the §14 audit.
