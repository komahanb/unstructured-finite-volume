# relational_graph migration map
Analysis only. Nothing in `src/` is deleted or edited by the commit that adds this file.

`graph_structure.f90` declares `relational_graph`, the second graph root. It is now migration
debt: `(S,P)` is a **view** of `G=(B1,B2)`, not a kind of graph (AGENTS.md, "The graph ontology").
This maps each member to what it becomes, using one question per item:

> Is this **ontology**, a **sequence representation**, a **relation view**, or **container
> storage**?

## Field and procedure map

| member | category | becomes |
|---|---|---|
| `held_set` | container storage | Deleted. A wrapper existing only because a Fortran array carries one dynamic type. A sequence of graphs needs no wrapper: `branch(1)` is the element, `branch(2)` the tail. |
| `held_relation` | container storage | Deleted, same reason. |
| `sets(:)` | sequence representation | A cons sequence of graphs, each read as a set. `S_i` is reached by traversal, not by subscript. The contiguous array is the **compiled representation** of that sequence, produced when an algorithm needs `O(1)` indexing. |
| `relations(:)` | sequence representation | The same, for relations. |
| `num_member_sets` | relation view | A traversal that counts the sequence, or `O(1)` off the compiled representation. Not ontology: the kernel does not count. |
| `member_set_at(k)` | relation view | Traversal to the k-th element, or indexing the compiled representation. Returns `type(graph), pointer`, which the kernel already gives via `branch(i) % known()`. |
| `num_relations` | relation view | As `num_member_sets`. |
| `relation_at(k)` | relation view | As `member_set_at`. This is the only member with production consumers today: `graph_profile.f90:162, 170, 478`. |
| `holds_set(S)` | relation view | A traversal comparing `same_as` against each element. Already expressible over the kernel with no new type. |
| `identity`, `label` | ontology, already present | `identity` is the kernel's; `label` is an attribute and belongs in an external map keyed on identity, as `graph_views.f90` already does for numbers, symbols and indices. |
| `create_graph` refusals | relation view | The signature-validity law — every slot of every relation names one of the graph's own sets — is a **property of the view**, checked when the view is constructed. It is not a kernel law, and the kernel must not learn `arity` or `domain`. |
| `declare`, `id`, `same_as`, `name` | ontology, already present | Provided by the kernel; `name` is an attribute as above. |

## What this means

Nothing in `relational_graph` is ontology. Every member is one of:

- **storage** that exists because of a Fortran array limitation (`held_set`, `held_relation`),
- a **sequence** the kernel already encodes with two branches (`sets`, `relations`),
- a **view** over that sequence (`num_*`, `*_at`, `holds_set`, the refusals),
- or an **attribute** that belongs in an external map (`label`).

So `relational_graph` collapses to a view plus a compiled representation. It does not become a
second graph type, and it does not need one.

## Cost, measured not assumed

`graph_structure` is named by **36 files** — 2 in `src/` (`graph_profile.f90`, and itself; the
third, `graph_state.f90`, is gone as of `8e4f4fe`) and 34 in `test/`. Only `graph_profile.f90`
consumes it in production, and only through `relation_at`.

That is the reason this is a map and not a deletion: the ratio is one production consumer to
thirty-four test call sites, so the work is retargeting suites, not rewriting algorithms.

## Next safest root to delete

`graph_profile.f90`'s dependency on `relational_graph` is a single procedure family reached
through `relation_at`. The order that keeps every step small:

1. Add a sequence view over `type(graph)` giving `num_relations` / `relation_at` semantics, in the
   `graph_views` line rather than in the kernel.
2. Retarget `graph_profile.f90` to it — one production file, three call sites.
3. Retarget the 34 test files, suite by suite, each independently green.
4. Delete `relational_graph`, `held_set`, `held_relation`.

Step 1 is the next commit-sized unit. Do not start it inside the pass that deletes `graph_state`.

Not addressed here, and still debt: the two abstract types named `graph` in `graph_grammar.f90`
and `interface_graph.f90`, each carrying a migration-debt banner as of `687678e`. They are
interfaces rather than containers, so they migrate differently and should be mapped separately.
