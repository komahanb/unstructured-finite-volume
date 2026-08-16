# Mapping member set and relation onto the fractal graph
Analysis pass, from `3b8e8d8`. **Nothing migrated, nothing renamed**, `src/` untouched.
```
da21e6a  M1  member-set-fractal-map.md      + set.f90, scale.f90
43950f6  M2  relation-fractal-map.md        + relation.f90
3d42056  M3  naming-and-transformation-algebra.md
```
Each builds independently, verified in a clean worktree. New suite `test/fractal-map`, 37
assertions, run on every pass; every other suite unchanged.

## The finding both maps turn on
**Members and tuple components are already integers, not graphs.** `counted_set` holds one
integer; `stored_relation` holds `integer :: entry(:,:)`; `csr_relation` holds four integer
arrays. The semantic/representation split this pass was asked to *design* **already exists in the
data** — unnamed, with the identity token on the wrong side of it. So the migration ahead is a
**re-rooting of identity**, not a redesign of storage, and §18's exascale requirement is already
met by the storage that exists.

## The measurement (`test/fractal-map/scale.f90`)
56 bytes per `type(graph)`, measured, not assumed.

| | design A · extensional graph | design B · semantic graph + representation |
|---|---:|---:|
| `\|S\| = 10^9` | **112 GB** of graph objects | **60 bytes** |
| `\|R\| = 10^10` tuples, arity 2 | **1680 GB** of graph objects | **88 GB** of flat integer arrays |

Design A is rejected **for extensions** and kept **for signatures and short sequences**, where
`graph_sequence_view` already uses it. A 10⁹-member extension is not control-plane data.

## The 18 answers
**1. Does `member_set` survive as ontology?** No. Its five deferred questions are one *storage*
contract; its identity block belongs to a graph. It survives as `graph_set_representation`.

**2. Is a finite set a view, a representation, or both?** Both: a **graph** answers *which set?*,
a **representation** answers *how are its members stored here?* Neither alone is a set.

**3. Can `counted_set` remain O(1)?** Yes — demonstrated. The prototype declares a 10⁹-member set,
allocates no member object, and asserts both branches of the set graph are NULL.

**4. Does `subset_set` need a subtype?** No — but it does need a second *representation*.
`inclusion_of` builds tuples `(s,s)`, the **same value on both sides**, while indexing CSR rows
through `S % local_index`, the *position*: measured, member 5 stands 2nd in S and 5th in A. Two
coordinate systems, not three. `S ⊆ A` is a property between two set graphs; `ambient`/`host`
become a map, `is_subobject_of` an algorithm, and `roll` becomes
`graph_set_listed_representation`.

**5. What becomes `local_index`?** The representation's numbering function, **not identity**.
**6. Does `relation` survive as ontology?** No — same shape as (1). It becomes
`graph_relation_representation`, with `arity`/`domain(k)` lifted into a view.

**7. Where does relation identity live?** On the relation's `type(graph)`. Demonstrated: one
relation graph bound to a table **and** a CSR of the same pairs — both answer the same extension,
and `same_as` between them is false. Representation identity is storage identity, a different
question. **8. Where does its tuple extension live?** In the representation, as integers —
already there.

**9. Does `binary_relation` need a type?** No. `source`/`target` already read
`domain(1)`/`domain(2)` verbatim in the production file that defines them; `image`/`preimage` are
algorithms over CSR. A *type* encoding "arity is 2" duplicates what the signature says.

**10. Does `materialized()` survive?** No. Its three inhabitants split **representation from
view**, not relation from relation. A boolean fixed by an object's *kind* is evidence the kind is
missing from the type system. **11. Is CSR a relation or a representation?** A representation —
with the table, one of two storages for one relation.

**12. What equality does the transpose involution preserve?** **Extension only, today.**
`transpose_of` mints a fresh token, so `T(T(R))` is a *third* identity; asserted both ways.

**13. Can double transpose return the original semantic object?** **Yes** — if transpose is a role
permutation over a relation graph rather than a constructor. The prototype builds it: two
transposes compose to the identity permutation and the pair denotes the same R. The classification
is involution **on views** — R has no transpose until a view reads its signature in an order.

**14. What is the exact law of partition/assembly?** **Class C, section/retraction** — by measured
evidence. `test/graph-contract` proves `assemble(partition(G)) == G` only at **nparts = 1**, and
calls that "the weakest case" itself. For nparts > 1 the partitioner returns *one* part, so a
single-part round trip cannot be the identity; what is proven there is **ownership** — every cell
owned exactly once, three rules, four part counts. Promoting the pair to class B needs a test
that does not exist: partition into k, assemble from all k, compare against G.

**15. What stores ownership / local-global / halo?** Today: four integer arrays and two scalars
**inside `stored_graph`** — ρ living inside the object it describes. By role it is four things
(locality map, ownership map, a view over that map, the cut's state); the partition and assembly
*algorithms* are two more, and a communication representation is a sixth this repository does not
have yet and must not conflate with ownership.

**16. What naming algebra?** `graph_<mathematical noun>_<role>`, roles closed to
view / map / representation / algorithm / algebra / calculus. One law follows: **only maps and
representations may be O(n).** Views and graphs are the semantic control plane and stay small.

**17. Does either proposal require O(n) semantic graph objects?** No — that is why design A is
rejected for extensions. Signatures stay O(k) with k small.

**18. What should migrate next?** `subset_set` → set graph + subset relation: one construction
gate, five overrides, and every consumer already asks `is_subobject_of` rather than reading
`host`. Then `slot`, `materialized()`, `binary_relation`, and partition/assembly last — gated on
the missing assemble-from-all-parts test.

## What the prototypes demonstrate rather than assert
- **Two distinct empty sets**: `set_equivalent` true, `same_as` false — the empty set of cells is
  not the empty set of faces.
- **Design A, built** at n = 4, so its cost is measured; **A × A, A × B, A × B × C** through one
  mechanism, seven cells for all three; a counted and a listed extent, positions disagreeing.
- **The tuple trap**: two independently built holders of `(a,b)` are *not* the same graph while
  being the same tuple; `(b,a)` is a different tuple. Holder identity must never be tuple equality.
- **One relation, two representations**, agreeing on extension and differing in object identity;
  **transpose twice**, in production (extension only) and as a role permutation (identity).

## Four things the maps refuse to say
- **Transpose is not adjointness.** One permutes signature slots; the other is an identity between
  inner products. This repository proves both, in different suites.
- **Dual is not transpose.** They share a law class — involution — and nothing else. **Nothing
  implements dual**; it exists only in AGENTS.md prose, and stays a stated law for a future view.
- **Assembly is not the inverse of partition** (answer 14), and **refine/coarsen is
  unclassified** — class D is the hypothesis, and no test settles it.

## Stopped, per §26
`member_set`, `relation` and partition/assembly are untouched. The next pass should be chosen from
the measured consequences above and should carry the missing partition round-trip test with it — a
law nobody has checked is a poor thing to migrate on.
