# Set identity re-rooting, the inclusion map, and the partition law
From `0a7a226`. `src/fractal_graph.f90` unchanged; no existing production module touched.
```
8979aac  G0  subset provenance + inclusion provenance + partition law
0b9bf75  A   set representation family, set map, set view
bfed9bc  B   inclusion map, declared-subobject algorithm
```
Commits C and D were **not** made; §15 forbids them and §1 is why. Suites: `graph-set-view` (26), `graph-inclusion` (21), `graph-partition` (19) new; `graph-carrier`
49 → 58, `graph-binary` 43 → 47. Nothing else moved; valgrind clean on both new suites.
**Both precision corrections adopted.** `|G| = O(N_semantic)`, `|representation| = O(N_extent)`,
and **a bulk extension must not require O(N_extent) semantic graph objects merely to exist** — a
10⁹-member set is one integer; a genuine 10⁶-node task graph may hold 10⁶ graphs. And a subset
carries **two** things, extent and declared inclusion: the subtype goes, the inclusion does not.

## 1. Can subset_set be deleted before relation signatures move? **No.**
Two blockers, both measured, neither removable without an adapter.

**Inheritance.** `src/graph_forms.f90:52` is `type, abstract, extends(subset_set) :: form`, with
`polynomial_form` and `harmonic_form` extending it and assigning the parent component directly. A
form *is* a subset, so deleting the parent re-roots three production types.

**The two signature gates.** All 19 `subset_set` constructions in `src/` hand their result out as
`class(member_set)`, and it flows into exactly two places:
```
src/graph_relation.f90:174     type :: slot ; class(member_set), allocatable :: carrier
src/class_graph_field.f90:86   class(member_set), allocatable, private :: on
```
plus `relation % domain(k)`, which *returns* `class(member_set), allocatable`. That is **105
declarations across 29 modules**. Until a relation slot and a field domain can be a graph plus a
set map, a subset cannot stop being a `member_set` without a compatibility ontology — which §15
forbids. So this pass stops additive. Every §14 bucket says the same: nothing sits in **A** or
**E** alone — each site hands its subset on as **C**, a relation signature or field domain; **B**
is `graph_forms` and the eight propositions; **D**, the dynamic type test, is **0 in `src/`, 5 in
`test/`**, test-side and moving with the cutover. **Commit C is empty**: no consumer moves alone.

## 2–4. What `host` carried, and what inclusion derives
`host` carried **the ambient's identity**, and nothing else `roll` lacked. G0 measures it: two
subsets over `[2,5,6]`, one carved from `A = 1..8` and one from `B = 1..8`, produce inclusions
with **identical tuples** and different signatures — one names A, one names B. So the extension,
`(s,s)` for `s in S`, is recoverable from S's roster alone, while the **ambient rides the
signature** and no reading of the tuples recovers it. `graph_inclusion_map` therefore stores an
association and no table: member values are shared, so `inclusion_value(s) = s` needs no storage.

## 5. Extensional subset vs declared subobject
An **extensional subset** asks whether every member of S belongs to A — a question about two
extents. A **declared subobject** asks whether an inclusion path `S -> ... -> A` exists — a
question about what was declared. `test/graph-inclusion` builds `S = {2,5,6}`, declares it into `A = 1..8`,
and builds a second set `B = {2,5,6,7}` also containing every member of S: S is extensionally
inside both and a declared subobject of A alone. **An inclusion edge is never inferred from
containment.**

The old law is preserved verbatim: `S <= S`, `S <= A`, `S` is **not** `<= B`; `T c--> S c--> A`
gives `T <= A` and not `T <= B`. The same eight propositions now hold in two places —
`graph-carrier` over `subset_set`, `graph-inclusion` over graph + map — so when the old type goes,
the law demonstrably did not go with it.

## 6–7. No semantic identity in the representation; counted still O(1)
`set_representation` has five deferred procedures — `size`, `member`, `members`, `has`,
`local_index` — and nothing else: no token, no `declare`, no `id`, no `same_as`, no `name`. Two
graphs over 1..8 are extensionally equal and are two sets; two **empty** sets described by
*different representation kinds* likewise. The representations take no part, having nothing to
take part with. The counted representation is **one integer**, asserted through the production
path at `|cells| = 10^9` with both branches of the set graph `GRAPH_NULL` — the extension is
provably not in the graph, and membership is one comparison.

## 8. Did any new API lend a pointer? **No — and that is the result.**
Every answer is by value; `ambient_of` returns a `type(graph)` by value, and a copy carries the
token, so it is the same declared set. `relational_binding` needed individually allocated objects
behind pointers, a `final`, a defined assignment and a refusal of assignment. `set_map` needs
none of it: plain `allocatable` rows, free `move_alloc` growth, and intrinsic assignment already
deep-copies. **Not lending is what removes the lifetime problem, not care.** §10's STOP rule
stands.

## 9. Partition / assembly, for nparts > 1
Current public interfaces only, no merge operation invented. Four things, four answers:

| | law | nparts 1 | 2 | 3 |
|---|---|---|---|---|
| **B** vertex data | `sum_k A_k(P_k(D)) = D` | ✓ | **✓** | **✓** |
| **C** edge data | the same sum law | — | **✓** | **✓** |
| **A** structure | coverage: assembled edges are G's, and G's edges are all assembled | ✓ | ✓ | ✓ |
| **A** structure | *identity* | ✓ | ✗ | ✗ |
| **D** ownership | every cell owned once; cut cells borrowed | — | ✓ | ✓ |

**The data half is stronger than the maps assumed** — a genuine multi-part law, because the
assembler writes only what a part owns and leaves the rest at zero. **The structure half is
weaker than "section" suggested**: union over parts is not a public operation, so coverage is the
strongest statable law. The one-part case is a **restricted identity on the one-part subdomain**,
and the suite measures the restriction: at nparts = 2 one part alone assembles back strictly
fewer edges than G.

## 10. The next atomic migration boundary
**One cutover: relation slot + field domain + the `form` inheritance.** `slot % carrier` and
`field % on` become a graph plus a set map; `relation % domain(k)` answers a graph; `form` stops
extending `subset_set`. Only then does `subset_set` delete, and the 19 construction sites and 105
declarations move with it — one commit or none, since an adapter in the middle would be the
compatibility ontology §15 exists to prevent.

Two debts stay owed and unstarted: the transpose borrower lifetime gate, and an
assemble-from-all-parts operation that would let structure be stated as strongly as data.
