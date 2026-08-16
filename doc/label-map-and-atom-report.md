# The label map, and the measured interior of the A/B/C/D atom

From `0f7e646`. **`b9494ca` (G1) is the only commit** — it lands
`graph_label_map` complete and green. A/B/C/D was then attempted and
**reverted**, not committed: the atom needs the whole repository to move at
once, and a half-migrated tree is worse than an unmigrated one. What the
attempt bought is measurement, below, including one architectural cost that
was not in the plan.

## G1 — done

**1. Did graph acquire any set/map/label storage? No.** Nothing was added to any
graph object. **2. Where does domain label live?** `src/graph_label_map.f90` —
`set graph -> character label`, keyed on `type(token)` by value, `bind` takes no
`TARGET`, and the file contains neither the word `pointer` nor `target`
(§19 carried forward from `400a2c3`).

The four §2 proofs are `test/graph-identity-map/labels.f90`, in the suite that
already runs valgrind: two graphs share one label and stay two sets; a copied
token retrieves the same label; the map answers after every binding graph is
deallocated; assignment deep-copies. Refusals: an unsigned token, naming twice.
There is deliberately **no lookup from label back to graph** — it would have to
answer for a string naming two objects. Naming is not addressing. Stated
decision (§2 left it open): an unbound graph's label is `''`, not a refusal.

**14. `src/fractal_graph.f90` is byte-identical** to `57d8c51`. Suites 32/33;
`graph-benchmark` was already broken at HEAD (`class_graph_support`, `a3817e3`).

## The unbudgeted cost: naming a domain `type(graph)` removes `pure`

`type(graph)` contains `type(graph_branch)`, which contains
`type(graph), pointer :: known_`. F2018 **C1594** forbids copying such a type
out of an `intent(in)` dummy inside a `pure` subprogram. So *every* accessor
that answers a domain or a field by value must drop `pure`. In the fraction of
`src/` the attempt reached, ten already had to:

    graph_carrier_interface, vertex_set, edge_set        (the domain answers)
    cell_volume, cell_center, face_area, face_delta,
    face_normal, face_center, face_weights               (mesh field accessors)

The mesh accessors are the surprise: they lose purity not because *they*
changed but because `type(field)` now transitively contains a pointer — which
spreads to any `pure` procedure returning a field by value, anywhere. It does
**not** reach the hot path, and that is §9 working as designed: CSR keeping
coordinate representations *by value* rather than reaching a graph is what keeps
`csr_has` and the fibre views `pure`. A real, permanent cost, priced nowhere in
the plan.

## §5/§6 — field storage, and why it no longer needs a map

**3. What exactly does field retain?** `type(set_graph) :: on` and
`integer :: nentries`, both private, plus values. No representation, no map, no
inclusion. `domain()` answers the graph by value; `num_entries()` answers the
snapshot; the five value-vector width checks use the snapshot.

**The snapshot changes no behaviour**: `class_graph_field.f90:161` did
`allocate(this % on, source=on)`, so a field has held a *copy* of its domain
since construction and `num_entries` never tracked later mutation. Freezing the
integer says out loud what copying said in private.

**4. Field domain duplication before/after.** Before, measured: 40 fields on one
200 000-member subset domain carry **28.7 MB** of duplicated extension (61.0 →
91.5 MB against 61.0 MB of values), matching the 30.5 MB predicted for K copies
of an N-integer roll. After: the field holds one graph and one integer, so the
excess is **0**. That follows from the type, not from a benchmark — the atom was
reverted, so it could not be re-measured against a built library.

**5. Which map dependencies became explicit?** Decisively, **almost none need
to be**, because §3's own escape applies — *S size → frozen snapshot where
valid* — and it is valid throughout the operation chain. Changing
`operation % domain(input_graph, dom)` to `domain(input_graph, dom, nentries)`
lets the count travel beside the identity, and that one change removes the map
from every operation, every field construction and every domain check in the
tower. Maps become explicit in exactly three places: **partition/assembly data
transforms** (§4), which genuinely ask M/E/L/O/N; **the carved-set accessors**
(below); and **relation/CSR construction** (§11), for membership validation.

`class_robin_condition` shows a fourth, smaller case: `faces()` takes the maps
because its answer *outlives the call*, while `measures_of()` carves and consumes
a set within one call and holds them as **locals** — a local temporary is not a
hidden dependency, since nothing needing it escapes.

## The grammar split nobody had named: stable vs carved

Measured in `class_graph.f90`: `all_vertices`/`all_edges` return a copy of the
stored carrier (**stable**, `:582`, `:698`), while `interior_*`, `boundary_*`,
`tagged_*`, `owned_*`, `borrowed_*`, `overlap_*` build a fresh `subset_set` per
call (**carved**, `:607` and eleven more) — so two calls have *always* been two
domains. **Zero tests assert `same_as` across two such calls**, so fresh-per-call
is preservable literally. That splits the 26 deferred interfaces cleanly:

| kind | count | signature |
|---|---|---|
| stable | 4 | `type(set_graph) function(this)` — identity only, no map |
| carved | 12 | `subroutine(this, [tag/part], sets, labels, inclusions, members)` |

A carved set must bind **three** things or it is unusable, in one place: a
missing representation stops the program at the first query and a missing label
answers `''`, but a **missing inclusion answers FALSE to `is_subobject_of` —
silently, and only on a real mesh** (assembler:219, partitioner:481). The
attempt routed all twelve through one private `carve` gate for that reason.

`test/graph-carrier/test.f90:400` asserts the carriers are named `'vertices'`
and `'edges'`. Preserved (§16) by `stored_graph % name_carriers(labels)`: the
graph tells the caller its own domain names, binding into the caller's map and
owning nothing.

## Not executed

**6.** Signature is still `type(slot)` over `class(member_set)`. **7.** slot
remains. **8/9/10.** CSR still holds slot-wrapped carrier copies and calls
`carrier % local_index` per query at `graph_binary_relation.f90:410, 428, 453`;
no coordinate representation was built, so no ratio was measured. **11.** `form
extends subset_set` still. **12.** subset_set remains. **13.** `member_set` 163
src / 396 test, `counted_set` 43 / 596, `subset_set` 52 / 289, `slot` 54 / 184.

**Why it is one atom, proven rather than asserted:** `class_graph % vertex_set()`
feeds *both* field domains and relation signatures. Retyping it to
`type(set_graph)` breaks relations in the same edit that fixes fields, so A and
B cannot be separated; `subset_set` deletion then requires C and D. Measured
surface: **134 files** need to move together — 213 field constructions (39 src /
174 test), 190 stable-accessor call sites, 92 carved-accessor call sites, 28
`graph_operation` implementations (9 src / 19 test fixtures), 4
`graph_transform`. The attempt rewrote 31 src files (+557/−277) and left one
compiling blocker at `graph_relation.f90` — the entrance to B, which is a
redesign rather than a mechanical edit.

**15. Next atomic boundary.** Unchanged in membership; now with its interior
mapped. Recommended order, tightest constraint first: (1) the 26 grammar
interfaces above — write them and let the compiler enumerate the rest;
(2) `class_graph`, `stored_graph` keeping `type(set_graph) :: vset, eset`
declared once at construction alongside `nv`/`ne`; (3) field; (4) relation,
slot, CSR coordinates; (5) form, then `subset_set`. One decision is still
yours and should be made before step 1: **whether the purity loss above is
acceptable**, since it is permanent and reaches every `pure` accessor that
returns a field by value.
