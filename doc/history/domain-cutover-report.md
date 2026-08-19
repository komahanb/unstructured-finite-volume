# The domain cutover: src crossed, tests not

Two branches. `tower-graph-as-sets-relations` stays **green at 32/33** (`3a87ab5`).
The cutover is `tower-domain-cutover` / **`e2f32ee`** — 36 files, +1302/−692, in the
ruled order. **The library builds clean and the three foundational suites pass 77
propositions. The 132 test files that construct domains are NOT migrated**, which
is why it is a branch and not a merge.

## 1. Did graph acquire any set/map/label storage? **No.**

`grep -nE "type\((set_map|label_map|inclusion_map)\)" src/class_graph.f90
src/fractal_graph.f90` is empty of stored components. `stored_graph` keeps
`type(set_graph) :: vset, eset` — identity beside the `nv`/`ne` it already had,
O(1) — and the extension is reconstructible by any caller as
`counted_set_representation(g % num_vertices())`.

**2. Where does domain label live?** `graph_label_map` (G1, `b9494ca`), keyed on
token by value. `stored_graph % name_carriers(labels)` binds `'vertices'`/`'edges'`
into the *caller's* map, preserving `test/graph-carrier:400`.

## 3–5. Field, and the dependencies that became explicit

**3.** A field retains `type(set_graph) :: on` + `integer :: nentries`, private,
plus values. **4. Duplication, same probe, before and after:**

| | values | batch | domain excess |
|---|---|---|---|
| before, 40 fields on a 200k **listed** domain | 61.0 MB | 91.5 MB | **28.7 MB** |
| after | 61.0 MB | 62.6 MB | **~0** (1.5 MB, what a *counted* domain always cost) |

**5.** Explicit map dependencies, and only these:

- **partition_data / assemble_data** take `sets, labels, inclusions` — the only
  consumers that genuinely ask membership, enumeration, position, embedding and
  label. `is_subobject_of` became `declared_subobject(dom, X, inclusions)`.
- **the twelve carved-domain accessors** take all three and bind into them.
- **relation and CSR construction** take `sets`, read it once, store nothing.
- **`class_robin_condition % faces`** takes them because its answer *outlives the
  call*; `measures_of` and `diffusion_statement` hold them as **locals**, because
  the set they carve is born and dies inside one call. A local temporary is not a
  hidden dependency.

Everything else needed none, because §3's escape held: `operation % domain(g,
dom, nentries)` carries the count beside the identity, and that removes the map
from the whole operation and field chain. **No set_map is threaded anywhere
merely to ask `size()`.**

**The stable/carved split** is what made this possible. `vertex_set`, `edge_set`,
`all_vertices`, `all_edges` answer a domain the graph already holds — identity
alone, no map. The other twelve mint a fresh set per call, as they always did.

## Carved creation is atomic

All twelve route through one private `carve`, and `graph_algorithms` repeats it
for `sources`/`sinks`:

    call members % declare()
    call sets       % bind(members, listed_set_representation(roll))
    call labels     % bind(members, label)
    call inclusions % include_in(members, ambient)

Bound in one place because the third is the one an author forgets: a missing
representation stops the program at the first query and a missing label answers
`''`, but a **missing inclusion answers false to the subobject question —
silently, and only on a real mesh**. Partition and assembly carve the same way
when transporting a subset, carrying the source domain's label across.

## 6–10. Relation and CSR

**6.** `type(set_graph), allocatable :: signature(:)` — an ordered array of graph
values. That **is** the ordered graph sequence; `graph_sequence_view`'s own header
says to compile to a contiguous representation for indexed access, and its
`sequence_element` hands back a *pointer* into cells the holder must keep alive —
the borrow `400a2c3` removed. Repeated domains repeat the identity. **7. slot is
gone** (remaining `grep` hits are message strings and an unrelated
`interface_graph` local).

**8.** `csr_relation` holds two `class(set_representation)` coordinate stores,
copied by value out of the map at construction, beside — and visibly apart from —
the semantic `signature(:)`. Identity answers `domain(k)`; coordinates answer
`local_index`. Neither answers the other's question. `set_map % extent_of`
returns a fresh allocatable, so this is a copy, not the borrowed accessor the map
header gates.

**9. No.** `image_view`, `preimage_view` and `has` reach `source_coords` /
`target_coords` directly — no map row scan, no graph traversal, no label lookup —
and all three remain `pure`. **10.** For N=200 000, nnz=1 200 000: CSR arrays
**10.681 MB**, coordinates **2 integers** (counted domains), ratio ~0. Were both
domains listed instead: 1.526 MB, ratio **0.143**.

## 11–12. Form, and subset_set

**11. No** — `type, abstract :: form`, holding a basis identity and a
`listed_set_representation`. The inheritance bought exactly one method,
`members()`, and charged for the whole carrier contract. `restrict` still
**mutates**, deliberately unchanged. **12. Not yet.** `subset_set` has no
remaining use in `src/` outside its own definition file, but the test suites
still construct it, so deleting it is part of the test migration, not separable
from it.

**13.** `member_set`/`counted_set`/`subset_set` survive only in
`src/graph_carrier.f90`; every other src module is clean. After the tests move,
their src surface is mechanically empty and `graph_carrier` can go whole.

**14. `src/fractal_graph.f90` is byte-identical** to `57d8c51`. The core learned
nothing about set, label, field, relation, CSR, inclusion or partition.

## The accepted PURE loss, enumerated (ruling §17)

Eleven procedures, all semantic accessors returning a graph or a field by value:

    graph_carrier_interface, vertex_set, edge_set          domain answers
    cell_volume, cell_center, face_area, face_delta,
    face_normal, face_center, face_weights                 mesh field accessors
    reduction % create                                     declares its home

Preserved where it matters, exactly as ruled: `csr_has`, `csr_num_tuples`,
`csr_tuples`, `view_has`, `view_num_tuples`, `view_tuples` are all still `pure`,
**because the coordinates are the relation's own**. The §8/§9 split is what keeps
the loss out of the numerical path.

## 15. Next atomic boundary

The test migration, which is now the *whole* remainder and is one commit with
`subset_set`'s deletion: 132 files across 33 suites — 7 towers at 11–17 files
each (their fixtures implement `graph_operation`, so each needs its domains
declared and its maps threaded), plus 27 unit suites. The pattern is mechanical
and fixed: declare a `type(graph)` where a `counted_set` was constructed, bind a
representation and a label beside it, carve where a `subset_set` was built, and
pass `sets` to relation construction. `test/graph-set-view`, `test/graph-inclusion`
and `test/graph-identity-map` are already written that way and pass on the branch
— they are the worked examples for the rest.
