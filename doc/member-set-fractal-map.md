# Member set, mapped onto the fractal graph
Analysis only. Nothing in `src/` is changed; the prototypes are `test/fractal-map/set.f90`
and `test/fractal-map/scale.f90`, run on every pass.

## 0. The finding that reframes everything
**A member of a member set is already an integer, not a graph.**

```fortran
pure integer function set_member_interface(this, local_index)   ! answers an integer
pure logical function set_has_interface(this, member)           ! takes an integer
integer, private :: n                                           ! counted_set, entire
```

The semantic/representation split this map is asked to design **already exists in the data**.
It is not named as a split, and the identity is rooted in the wrong place: `member_set` carries
the `token`, so the thing that answers *which set?* and the thing that answers *how are its
members stored?* are one Fortran type. The migration is a **re-rooting of identity**, not a
data-structure redesign, and §18's scale requirement is already satisfied by the storage that
exists.

## 1. Two designs, both built
`test/fractal-map/set.f90` builds design A at n = 4 so its cost is measured, not hypothesised.

| | design A · extensional graph | design B · semantic graph + extent |
|---|---|---|
| what the graph holds | `branch(1)` = member sequence | nothing; both branches NULL |
| objects per member | 2 (one element, one cell) | 0 |
| `|S| = 10^9` | **112 GB** of graph objects | **60 bytes** |
| membership | O(n) traversal | one comparison |
| device-resident | no — pointer chasing | yes — the extent is integers |

56 bytes per `type(graph)`, measured. Design A is **rejected for extensions**: it forces
semantic allocation proportional to a mesh's cell count merely for the set to exist. It is
retained for **small** sequences, where it is the right tool — signatures, member-set lists,
relation lists (`graph_sequence_view` already does exactly this).

Design A is not wrong; it is *the semantic control plane*, and a 10⁹-member extension is not
control-plane data.

## 2. Every symbol in graph_carrier, classified
| symbol | role | disposition |
|---|---|---|
| `member_set` (abstract type) | **representation contract** — its five deferred questions are all about storage | does **not** survive as ontology; becomes `graph_set_representation` |
| `identity` (token) | **graph identity** | moves to the set's `type(graph)`; leaves the representation |
| `label` | metadata | moves with identity; already "no part of the mathematics" |
| `declare` / `id` / `same_as` | **graph identity** | the kernel's, unchanged |
| `size` / `member` / `members` / `has` / `local_index` | **representation** | stay together — they are one storage contract |
| `name` | metadata | external map, keyed on identity |
| `counted_set` | **representation** | `graph_set_counted_representation` |
| `n` | representation state | unchanged, O(1) |
| `subset_set` | **redundant as a type** | a set graph + a subset relation |
| `host` | **map** (subset → ambient) | external, or read off the subset relation |
| `roll` | **representation** | `graph_set_listed_representation` — an explicit-index extent |
| `ambient` | **map** | the same association as `host`, read forward |
| `is_subobject_of` | **algorithm** over the subset relation | derived, not stored |

Nothing is unclassified. The one type that survives as a *type* is the representation family;
`member_set`'s identity block is the part that leaves.

## 3. Counted set stays O(1) — demonstrated
`test/fractal-map/set.f90` §3 declares a set of 10⁹ members, binds a counted extent, and
allocates **no member object**. `has(999999999)` is one comparison; `member` and `local_index`
are the identity map. The set graph's branches are both NULL, asserted — the extension is
provably not in the graph.

`size`, `member`, `has`, `local_index` belong to the representation, all four. They are the
questions "how are the members stored and enumerated here?", and a bitset, a device array or a
run-length extent would answer them differently for the same set.

## 4. Subset needs no subtype — and the reason is two coordinate systems
Current `subset_set` conflates them; the prototype separates them, and both are measured against
production for `S = {2,5,6}` inside `A = 1..8`:

```
                  member VALUE      local POSITION
      in S              5                 2
      in A              5                 5
```

Exactly **two** coordinate systems, not three: the value is shared, the position is private to
each representation. `inclusion_of(S)` is the proof in production code — it builds tuples
`(s % member(k), s % member(k))`, **the same value on both sides**, while indexing its CSR rows
through `S % local_index`, which is the *position*. So this repository is already in §5's case
**A**: a subset reuses the ambient's member values, and inclusion needs no translation map.

The prototype builds a **listed extent** beside the counted one and asserts the enumeration law
`member(local_index(v)) = v` holds inside *each*, with the two positions disagreeing. That is the
whole of what `roll` does. The claim is therefore precise: a subset needs **no subtype**, and it
does need a second **representation**, `graph_set_listed_representation` — a type deletion plus a
representation, not a type deletion alone.

Therefore:
```
S and A are both set graphs
S ⊆ A is a property between them, decided from the two representations
```
No third type, no inheritance, no `ambient` component. `is_subobject_of` becomes transitive
closure over the subset relation — an algorithm, and the place where a chain of subsets is
walked once rather than stored in every link.

**`local_index` becomes** the representation's numbering function, and it is *not* identity.
The integer 7 names a member within one representation of one set; two sets may both hold a
member numbered 7 and mean nothing by it.

## 5. Identity is not extensional equality
Preserved exactly. Two independently declared empty sets are:
```
set_equivalent(cells, faces) = .true.      equal members
cells % same_as(faces)       = .false.     different sets
```
Both asserted. `set_equivalent` is defined **separately**, in the prototype, over extents — it
is a view's question, never the kernel's. The empty set of cells is not the empty set of faces,
and no amount of shared emptiness makes it so.

## 6. Unknown sets: do not spend the branch states
A set graph under design B has both branches NULL and its extension elsewhere. That leaves
UNKNOWN unused, and it should stay unused **until a view gives it an exact meaning**. The
candidate meaning is real but not yet earned: a set whose *extent is not yet known here* — a
halo region before exchange, a partition before it is cut. Record it; do not encode it.

The one thing to avoid is spending UNKNOWN on "the extension is stored elsewhere", which is
true of every set under design B and therefore says nothing.

## 7. What this map does not decide
- Which representation family is canonical (counted, listed, bitset, device) — that is a
  measurement, once a caller needs a second one.
- Whether the set graph's branches carry anything at all, or stay (NULL, NULL) forever. A set
  is a leaf in the control plane; the honest default is that it holds nothing.
- Where `name` lives. It is metadata either way, and the choice costs nothing.

## 8. Next
The subset collapse is the highest-value single change and the most testable: `subset_set` has
one construction gate, one ambient chain and five overrides, and every consumer already asks
`is_subobject_of` rather than reading `host`. It should not move before the relation map's
signature decision, because a subset relation is a relation.
