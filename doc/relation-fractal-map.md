# Relation, mapped onto the fractal graph
Analysis only. Nothing in `src/` is changed; the prototype is `test/fractal-map/relation.f90`,
run on every pass. Read after `member-set-fractal-map.md`, whose finding this one repeats.

## 0. The same finding, one level up
**A tuple component is already an integer, not a graph.** `stored_relation` holds
`integer :: entry(:,:)`; `csr_relation` holds four integer arrays. No relation in this
repository allocates one object per tuple, and none ever did.

So the three things §9 asks to separate are separable **today**, and the split is again one of
naming and identity rather than of storage:

```
R                 the semantic relation      -> a graph, identity, O(1)
(A_1,...,A_k)     the ordered signature      -> small, a sequence over set graphs
extension(R)      the tuples                 -> large, integers, already flat
```

## 1. Every symbol, classified
| symbol | role | disposition |
|---|---|---|
| `relation` (abstract type) | **representation contract** | does not survive as ontology |
| `identity`, `declare`, `id`, `same_as` | **graph identity** | move to the relation's `type(graph)` |
| `label`, `name` | metadata | external map |
| `arity` | **view** — the signature's length | `sequence_size` of the signature branch |
| `domain(k)` | **view** — the k-th slot | `sequence_element` of the signature branch |
| `slot` | **redundant** | existed because a Fortran array carries one dynamic type; a sequence of set *graphs* is homogeneous |
| `has`, `num_tuples`, `tuples` | **representation** | the extension's questions |
| `stored_relation` | **representation** | `graph_relation_table_representation` |
| `csr_relation` | **representation** | `graph_relation_csr_representation` |
| `binary_relation` (abstract type) | **redundant as a type** | arity-2 views and algorithms |
| `source` / `target` | **view** | `domain_at(R,1)` / `domain_at(R,2)`, already literally that |
| `image` / `preimage` | **algorithm** over a representation | keep; they are the hot path |
| `image_view` / `preimage_view` | **representation** primitive | the borrow the algorithm rides |
| `transposed_view` | **view** | `graph_relation_transpose_view`, a role permutation |
| `transpose_of` | **algorithm/view constructor** | returns the permuted view |
| `inclusion_of` | **algorithm** producing a **representation** | derivable; see §7 |
| `materialized` | **redundant** | a role marker; see §6 |

`source`/`target` already read `domain(1)`/`domain(2)` verbatim in `graph_binary_relation.f90`.
The evidence that arity two needs no ontology is in the file that defines it.

## 2. Signature: a sequence, and one mechanism for every arity
`test/fractal-map/relation.f90` §1 builds three signatures over `graph_sequence_view`:

```
R ⊆ A × A     repeated domain — the SAME set graph in both slots, by identity
R ⊆ A × B     mixed domains
R ⊆ A × B × C arity three — a longer sequence, not a special case
```

All three read through `sequence_size` / `sequence_element`. Cost is O(k) cells with k small —
seven cells for all three relations together. This is design A from the set map, used where it
belongs: the signature *is* control-plane data.

## 3. Where the extension lives — and what a representation is not
`test/fractal-map/relation.f90` §3 binds **one** relation graph R to **two** representations, a
table and a CSR of the same three pairs. Measured:

- both answer the same membership questions, tuple for tuple, including a non-member;
- `pt % same_as(pc)` is **false** — they are different objects;
- R is one relation, whichever storage answers.

That is the split stated as a law: **relation identity belongs to R; representation identity, if
any, is storage identity and is a different question.** The current types answer `same_as`
because they inherited an identity block they should not have had; a representation does not
need to be told apart from another representation of the same relation, and if it ever does,
that is a storage question with a storage answer.

## 4. Tuple equality — the trap, built
Load-bearing, so it is demonstrated rather than argued. `relation.f90` §2 builds two tuple-holder
graphs independently, both spelling `(a, b)`:

```
h1 % same_as(h2)   = .false.     different holders
same_tuple(h1, h2) = .true.      same tuple, by component identity
same_tuple(h1, h3) = .false.     (b,a) is a different tuple
```

If a relation ever deduped on **holder** identity it would hold `(a,b)` twice and stop being a
set. The law:

> `(a_1..a_k) = (b_1..b_k)` iff `a_i % same_as(b_i)` for every i.

The integer table is the representation that makes this automatic — integer equality within a
signature slot *is* component identity, because a member value is meaningful only relative to
the set the slot names. Any future tuple-as-graph representation owes this dedupe explicitly.

## 5. Binary relations need no type
Arity two earns four questions — `source`, `target`, `image`, `preimage` — and nothing else.
Two of them are already `domain(1)`/`domain(2)`. The other two are algorithms over whichever
representation is present, and CSR is the representation that makes them O(degree).

Proposed decomposition, no graph subtype anywhere:

```
graph_relation_view              arity, domain_at, over any relation graph
graph_binary_relation_view       source, target, image, preimage — applicable when arity == 2
graph_relation_csr_representation   the storage those algorithms ride
```

The applicability test is `arity(R) == 2`, answered by the signature. A view that refuses at
arity ≠ 2 is a view with a domain, which is ordinary; a *type* that exists to encode "arity is
2" duplicates what the signature already says.

## 6. materialized() does not survive
It has exactly three inhabitants and the previous pass proved the partition exact:

```
stored_relation   .true.     a representation
csr_relation      .true.     a representation
transposed_view   .false.    a view
```

The predicate does not distinguish two kinds of relation. It distinguishes **a representation
from a view** — two architectural roles sharing one type hierarchy. Once they are separate
modules the question is answered by which module the object came from, and the boolean has
nothing left to say. A boolean whose value is determined by an object's *kind* is evidence that
the kind is not in the type system yet.

Note what it must not become: the binding's storage refusal ("a view cannot be bound") is a real
law and survives — but it will be stated as *only a representation may be bound*, which is a
statement about roles, not a flag consulted at run time.

## 7. Inclusion is derivable, and mostly need not exist
`inclusion_of(S)` currently builds a whole `csr_relation` with `size(S)` tuples of the form
`(s, s)`. Under the set map's finding — a subset reuses the **ambient's member values** — that
relation carries no information the two representations do not already have: `I_S` holds `(s,s)`
exactly when `S % has(s)`.

So it is a **compiled representation of a derivable view**, and it should be built only where a
caller needs a relation *object* to hand to the algebra (composition-with-inclusion is
restriction). Where the question is membership, the view answers without materialising anything.

## 8. Transpose: what its involution preserves, measured
`relation.f90` §4 measures the production view first:

| level of equality | `transpose_of(transpose_of(R))` today |
|---|---|
| graph identity | — (no graph) |
| relation identity | **no** — `transpose_of` mints a fresh token, so T(T(R)) is a *third* identity |
| extension | **yes** — asserted, tuple for tuple |
| representation | no — the view stores nothing |

So the current law is **extensional involution only**, and the report should say so rather than
claim more.

The prototype then shows what a role permutation buys. A view is `(base, role)` where `role` is
a permutation of the slots; transpose swaps it; two transposes compose to the identity
permutation, and the pair denotes the base itself:

```
apply_role(T(v), 1)    = 2         the slots swap
apply_role(T(T(v)), 1) = 1         the permutation is the identity
T(T(v)) % base same_as R           the SAME semantic relation
```

**Double transpose can return the original semantic object** — asserted — provided transpose is
a role permutation over a relation graph rather than a constructor of a new relation. Nothing in
Fortran makes this artificial: the view is a pointer and a two-element permutation.

That pointer is not free. A role view **borrows** its base and inherits exactly the lifetime
obligation `relational_binding` was gated on: the base must outlive every view taken of it, and
the borrow must survive whatever the base does in between. Whoever implements this owes that
gate a test before the view is public.

The classification is **involution**, and it is an involution on *views*, not on relations: R has
no transpose until a view reads its signature in an order.

## 9. Storage order is not identity
A relation is a set of tuples. CSR row order, table column order and declaration order are
**representation properties**. The current code already honours this — `stored_relation`
collapses duplicates keeping first appearance, and `graph_profile` canonicalises output by the
carrier's `local_index` rather than by numeric value. Nothing in this map changes it; it is
recorded so that no future compiled representation is allowed to smuggle its order into the
relation's meaning.

## 10. Scale
`test/fractal-map/scale.f90`, at 10^10 tuples of arity 2:

| | semantic storage |
|---|---|
| A fully extensional (holder + two component cells per tuple) | **1680 GB** of graph objects |
| B semantic graph + CSR (targets + offsets, both directions) | **88 GB** of flat integer arrays |

Design A is rejected. Design B is what exists. No hot CSR operation traverses a pointer, and the
flat arrays are device-transferable as they stand.

## 11. What this map does not decide
- Whether the signature branch is `branch(1)` of the relation graph, and what `branch(2)` holds
  — that is the (S,P)-style choice a *view* makes, and there is more than one defensible answer.
- Whether `image`/`preimage` keep the borrowed-fibre pair or gain a device form.
- Whether relation representations get identities of their own (§3 says only if a caller earns
  one).
