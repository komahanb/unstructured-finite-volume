# Mathematical Architecture

This document is normative. It states the mathematics the code must implement.
When code and mathematics disagree, the code changes.

It covers Levels 0 through 3 — set, relation, relation algebra, related graph.
Levels 4 and above (interpretation, field calculus, discretization, minimization,
constitution, statement) are governed by `AGENTS.md` until they are stated here.

## Level 0 -- Set

A graph is a pair

```text
G = (S, R)
```

where `S` is a collection of declared sets and `R` is a collection of declared
relation objects over those sets.

The root hierarchy is:

```text
graph
├── unrelated_graph        |R| = 0
│   └── set                |S| = 1, |R| = 0
│       └── subset         set with host-inclusion law
└── related_graph          |R| > 0
```

`|R|` counts declared relation objects, not tuples. A `related_graph` may hold
a declared relation whose tuple set is empty.

### Constructions

```text
graph           = (S, R)
unrelated_graph = (S, empty)
set             = ({A}, empty)
subset          = (S <= A), where S is a set and every member of S is in A
related_graph   = (S, R), where at least one relation object is declared
```

`unrelated_graph` is **constructible for arbitrary `S`**:

```text
(empty, empty)          the empty graph
({A}, empty)            one domain
({A,B,C}, empty)        three domains, no structure between them
```

`set` is the `|S| = 1` specialization of that same object, and **no node exists
between them**. An empty set is a different thing again: one declared set with
zero members, which answers `num_sets() = 1`.

Because `set` extends a concrete parent, every set carries the parent's
collection component. It is never allocated there — `set` overrides `num_sets` and
`set_at` to answer from itself — so the cost is one unallocated array
descriptor per set and no heap traffic. That is the price of keeping `set` a
subtype of the object it specializes rather than a sibling of it, and it is
paid deliberately.

### The root contract

Every graph answers exactly six questions, and no more:

```text
id(G)              declared identity
equals(G, H)       identity equality
name(G)            metadata; no law reads it
num_sets(G)        |S|
set_at(G, i)       the i-th declared set
num_relations(G)   |R|
```

Beside these the root carries one **act**, not a question:

```text
declare(G, name)   sign once; a second signing is refused
```

so a graph type exposes seven bindings and answers six questions.

`relation_at` is **not** a root question. Naming a relation requires the Level-1
relation type, which is built over Level 0; a root that could name one would
close the cycle `graph -> set -> relation -> graph`. Counting requires no such
type. `relation_at` is the structure Level 3 adds.

The vertex/edge vocabulary — `edge_tail`, `edge_head`, `incident_edges`,
`adjacent_vertices`, `outgoing_edges`, `incoming_edges` — is not a root
question and never will be. Those are readings of declared relations.

### The `set_at` law

```text
set_at : {1, ..., num_sets(G)} -> S
```

Total on that range and refused outside it, on every graph. When `num_sets(G)`
is zero the range is empty and `set_at` is the empty map; admitting the empty
graph therefore costs no axiom about `S` being occupied.

`set_at` **borrows**, and what it borrows depends on who owns the domain:

```text
set             set_at(A, 1) is A itself          the set IS its one domain
container graph set_at(G, i) is G's own declared set   declared at construction,
                                                  stable for the graph's life
```

A container copies its domains into owned storage at construction, so what it
holds is the declared domain by the token law — a whole-object copy *is* the
identity — and the borrow is into storage the graph controls rather than into
the caller's variable, which may not outlive it. Asked twice, `set_at` answers
the same storage; it never manufactures a fresh copy per call.

The borrower's obligation is the ownership policy used at every level: **the
graph must outlive the borrow**.

### Level-0 laws

Every set has declared identity. A set signs once; a second signing is refused.

Equality is identity equality:

```text
equals(A, B) iff A and B are the same declared mathematical object
```

Equal cardinality does not imply equality. Equal member lists do not imply
equality.

Every set exposes membership, enumeration, and local standing:

```text
size(A)
member_A(i)
members(A)
has_A(a)
local_index_A(a)
```

with the inverse laws:

```text
member_A(local_index_A(a)) = a       for a in A
local_index_A(member_A(i)) = i       for i in 1..size(A)
```

which force enumeration to be injective: a set holds each member once.

A subset is a set with one additional law:

```text
s in S implies s in A
```

where `A` is the declared host set, sealed at construction. Inclusion is
structural; it is not a boolean label on a member list. A subset signs its own
identity — a subobject is its own declared domain, never a disguise of its host.

Embedding is transitive and is asked, never inferred:

```text
is_subobject_of(A, A)                            always
is_subobject_of(S, X)  if  S <= A and is_subobject_of(A, X)
```

### Relation boundary

Level 0 exposes no relation tuples. It counts declared relation objects and
nothing else.

Relations are declared above sets, and the two branches are joined by a
complementary pair of maps:

```text
declare_relations : (S, empty), R  ->  (S, R),  |R| > 0
forget_relations  : (S, R)         ->  (S, empty)
```

```fortran
type(related_graph)   function declare_relations(over, relations, name)
    class(unrelated_graph) , target, intent(in) :: over
    type(declared_relation)        , intent(in) :: relations(:)
    character(len=*)               , intent(in) :: name

type(unrelated_graph) function forget_relations(over, name)
    class(related_graph)   , target, intent(in) :: over
    character(len=*)               , intent(in) :: name
```

**S is preserved by declared set identity**, not merely by extension: each set
of the result `equals` the set it came from, because a whole-object copy is the
declared domain. **Both maps create fresh graph identities** — the result signs
its own name, and `equals(source, result)` is false in every case, including a
round trip, which lands on a third graph and never on the first.

Neither map mutates. The source stands untouched, so borrows already taken from
it remain valid.

They are **module procedures, to preserve level direction**. `declare_relations`
returns a Level-3 type from a Level-0 argument, so binding it to
`unrelated_graph` would make the ground level name the container built over it.
Its partner is written the same way rather than bound to `related_graph`,
because a complementary pair that reads two ways at the call site has stopped
being a pair.

They carry no validation of their own. Both rebuild the collection from the
source's root contract and delegate to the constructor that owns the laws, so
the empty-family, view-ownership and signature-closure refusals are the same
refusals raised in the same place.

Incidence and adjacency are interpretations of declared relations. They are not
methods of the canonical `graph` root.

## Level 1 -- Relation

A relation is a named finite-arity subset of a Cartesian product:

```text
R <= A_1 x ... x A_k,     k >= 1
```

It has identity, arity, an ordered signature, and a membership law:

```text
id(R)             declared identity
equals(R, S)      identity equality
arity(R)          k
domain(R, i)      the set at slot i
has(R, tuple)     membership
num_tuples(R)     |R|, the tuple count
tuples(R)         the tuple table, t(arity, count)
materialized(R)   whole unto itself, or a borrowing view
```

A relation is first-class. It is constructible, queryable and testable without
any graph, and a graph is not required for it to exist.

### Level-1 laws

```text
k >= 1                          a relation relates at least one domain
every slot names a declared set
every tuple has exactly k parts
every part belongs to its slot's declared domain
```

All four are refused at construction, not on access.

A relation signs once, under the same token law as a set. Relations and sets
sign the same roll; the third signer is the related graph.

### Arity readings

For `k = 2` the signature decides the traditional word, and neither is a
separate primitive:

```text
R <= A x A     adjacency
R <= A x B     incidence
```

A binary relation additionally answers `source`, `target`, `image`, `preimage`,
and admits `transpose_of`. Fibre views (`image_view`, `preimage_view`) borrow
their base rather than materializing it.

The relational face of a subset is a Level-1 object, not a Level-0 method:

```text
inclusion_of(S) = I_S <= S x A,   s |-> s
```

total, functional and injective by construction.

### Multiplicity

A relation is a set of tuples; duplicate tuples have no separate identity. When
repeated connectors must be distinguished, the connector is a member: introduce
an edge set and relate through it. Relation-first does not mean edge-free.

## Level 2 -- Relation algebra

Relations generate relations. Every construction here is a **view**: it
references its source rather than duplicating topology, and it is not
materialized.

```text
sigma R           slot permutation
R^T = tau R       transpose, tau = (12); a binary specialization of permutation
project(R, s)     projection onto selected slots
restrict(R, i, P) restriction of slot i to P
S o R             binary composition
```

### Level-2 laws

```text
sigma^-1(sigma R) = R
(R^T)^T = R
R o 1_A = R,    1_B o R = R
T o (S o R) = (T o S) o R
```

Every tuple of a restriction satisfies the restricting predicate. Projection
returns exactly the selected coordinates under the documented multiplicity
semantics.

### Level-2 refusals

```text
a slot index must name a slot of the relation
a restriction domain must embed in the slot it restricts   (is_subobject_of)
a projection selects at least one slot
a projection selects each slot at most once
composition takes two binary relations
composition requires one shared middle domain
```

The restriction law is the reason `is_subobject_of` is a Level-0 question: the
algebra asks it, and a cardinality test would not do.

## Level 3 -- Related graph

```text
related_graph = (S, R),   |R| > 0
```

`related_graph` extends `graph` and adds exactly one question the root cannot
answer:

```text
relation_at(G, i)   the i-th declared relation
```

together with `holds_set(G, A)`, which is `equals` against the declared
collection.

### Level-3 laws

**Signature validity.** Every slot of every declared relation `equals` one of
the graph's own declared sets, or the graph refuses to exist.

**Identity, not signature.** Two relations of one signature coexist freely.
Relation identity is the address, never the signature.

**Collections hold each member once.** `S_i /= S_j` and `R_i /= R_j`, by
identity.

**Only whole relations may be owned.** A borrowing view is refused at the door:
copying a view into owned storage would copy a pointer to a base the graph does
not keep alive. Views ride above graph-owned relations, never inside them.

**The relation family is non-empty.** `|R| > 0` is what `related_graph` means.
A graph over which no relation is declared is an `unrelated_graph`, and
`related_graph` is not that type. The count is of declared relation objects,
never of tuples.

**Immutable after construction**, so that borrowers stay valid for the graph's
life.

### Refusal order

Domains are prior to relations, because relations are declared *over* sets:

```text
1  an empty or unsigned element   a graph holds declared domains only
2  the same domain twice          S_i /= S_j
3  an empty relation family       |R| > 0 is what related_graph means
4  a hollow or unsigned relation   a graph holds declared relations only
5  a borrowing view                a graph owns whole relations
6  the same relation twice        R_i /= R_j
7  an unheld domain               every slot names one of the graph's own sets
```

### Ownership

The graph owns stable relations; interpretations and fibre views borrow them.
Accessors answer pointers into owned storage, never copies, so a borrower may
hold views for as long as the graph lives.

## Legacy compatibility -- `ordinary_graph`

`ordinary_graph` is **not canonical mathematics**. It is the vertex/edge
compatibility contract that predates this hierarchy, retained so that existing
solver, mesh and partition code keeps working during migration.

```text
ordinary_graph      legacy vertex/edge contract        src/graph_grammar.f90
stored_graph        its concrete realization           src/class_graph.f90
mesh_graph          the quarantined mesh-loader stack  src/interface_graph.f90
```

Three rules govern it:

1. **It does not own the word `graph`.** The canonical root is `graph`, in
   `graph_set`. `ordinary_graph` is a separate name for a separate, weaker idea.
2. **Its vocabulary is not promoted.** `edge_tail`, `edge_head`,
   `incident_edges`, `adjacent_vertices`, `outgoing_edges`, `incoming_edges`
   stay where they are. Nothing in Levels 0–3 answers them.
3. **The canonical reading already exists.** `graph_interpretation` derives that
   entire vocabulary from two declared binary relations over a `related_graph`,

   ```text
   T <= E x V     the tail: every edge, exactly one
   H <= E x V     the head: at most one — a boundary half-edge is an ABSENCE
   ```

   storing no tail array, no head array and no adjacency of its own. Migration
   means routing callers to that interpretation, not reimplementing it.

`ordinary_graph` may keep concrete types in its own contract where a canonical
contract could not — its `vertex_set`/`edge_set` answer `index_set` by value,
which avoids a polymorphic allocation at every one of its call sites. That is a
legacy performance affordance, not a statement about the taxonomy.

## Canonical names

```text
graph_set           module owning graph, unrelated_graph, set, subset
declared_set        one declared element of S; a constructor wrapper, not a rung
declared_relation   one declared element of R; a constructor wrapper, not a rung
graph_relation      module owning relation and its stored realization
graph_relation_algebra   module owning the views
graph_structure     module owning related_graph
graph_interpretation     module owning the derived graph-theoretic readings

graph               root of G = (S, R)
unrelated_graph     |R| = 0
set                 declared finite set; |S| = 1, |R| = 0
subset              set with host-inclusion law
related_graph       |R| > 0
relation            finite-arity declared relation
index_set           a REALIZATION of set, A = {1, ..., n} — not a rung
equals              identity equality
is_subobject_of     identity-or-inclusion through host ancestry
```

`index_set` is a concrete realization of `set`, not a level of the taxonomy. A
domain that must list its members arrives as a second realization beside it,
neither above nor below.

Legacy names may exist only as temporary adapters. No canonical module should
introduce or preserve a non-mathematical synonym for `set`.

## Realization gaps

**Closed.** Every graph the laws admit is now constructible. `unrelated_graph`
is concrete and takes an arbitrary collection of declared sets, so
`(empty, empty)` and `(S, empty)` with `|S| >= 2` are declarable objects rather
than admissible fictions. No type was added to reach this: the node that
already meant `G = (S, empty)` was made constructible, and `set` remains its
`|S| = 1` specialization.

**No gap stands.** The relation-family maps are implemented in
`graph_structure` and specified under *Relation boundary* above. Basis pruning
and enrichment are to be described through them — relation removal, addition
and refinement — never through graph mutation.
