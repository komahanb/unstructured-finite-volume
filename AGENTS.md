# AGENTS.md — Relation-Centered Mathematical Tower

## Mission

Redesign the structural core of `unstructured-finite-volume` from first principles.

Do not preserve an old abstraction merely because the repository already depends on it.
The architecture is allowed to evolve when a more general mathematical object is discovered.

The present redesign supersedes the earlier incidence-first proposal.

The fundamental idea is no longer

\[
G=(A,B,I),
\]

because that still assumes a binary relation between exactly two member sets.

The new foundation is a family of typed finite-arity relations:

\[
R_\rho \subseteq
A_1\times A_2\times\cdots\times A_k,
\qquad k\ge 1.
\]

A graph is then a structured collection of member sets and relations:

\[
\boxed{
G=(\mathcal S,\mathcal R)
}
\]

where

\[
\mathcal S=\{S_1,\ldots,S_n\}
\]

and every relation \(R_\rho\in\mathcal R\) has a declared signature

\[
\operatorname{sig}(R_\rho)
=
(S_{\rho,1},\ldots,S_{\rho,k_\rho}).
\]

The framework must therefore treat **member sets and relations as first-class mathematical objects**.

---

# 1. The governing principle

The code should follow the mathematics that has survived scrutiny, not the history of previous implementations.

Existing principles such as absorption, generation, inheritance, or the current number of tower levels are **heuristics**, not constitutional law.

Keep them when they remain useful.
Discard or modify them when a stronger mathematical structure requires it.

In particular:

> Do not avoid a type hierarchy merely because an older philosophy preferred finite choices as values.

A type or level is justified when it represents a mathematically distinct object with its own laws, valid constructions, or implementations.

---

# 2. The new tower

The intended long-term stratification is:

```text
level 0   carriers             what members may exist
level 1   relations            how members may be related
level 2   relation algebra     how relations generate relations
level 3   relational graph     how relations coexist as one structure
level 4   graph calculus       graph-theoretic interpretations and transforms
level 5   field calculus       values and transport over graph domains
level 6   discretization       how continuous/local laws become discrete operators
level 7   minimization         how residual equations are solved
level 8   constitution         what material/model laws say
level 9   statement            what physical/mathematical problem is asked
```

These level numbers are architectural documentation, not names to encode into public type names.

A future level may be inserted between two existing levels if and only if it corresponds to a genuine mathematical dependency boundary.

Do not preserve consecutive numbering at the expense of a cleaner stratification.

---

# 3. Law of stratification

For every proposed type or module at level \(n\):

1. It may depend only on objects from levels \(0,\ldots,n\).
2. It must introduce one identifiable mathematical commitment not present below.
3. Moving downward must forget structure.
4. Moving upward must add structure, interpretation, or law.
5. A level must not exist solely because one source file became large.
6. A level must not combine two independent mathematical commitments merely to preserve the old tower.

The question for placement is:

> What additional mathematical statement becomes true at this level that was not true below it?

---

# 4. Level 0 — carriers

## 4.1 Fundamental object

The ground object is a **member set**.

Conceptually:

\[
A=\{a_1,\ldots,a_n\}.
\]

A member set has identity independent of its storage representation.

Two member sets are equal only when they are the same declared structural domain, not merely because:

- they contain equal integers;
- they have equal cardinality;
- they use the same Fortran type;
- they happen to enumerate identical values.

Examples:

```text
cells
faces
mesh_edges
points
operations
values
ports
parts
tags
```

are distinct member sets even if every one is represented by integer indices.

## 4.2 Candidate contract

Conceptually:

```fortran
type, abstract :: member_set
contains
    procedure :: id
    procedure :: size
    procedure :: members
end type
```

The exact interface may evolve.

Do not put graph-specific words such as `vertex` or `edge` here.

## 4.3 Identity

Member-set identity is structural.

Do not infer equality from array equality.

A relation signature must refer to member-set identities.

---

# 5. Level 1 — relations

## 5.1 Fundamental object

A relation is a named finite-arity subset of a Cartesian product:

\[
R
\subseteq
A_1\times\cdots\times A_k.
\]

The integer \(k\) is the relation's arity.

A relation therefore has:

```text
identity
arity
ordered signature
tuples / membership law
```

Conceptually:

```fortran
type, abstract :: relation
contains
    procedure :: id
    procedure :: arity
    procedure :: domain          ! slot -> member_set
    procedure :: contains
    procedure :: tuples          ! non-hot generic access
end type
```

Do not require this exact API if a more efficient Fortran contract emerges, but preserve the mathematics.

## 5.2 Relations are first-class

`relation` is not merely an implementation detail of `graph`.

Relations must be independently constructible, queryable, transformable, testable, and composable.

A graph contains relations.

A relation does not need a graph in order to exist mathematically.

This is intentional.

## 5.3 Arity specializations

Specialized abstract descendants are allowed when they add genuine laws or useful compile-time contracts.

Reasonable candidates include:

```text
relation
├── predicate          ! arity 1
├── binary_relation    ! arity 2
└── polyadic_relation  ! arity >= 3
```

Do not create these merely for naming aesthetics.

They are justified only if their interfaces or laws differ meaningfully.

For example, `binary_relation` may legitimately expose:

```text
source
target
image
preimage
transpose
```

because these notions are canonical for arity two.

A generic relation remains the common root.

---

# 6. Unary relations — predicates and supports

A unary relation

\[
P\subseteq A
\]

is a predicate or subset of \(A\).

This is the natural home of the concept currently called `support`.

Therefore:

\[
\boxed{
\text{support is not an edgeless graph}
}
\]

Conceptually, a support is:

```text
predicate over a member set
```

or equivalently an inclusion

\[
S\hookrightarrow A.
\]

Examples:

```text
all cells
boundary faces
owned cells
borrowed cells
wall faces
active degrees of freedom
```

are predicates/subsets.

Do not force a support to answer meaningless graph questions.

This removes the current need for degenerate implementations of:

```text
edge_tail
edge_head
adjacent_vertices
incoming_edges
outgoing_edges
...
```

on a pure membership set.

---

# 7. Binary relations — adjacency and incidence

For

\[
R\subseteq A\times B,
\]

the mathematical primitive is simply a binary relation.

The traditional graph-theoretic word depends on the signature.

If

\[
A=B,
\]

the relation is an **adjacency** relation.

If

\[
A\ne B,
\]

the relation is an **incidence** relation.

Thus:

\[
\boxed{
\text{adjacency and incidence are interpretations of binary relation}
}
\]

not separate primitives.

Examples:

\[
R_{VV}\subseteq V\times V
\]

is vertex adjacency.

\[
R_{EE}\subseteq E\times E
\]

is edge adjacency.

\[
R_{VE}\subseteq V\times E
\]

is vertex-edge incidence.

\[
R_{CF}\subseteq C\times F
\]

is cell-face incidence.

Do not create unrelated storage mechanisms merely because one receives the word `adjacency` and another receives the word `incidence`.

---

# 8. Higher-arity relations

The framework must not assume that every structural fact can or should be reduced to pairwise incidence.

A ternary relation is valid:

\[
R\subseteq A\times B\times C.
\]

Likewise for arbitrary finite arity.

This is the main correction to the earlier incidence-first design.

## 8.1 Ordinary directed edge example

Instead of storing an endpoint role as opaque metadata:

```text
(edge, vertex, role)
```

may itself be a ternary structural relation:

\[
R_{\mathrm{endpoint}}
\subseteq
E\times V\times P,
\]

where

```text
P = {tail, head}
```

is a member set of endpoint roles.

Then an interior edge may contain:

```text
(e, v_tail, tail)
(e, v_head, head)
```

while a boundary face may contain only:

```text
(e, v, tail)
```

No imaginary head is required.

## 8.2 Calculator example

A computation may use:

```text
O = operations
V = values
P = ports / roles
```

with

\[
R_{\mathrm{flow}}
\subseteq O\times V\times P.
\]

Examples:

```text
(+, 2, input)
(+, 3, input)
(+, 5, output)
(*, 5, input)
(*, 4, input)
(*, 20, output)
```

This preserves operation/value/role as structural dimensions rather than encoding role as an incidental integer attribute.

## 8.3 Reification remains optional

Any \(k\)-ary relation can often be represented through an auxiliary relation-object set and binary incidences.

Do not assume this factorization is always the canonical representation.

Choose reification only when the relation instance itself deserves identity or carries independent state.

---

# 9. Level 2 — relation algebra

Relations generate new relations.

This level contains the universal structural algebra, independent of graph semantics.

## 9.1 Slot permutation

For a \(k\)-ary relation

\[
R\subseteq A_1\times\cdots\times A_k
\]

and permutation \(\sigma\),

\[
\sigma R
\subseteq
A_{\sigma(1)}
\times\cdots\times
A_{\sigma(k)}.
\]

For a binary relation, transpose is the special permutation

\[
R^T=\tau R,
\qquad \tau=(12).
\]

Therefore:

\[
\boxed{
\text{transpose is a binary specialization of slot permutation}
}
\]

and should not be the deepest primitive.

## 9.2 Projection

A relation may be projected onto selected slots.

If

\[
R\subseteq A\times B\times C,
\]

then

\[
\pi_{AB}(R)\subseteq A\times B.
\]

Projection may lose information.

That loss must be explicit.

## 9.3 Selection / restriction

A relation may be restricted by predicates on one or more slots.

Conceptually:

\[
R|_{P}
\]

or relational selection.

## 9.4 Join

Relations sharing compatible domains may be joined.

Natural join is the general operation from which many binary compositions can be expressed.

## 9.5 Binary composition

For binary relations

\[
R\subseteq A\times B,
\qquad
S\subseteq B\times C,
\]

define:

\[
S\circ R
=
\{(a,c):\exists b,\;aRb\land bSc\}.
\]

Then incidence generates adjacency:

\[
R_{AB}\circ R_{AB}^T
\subseteq A\times A,
\]

and

\[
R_{AB}^T\circ R_{AB}
\subseteq B\times B.
\]

Do not require adjacency to be derived this way.

A same-set adjacency relation may also be primitive.

## 9.6 Identity relation

Each member set has an identity relation:

\[
1_A=\{(a,a):a\in A\}.
\]

## 9.7 Set-like operations

For compatible signatures, relations may support:

```text
union
intersection
difference
```

when useful.

Do not add operations to the fundamental runtime API merely because they exist mathematically; distinguish the mathematical algebra from the hot-path programming interface.

---

# 10. Relation views and relation implementations

The hierarchy should distinguish **what relation means** from **how it is produced**.

Possible concrete/view families include:

```text
relation
├── stored_relation
├── permuted_relation_view
├── projected_relation_view
├── restricted_relation_view
├── joined_relation_view
├── composed_binary_relation_view
└── generated_relation
```

This hierarchy is encouraged because these are genuinely different realization mechanisms.

A relation view should generally reference its source relation rather than duplicate topology.

Derived relations should be lazy unless materialization is explicitly required.

---

# 11. Relation laws

The test suite must encode mathematical laws.

Examples:

## 11.1 Permutation inverse

\[
\sigma^{-1}(\sigma R)=R.
\]

## 11.2 Binary transpose involution

\[
(R^T)^T=R.
\]

## 11.3 Identity

\[
R\circ 1_A=R,
\qquad
1_B\circ R=R.
\]

## 11.4 Associativity of binary composition

\[
T\circ(S\circ R)
=
(T\circ S)\circ R.
\]

## 11.5 Projection semantics

Projected tuples must equal the mathematical projection of the original tuples, including duplicate policy if the implementation stores multiplicity.

Any departure from set semantics must be stated explicitly.

---

# 12. Relation multiplicity

A pure mathematical relation is a set of tuples.

Duplicate identical tuples have no separate identity.

This is insufficient when repeated connectors themselves matter.

For example, a multigraph may contain two distinct edges between the same vertices.

Do not collapse first-class edge objects into duplicate adjacency tuples.

Instead retain an edge member set:

```text
V = vertices
E = edges
```

and represent endpoint structure through relations involving `E`.

Relation-first does not mean edge-free.

It means edges are members when edge identity matters.

---

# 13. Functional relations

Some relations satisfy stronger laws.

A binary relation

\[
R\subseteq A\times B
\]

may be:

```text
functional
total
injective
surjective
bijective
```

These properties may justify refined interfaces or validators.

Examples in this repository include concepts analogous to:

```text
global index map
part index map
owner part
```

Do not immediately create a subtype for every property.

But do not forbid such subtypes a priori either.

Use a refined type when downstream algorithms genuinely require the stronger law at construction or compile time.

---

# 14. Level 3 — relational graph

## 14.1 Graph definition

A graph is a **relational structure**:

\[
\boxed{
G=(\mathcal S,\mathcal R)
}
\]

where \(\mathcal S\) is a collection of member sets and \(\mathcal R\) is a collection of typed relations over them.

A graph is therefore composition, not inheritance from relation.

Conceptually:

```fortran
type, abstract :: graph
contains
    procedure :: num_member_sets
    procedure :: member_set
    procedure :: num_relations
    procedure :: relation
end type
```

The exact API may be refined.

The key architectural statement is:

\[
\boxed{
\text{Graph contains Relations}
}
\]

not:

\[
\boxed{
\text{Relation is hidden inside Graph}
}
\]

and not:

\[
\boxed{
\text{Graph is itself one Relation}.
}
\]

## 14.2 Relation names and identities

A graph may contain multiple relations with the same signature.

For example:

\[
R_1,R_2\subseteq A\times A
\]

may represent:

```text
physical adjacency
dependency adjacency
coarsening adjacency
communication adjacency
```

Relation identity is therefore independent of signature.

Do not address a relation solely by `(source_set, target_set)`.

---

# 15. Graph signatures and schemas

A graph has a relational signature/schema:

```text
member-set identities
relation identities
relation arities
relation slot domains
required relation laws
```

This schema is structural.

A graph instance supplies the members and tuples satisfying that schema.

This separation may become useful for:

```text
ordinary graph
directed graph
hypergraph
finite-volume topology
finite-element topology
computation graph
time-integration graph
```

A schema may state required relations without prescribing storage.

---

# 16. Graph views / profiles

Traditional graph objects should be expressed as **views or profiles over relational graphs**.

Examples:

```text
ordinary_graph_view
directed_graph_view
hypergraph_view
mesh_topology_view
computation_graph_view
```

A view earns a type when it guarantees additional laws and therefore permits additional meaningful queries.

For example, an ordinary directed graph view may guarantee the existence of:

```text
V = vertex set
E = edge set
endpoint relation
```

and may then expose familiar convenience queries:

```text
num_vertices
num_edges
tail
head
incident_edges
adjacent_vertices
incoming_edges
outgoing_edges
```

These are not Level-1 primitives.

They are graph-view vocabulary derived from the relational structure.

---

# 17. Duality

Duality is no longer foundational.

For a particular two-sort structure with

\[
R\subseteq A\times B,
\]

the familiar dual exchanges the two domains using transpose:

\[
(A,B,R)
\mapsto
(B,A,R^T).
\]

For higher-arity relations the analogous structural operation is a permutation of slots.

Therefore the generic concept is:

\[
\boxed{
\text{reorientation / permutation of relation slots}
}
\]

while `dual()` remains a meaningful higher-level graph view when a specific graph schema defines it.

Do not force every graph to have a dual.

---

# 18. Level 4 — graph calculus

This level introduces graph-theoretic interpretations and structural algorithms over relational graphs.

Examples include:

```text
walk
reachability
connected components
colouring
partition
assemble
coarsen
refine
graph quotient
graph condensation
graph refinement
```

These algorithms consume graphs or graph views.

They do not become methods on the fundamental relation or graph merely because they use them.

## 18.1 Adjacency queries

Adjacency may be:

```text
primitive relation
derived composition
cached/materialized view
```

depending on graph schema and performance needs.

The semantic source of truth must be explicit.

## 18.2 Boundary

Boundary is not a universal relation concept.

For a finite-volume topology, a boundary face can be defined through cell-face incidence degree.

For another graph schema, `boundary` may mean something else or nothing at all.

Keep boundary semantics in the relevant graph view/schema.

## 18.3 Direction

Incoming/outgoing are interpretations of orientation in a directed graph schema.

They are not universal methods of `relation`.

---

# 19. Graph transforms

Partition, assemble, coarsen, and refine remain graph-to-graph constructions.

Their contracts should consume relational graphs and preserve or transform relations according to explicit laws.

Do not hard-code vertex/edge as the only transferable member sets.

A transform should describe:

```text
which member sets are transformed
how each relation is mapped
which relation laws are preserved
how fields are transported
```

Graph transforms may themselves induce relations between source and target graphs.

---

# 20. Level 5 — field calculus

## 20.1 Field definition

A field is a function over a domain.

If \(A\) is a member set and \(V\) a value space:

\[
f:A\to V.
\]

If the field exists only on a support \(P\subseteq A\), then:

\[
f:P\to V.
\]

A field domain is therefore a member set or predicate/support, not an edgeless graph.

## 20.2 Field contract

Conceptually retain:

```text
name
units
domain
value kind
number of components
get/set vector
```

but make `domain` refer to the new domain abstraction.

Do not duplicate structural type information such as:

```text
vertex_field
edge_field
cell_field
face_field
```

when the domain already states where values live.

Introduce specialized field types only if their value algebra differs.

---

# 21. Relation-induced field transport

Relations naturally induce movement of data between domains.

For binary

\[
R\subseteq A\times B,
\]

field operations may include abstract forms of:

```text
gather
scatter
push
pull
reduce-over-fibre
broadcast-over-fibre
```

depending on coefficient/aggregation laws.

This is the proper mathematical neighborhood for many current graph operations.

Do not confuse relation structure with the numerical rule used to combine values.

The relation says **which values are related**.

The field operation says **what arithmetic is performed over those related values**.

---

# 22. Reduction and broadcast

Reduction and broadcast remain field-calculus operations.

A reduction maps many field entries to a functional.

A broadcast maps a functional to many field entries.

Their domains should be expressed using the new member-set/support model.

Do not couple reduction semantics to vertices.

---

# 23. Level 6 — discretization

This level binds structural relations and field arithmetic into discrete numerical operators.

Examples:

```text
differential operator
stencil operator
discretization operator
linearization operator
balance
interpolation
gradient
divergence
laplacian
```

The distinction from Level 5 is:

> Level 5 knows values move and combine over relations.
> Level 6 knows a particular discrete mathematical law.

## 23.1 Structural kernel

Current difference/incidence operations should be recognized as algebra over oriented relations.

For an oriented binary representation \(C\):

\[
C
\quad\text{and}\quad
C^T
\]

are structural actions.

The current formulas such as

\[
Dz=M^{-1}Cz
\]

and

\[
Gq=-H^{-1}C^Tq
\]

should reuse common relation traversal where possible.

## 23.2 Higher-arity relations

Do not assume every discretization must first collapse its structure to an ordinary edge graph.

A discretization may operate directly on higher-arity relations if that is the natural stencil structure.

---

# 24. Dependency structure

A discretization's dependency pattern is itself a relation or relational graph.

Do not require dependency structure to be an ordinary graph unless the algorithm requires one.

For a field on \(A\), a dependency relation may naturally be:

\[
D\subseteq A\times A.
\]

For multi-field operators, dependencies may span several member sets.

This should become more expressive than the current single ordinary-graph dependency pattern.

---

# 25. Linearization and adjoints

Reverse-mode traversal is fundamentally related to reversal/permutation of structural relations plus transpose of numerical actions.

For binary linear actions:

\[
\langle Cx,y\rangle
=
\langle x,C^Ty\rangle.
\]

For general relational structures, the reverse pass should follow the opposite data-dependency relation/view.

Do not duplicate topological traversal logic for forward and reverse modes when relation algebra can express the reversal.

---

# 26. Level 7 — minimization

Minimization remains the level where a residual equation is driven toward a solution.

It may consume:

```text
discrete operation
dependency relations
field inner products
norms
colourings
preconditioner structures
```

but must not assume every unknown lives on "vertices".

Unknown domains are member sets/supports.

Replace vertex-specific assumptions in minimizers with domain-driven logic as the migration reaches this level.

---

# 27. Level 8 — constitution

Constitution binds material/model behavior.

It may introduce physical vocabulary and constitutive laws.

It should consume discretization/field abstractions rather than graph-storage details.

No relation primitive should learn words such as:

```text
conductivity
viscosity
density
pressure
temperature
```

unless those words are themselves user-defined data carried by fields.

---

# 28. Level 9 — statement

The statement level expresses the problem being asked.

It combines:

```text
domains
fields
constitutions
sources
boundary conditions
discretizations
objectives/residuals
```

without reaching downward into raw relation storage.

---

# 29. Reinterpretation of current `graph_grammar`

The existing `graph_grammar` currently places four roots at Level 0:

```text
graph
graph_field
graph_operation
graph_transform
```

and defines graph structurally through vertices, edges, tail, and head.

This structure is no longer the desired ground level.

Do not mechanically patch more relation methods into the existing `graph` root.

Instead migrate toward the new tower:

```text
member_set
relation
relation algebra
graph
field
operation/transform
...
```

The current four-root grammar may remain temporarily as a compatibility layer during migration.

It is not the target ontology.

---

# 30. Reinterpretation of current `stored_graph`

The current `stored_graph` uses:

```text
nv
ne
tail(:)
head(:)

incident-edge CSR
adjacent-vertex CSR
outgoing-edge CSR
incoming-edge CSR
```

as its source representation.

The target should separate:

```text
member sets
stored relations
derived/cached relation views
graph profile
```

For the current ordinary directed graph, a viable transitional schema is:

```text
V = vertices
E = edges
P = endpoint roles

endpoint relation:
    R_endpoint ⊆ E × V × P
```

or an equivalent pair of binary endpoint relations if benchmarking or simplicity strongly favors them.

Do not decide this representation solely by aesthetic purity.
Preserve the general graph model while allowing optimized schema-specific storage.

---

# 31. Hot-path specialization is allowed

The abstract mathematics must be general.

The implementation may be specialized.

For example:

```text
generic relation interface
    ↓
stored_binary_relation_CSR
    ↓
ordinary_graph_view
```

is acceptable.

Likewise a ternary relation may use storage optimized for one slot ordering.

The requirement is semantic equivalence, not one universal storage layout.

---

# 32. Query complexity

Hot graph traversal must remain asymptotically efficient.

For a stored binary relation:

```text
image/member lookup      O(degree)
preimage lookup          O(degree)
transpose view creation  O(1)
size/count               O(1)
```

where indexes are available.

For higher-arity relations, explicitly document which slot projections are indexed.

Do not promise O(degree) lookup on every slot unless the required index exists.

A relation may expose capabilities or indexed-slot metadata if needed.

---

# 33. No hidden allocation in hot loops

Do not introduce polymorphic allocation or tuple materialization inside million-cell loops merely to satisfy the generic interface.

Provide low-level views/slices/iterators for hot traversal.

Generic relation algebra may allocate outside hot loops.

Compiled/materialized derived relations are permitted when repeated interpretation would be expensive.

---

# 34. Interpretation versus compilation

Preserve the useful distinction already emerging in the codebase:

```text
interpreted relation/operator
compiled/materialized relation/operator
```

An interpreted view reads source relations dynamically.

A compiled view materializes indexes, adjacency, stencil weights, or other derived structure for repeated use.

Both represent the same mathematics.

Do not conflate semantic type with execution strategy.

---

# 35. Tags

Tags should be reconsidered relationally.

A tag system may be represented as:

```text
T = tag names
R_tag ⊆ A × T
```

rather than one fixed character array attached separately to every member-set kind.

This need not be migrated immediately.

The important point is that tagging is a relation between members and labels, not a special property of vertices and edges.

---

# 36. Partition frame

Partition metadata should also be reconsidered relationally.

Possible member set:

```text
P = parts
```

Possible relations:

\[
R_{\mathrm{owner}}\subseteq A\times P
\]

for ownership, with a functional-law constraint.

Borrowed/overlap membership may be predicates or relations between local structures and global structures.

Global/local index correspondence may be a functional or bijective relation.

Do not force all of these into the first refactor.
But design the new graph core so they can become relations naturally.

---

# 37. Support migration

The existing `support` class is a prime early migration target.

Target:

```text
support = predicate/subset over a member set
```

not:

```text
support extends graph with zero edges
```

Migration rule:

1. introduce the new predicate/support abstraction;
2. adapt field domains to accept it;
3. migrate named-set queries;
4. remove degenerate graph procedures from support;
5. only then retire the old graph-support inheritance.

---

# 38. Field migration

Fields currently reference graph-shaped domains.

Target:

```text
field.domain -> member_set or support
```

A field on all cells references the cell member set.

A field on boundary faces references a boundary-face predicate/support.

The field should not need to know whether its domain is called vertex, edge, cell, face, operation, or value.

---

# 39. Partition / assembly migration

Current vertex-field and edge-field branches should eventually become domain-parametric.

Prefer:

```text
carry_field(domain, ...)
assemble_field(domain, ...)
```

over:

```text
carry_vertex_field
carry_edge_field
assemble_vertex_field
assemble_edge_field
```

But the new design should go beyond a two-valued `side`.

The domain is a first-class member set/support.

Do not replace `vertex/edge` with merely `A/B` and declare the problem solved.

---

# 40. Differential-operator migration

The current differential operator explicitly alternates:

```text
vertices -> edges
edges    -> vertices
```

and already reasons in transpose pairs.

Migration should:

1. identify the structural relation used by each elementary step;
2. express forward and reverse traversal through relation views/permutations;
3. separate relation topology from numerical coefficients;
4. preserve current signs and boundary semantics;
5. test adjoint identities;
6. permit future operators over relations beyond vertex-edge incidence.

---

# 41. Graph calculus migration

The current walk, colouring, partitioning, coarsening, and refinement algorithms may initially consume an ordinary-graph view.

Do not rewrite every algorithm against arbitrary \(k\)-ary relations immediately.

Instead:

```text
general relational graph
    ↓ profile/view
ordinary graph view
    ↓
existing algorithms
```

This gives generality without destroying mature graph algorithms.

---

# 42. Ordinary graph compatibility

The ordinary graph API remains useful.

A compatibility profile should continue to provide:

```text
num_vertices
num_edges
edge_tail
edge_head
edge_has_head
incident_edges
adjacent_vertices
outgoing_edges
incoming_edges
outgoing_vertices
incoming_vertices
```

while these are derived from or backed by the new relational graph.

Do not remove user-facing graph vocabulary merely because the internal ontology changed.

---

# 43. Hypergraph compatibility

Hypergraphs fit naturally.

Use:

```text
V = vertices
E = hyperedges
R_incident ⊆ V × E
```

or a higher-arity relation when roles/order require it.

Arbitrary relation degree must be permitted.

No root contract may assume an edge has exactly two incident vertices.

---

# 44. Mesh topology

A geometric mesh naturally has several member sets:

\[
X_0,X_1,X_2,X_3
\]

typically:

```text
X0 = points
X1 = edges
X2 = faces
X3 = cells
```

and multiple relations:

\[
R_{01},R_{12},R_{23},R_{02},R_{13},\ldots
\]

as primitive or derived structure.

Same-dimensional adjacency relations may also exist:

\[
R_{00},R_{11},R_{22},R_{33}.
\]

Do not force the entire mesh into one ordinary graph.

A finite-volume solver may choose the cell-face relation as its primary operational view without making it the universal mesh ontology.

---

# 45. Composition of graphs

Since a graph is a collection of member sets and relations, graph composition should eventually be expressible as composition/union of compatible relational structures.

Potential operations include:

```text
merge schemas
join structures on shared member sets
restrict to substructure
quotient
refine
coarsen
compose dependency structures
```

Do not define these all now.

But do not design `graph` so narrowly that composition becomes impossible without rebuilding the object.

---

# 46. Subgraphs and substructures

A subgraph should be reconsidered as a relational substructure:

```text
selected member subsets
relations restricted to those subsets
```

This is more general than selecting vertices and induced edges only.

Induced ordinary subgraphs become one specialization.

---

# 47. Quotients and refinement

A quotient graph should be formulated through maps/relations from fine member sets to coarse member sets.

Refinement is the complementary direction only when explicit laws make it so.

Do not hard-code quotient/refinement semantics into the relation root.

They belong to graph calculus and transforms.

---

# 48. Dependency graphs

Dependencies are relations.

A computation dependency structure should not automatically introduce fake edge objects.

If edge identity is irrelevant:

\[
D\subseteq A\times A
\]

is sufficient.

If dependency instances require identity, introduce a separate member set.

Choose representation according to what is first-class in the problem.

---

# 49. Structural roles as domains

Before adding an integer or enum attribute to a relation tuple, ask:

> Is this attribute merely implementation metadata, or is it itself a structural axis whose values participate in laws?

If structural, prefer giving it a member set and relation slot.

Examples that may deserve domains:

```text
tail/head
input/output
port number
time level
derivative order
part
material region
```

Do not over-reify trivial flags.
Use mathematical judgment.

---

# 50. Type hierarchy admission law

The new hierarchy rule is:

A type earns existence if at least one of the following is true:

1. it has additional mathematical laws;
2. it admits operations not meaningful on its parent;
3. it restricts valid constructions in a useful way;
4. it has multiple concrete realization strategies;
5. downstream code benefits materially from demanding that stronger contract.

A type does **not** earn existence merely because:

- a noun exists in English;
- one enum value has a convenient name;
- a source file would otherwise be long;
- old code already used that name.

This replaces any blanket rule that finite distinctions must always be absorbed as values.

---

# 51. Inheritance versus composition

Use inheritance for:

```text
is-a mathematical refinement
```

Use composition for:

```text
has-a mathematical constituent
```

Therefore:

```text
binary_relation is-a relation
stored_relation is-a relation realization
graph has member_sets
graph has relations
field has a domain
```

A graph is not a relation.

A support is not a graph.

A field is not a graph.

Do not use inheritance merely to reuse procedures.

---

# 52. Laws before convenience

Every abstraction should state its laws before its convenience API.

Examples:

```text
relation transpose involution
slot permutation inverse
composition associativity
support subset validity
ownership functionality
partition/assembly round trip
adjoint transpose identity
```

A law belongs in tests if Fortran's type system cannot enforce it.

---

# 53. Capability over fiction

Do not make a generic relation pretend it can answer an efficient query it cannot support.

If only some slot orderings are indexed, expose that truth.

Possible designs include:

```text
supports_lookup(slot)
supports_projection(slots)
supports_transpose()
```

or concrete view requirements.

Avoid hidden O(N) scans behind methods that look like O(degree) neighborhood queries.

---

# 54. Immutability

Structural objects should remain immutable after construction by default:

```text
member sets
relations
graphs
schemas
```

Derived cached indexes may be constructed eagerly or lazily only if externally observable semantics remain immutable and thread-safe.

If mutable graphs are needed later, introduce a distinct mutable-builder or mutable-graph abstraction rather than weakening the immutable contract silently.

---

# 55. Builders

Construction complexity should not contaminate immutable runtime objects.

Consider builder types for:

```text
member sets
stored relations
relational graphs
mesh topology
```

if validation and indexing become substantial.

A builder may mutate.

The built structural object should not.

---

# 56. Validation

Validate structural laws at construction or at explicit validation boundaries.

Examples:

```text
tuple indices belong to declared slot domains
functional relations have at most one image where promised
ordinary edge profile satisfies endpoint laws
partition ownership is total where required
relation signatures reference graph member sets
```

Do not repeatedly revalidate these invariants in hot loops.

---

# 57. Naming

Prefer mathematical names at fundamental levels:

```text
member_set
relation
predicate
signature
slot
permutation
projection
restriction
join
composition
graph
field
```

Use application vocabulary in profiles and higher layers:

```text
vertex
edge
cell
face
tail
head
wall
flux
gradient
```

Do not replace precise established mathematical terminology with vague generic words merely to sound abstract.

---

# 58. Repository module direction

The long-term module layout may evolve toward something conceptually like:

```text
graph_carrier.f90
graph_relation.f90
graph_relation_algebra.f90
graph_structure.f90
graph_calculus.f90
graph_field_calculus.f90
graph_discretization.f90
graph_minimization.f90
...
```

These names are illustrative, not mandatory.

Do not rename modules until their responsibilities are actually separated.

Architecture first, filenames second.

---

# 59. Migration strategy

Proceed by semantic seams, not by mass rewrite.

## Phase 0 — characterization

Before changing behavior:

- preserve current tests;
- add tests for parallel edges, boundary half-edges, directed traversal, supports, partition round trips, and differential adjoints;
- benchmark current hot graph traversals.

## Phase 1 — introduce `member_set`

Create a first-class domain/member-set abstraction.

Use it initially beside current vertex/edge APIs.

## Phase 2 — introduce generic `relation`

Implement a first-class finite-arity relation contract.

Provide at least one stored implementation and slot-signature validation.

## Phase 3 — implement binary specialization

Provide efficient binary relation traversal and transpose view.

Use CSR or equivalent for the current ordinary graph path.

## Phase 4 — migrate ordinary graph structure

Represent current vertex/edge topology through member sets + relation(s).

Keep the old `graph` API as a compatibility view.

## Phase 5 — migrate support

Replace edgeless-graph support with unary predicate/subset semantics.

Adapt fields.

## Phase 6 — relation algebra

Introduce permutation/transpose, restriction, projection, and binary composition as views/constructions as needed by real callers.

Do not implement a complete database engine.

## Phase 7 — migrate partition/assembly

Generalize data transport by actual domains, not vertex/edge branches.

## Phase 8 — migrate differential topology

Unify forward/reverse structural traversal through relation algebra.

## Phase 9 — split graph calculus and field calculus

Move graph-theoretic algorithms and field/value operations into distinct mathematical levels if current module crowding persists.

## Phase 10 — migrate minimization assumptions

Remove vertex-specific unknown-domain assumptions.

## Phase 11 — consolidate legacy graph stacks

Only after the relational architecture is exercised by real code should duplicate old/new graph infrastructures be retired.

---

# 60. What not to do

Do not:

- simply rename `vertex/edge` to `A/B`;
- simply rename `incidence` to `relation` while retaining binary-only assumptions;
- make `relation` a private helper inside `graph`;
- force all relations to be binary;
- force all higher-arity relations through reified edge objects;
- force all adjacency to be derived from incidence;
- force all incidence to be derived from adjacency;
- make support an edgeless graph;
- make graph inherit from relation;
- make relation own physics;
- make graph own fields;
- make fields own topology;
- expose O(N) scans as cheap neighborhood queries;
- remove useful ordinary-graph views merely for abstraction purity;
- refactor every solver layer in the same commit;
- preserve an old law solely because it was once documented as fundamental.

---

# 61. Acceptance laws — relation layer

The relation layer must test at least:

### Domain validity

Every tuple component belongs to its declared slot member set.

### Arity

Each tuple has exactly the relation's declared arity.

### Permutation inverse

\[
\sigma^{-1}(\sigma R)=R.
\]

### Binary transpose

\[
(R^T)^T=R.
\]

### Binary composition identity

\[
R\circ 1=R,
\qquad
1\circ R=R.
\]

### Binary composition associativity

\[
T\circ(S\circ R)
=
(T\circ S)\circ R.
\]

### Restriction

Every tuple in a restricted relation satisfies the restricting predicate.

### Projection

Projection returns exactly the selected tuple coordinates under the documented multiplicity semantics.

---

# 62. Acceptance laws — graph layer

### Signature validity

Every graph relation references member sets belonging to that graph or explicitly imported shared domains.

### Relation identity

Two relations with the same signature may coexist without collision.

### Ordinary graph compatibility

Current ordinary graph queries must reproduce existing behavior through the compatibility profile.

### Parallel-edge identity

Distinct edge members connecting the same vertex members remain distinct.

### Boundary half-edge

A boundary face requires no imaginary remote member.

### Derived adjacency

Where adjacency is defined by composition, it agrees with the ordinary graph profile's adjacency semantics.

---

# 63. Acceptance laws — support/field layer

### Support validity

Every support member belongs to its host member set.

### Field domain validity

Field storage shape agrees with domain cardinality and component count.

### Domain independence

The same field implementation works on cells, faces, vertices, edges, operations, values, or future domains without separate structural subclasses.

---

# 64. Acceptance laws — partition layer

### Ownership validity

Where ownership is declared functional, each owned member has exactly one owner.

### Round trip

\[
\operatorname{assemble}(\operatorname{partition}(G))=G
\]

under existing exact-partition semantics.

### Field round trip

Fields on every supported member set/domain survive partition and assembly.

---

# 65. Acceptance laws — discretization layer

### Structural transpose

For compatible linear structural actions:

\[
\langle Cx,y\rangle
=
\langle x,C^Ty\rangle.
\]

### Adjoint consistency

Existing tangent/adjoint operator tests must continue to pass.

### FV conservation

Interior face contributions cancel with opposite signs across incident cells according to the balance convention.

Boundary faces contribute exactly according to their actual relation tuples.

---

# 66. Performance acceptance

Before deleting old storage, benchmark:

```text
incident query
adjacent query
incoming/outgoing query
field assembly
differential operator apply
partition construction
```

The new abstraction must not impose material regression on hot paths without a documented reason and an optimization plan.

Semantic generality should be obtained primarily through interfaces and views, not through repeated dynamic reconstruction.

---

# 67. Design questions agents must ask

Before implementing a structural feature, ask in this order:

1. **What are the member sets?**
2. **What is the arity of the relation?**
3. **What is its ordered signature?**
4. **Does the relation itself need identity?**
5. **Is it primitive or derived?**
6. **If derived, by permutation, projection, restriction, join, or composition?**
7. **Does a stronger subtype have additional laws?**
8. **Does the graph merely contain this relation, or does a graph profile interpret it?**
9. **What queries are hot and therefore need indexing/materialization?**
10. **What mathematical law will test the implementation?**

Do not start from:

> Is this a vertex thing or an edge thing?

That question belongs much later.

---

# 68. The new conceptual endpoint

The previous endpoint was:

\[
\text{incidence is fundamental; primal and dual are views}.
\]

That is superseded.

The new endpoint is:

\[
\boxed{
\text{members are organized into sets}
}
\]

\[
\boxed{
\text{relations connect tuples of members}
}
\]

\[
\boxed{
\text{relation algebra transforms and composes relations}
}
\]

\[
\boxed{
\text{a graph is a structured collection of member sets and relations}
}
\]

\[
\boxed{
\text{adjacency, incidence, direction, duality, hyperedges, supports, and dependencies are interpretations or specializations of that structure}
}
\]

The implementation should evolve toward this architecture without sacrificing the repository's existing numerical correctness, immutability, or hot-loop performance.
