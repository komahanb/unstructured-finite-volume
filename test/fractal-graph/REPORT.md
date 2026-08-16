# Two-type fractal graph kernel

Spike only. `src/` is unchanged, the tower is unrelevelled, no module is deleted.
Built and run by `test/fractal-graph/run.sh` (gfortran-15, `-std=f2023`).

**Kernel size**, `fractal_graph.f90`, by `grep -c -v -E '^[[:space:]]*(!|$)'` and its
complement — **54 code lines**, 35 comment, 43 blank.

## 1. Declarations of graph and graph_branch

```fortran
integer, parameter :: GRAPH_NULL = 0, GRAPH_UNKNOWN = 1, GRAPH_KNOWN = 2

type :: graph_branch
   integer              :: status = GRAPH_NULL
   type(graph), pointer :: known  => null()
end type graph_branch

type :: graph
   type(graph_branch)   :: branch(2)
   type(token), private :: identity
 contains
   procedure :: declare, id, same_as
end type graph
```

Two types; no third type stands between them. `branch` is public, so
`g % branch(1) % known % branch(2) % known` traverses directly at any depth (test C).

## 2. Does direct recursive referencing work?

Yes, subject to a declaration-order constraint. The three representations are in
`fortran-recursion/` and `run.sh` compiles all three on every run with the compiler and
standard the kernel is built with:

| representation | gfortran-15 `-std=f2023` |
|---|---|
| `branch_before_graph` — `type(graph), pointer :: known` | **compiles**; used by the kernel |
| `polymorphic_known` — `class(graph), pointer :: known` | compiles; unused, nothing extends `graph` |
| `graph_before_branch` | **rejected** — *Derived type at (1) has not been previously defined and so cannot appear in a derived type definition* |

The constraint is on declaration order, not on the ontology. One further constraint:
`known_branch` cannot be `pure`, since pointer assignment to an `INTENT(IN)` target is
prohibited in a pure procedure (F2018 C1594). `null_branch` and `unknown_branch` are pure.

The obstruction found by the previous spike is an obstruction to **ownership**
(`class(graph), allocatable`), not to the two-type ontology. Reference removes it.

## 3. Lifetime rules

1. A branch references; it does not own its target. A graph must outlive every branch that
   references it. The kernel allocates and deallocates nothing.
2. `TARGET` must be on the actual argument, not only the dummy: pointers set into a `TARGET`
   dummy remain associated after return only if the actual argument is also `TARGET` and is
   not an array section with a vector subscript (F2018 15.5.2.4). Every graph in this spike is
   declared `type(graph), target` in the scope that owns it, and no procedure returns a pointer
   to a graph it constructed.
3. A `graph` is never returned by value from a procedure that assigned branches referencing its
   own local variables: intrinsic assignment copies `branch(2)` shallowly, so the copy would
   reference deallocated storage. Construction is `call g % declare()` in place; there is no
   graph-returning constructor.
4. Storage reclamation is by region, not by ownership. A cyclic reference structure admits no
   ownership discipline; the caller releases its storage as a whole. The kernel names no storage
   type, so this stays outside the ontology (§9).

## 4. Can sharing be represented?

Yes (test B). Both branches of `a` are assigned `known_branch(b)`, and three propositions hold:
`associated(a%branch(1)%known, a%branch(2)%known)`;
`a%branch(1)%known % same_as(a%branch(2)%known)`; and a mutation of `b` is observed on both
branches, which excludes a copy.

## 5. Can cycles be represented?

Yes (test C). `a -> b -> a`, verified by identity (`same_as`) and by reference
(`associated(a%branch(1)%known%branch(1)%known, a)`). The structure is a graph, not a tree.

## 6. Is mutation required for cycles?

**Required, and it is construction mutation.** Each of the two pointer assignments requires the
address of the other graph, so no order exists in which both are established when their graph is
created; the second assignment writes into an existing object. Fortran provides no recursive
binding form and no laziness. Mutation *after publication* is not required: in test C the cycle
closes before either graph is read.

## 7. Recommended realization model

**Construction-mutable, publication-immutable.** Mutation is confined to construction; after
publication the branch states are fixed and identity is stable. Not implemented — §10 says do not
legislate — but the enforcement point is a single `seal` operation shaped like `declare`.

| criterion | A. mutable | B. persistent (new identity per change) |
|---|---|---|
| sharing | holds — one graph, many references (test B) | **degrades** — untouched substructure stays shared, but every path to a changed graph is rebuilt, and until it is, predecessors reference the previous graph (test J) |
| cycles | holds (test C) | **no fixpoint** — rebuilding `q` requires `p1`, which requires `q1` (test J) |
| caching on identity | **unsound after publication** — the view value changes, identity does not (test I) | sound — a changed graph is a new key |
| immutable views | unsound while mutation is admitted | sound |
| CSR compilation | valid only against fixed branch states | sound; recompiled per version |
| checkpoint/restart | *argued, not tested* — one version to serialize; restart re-mints identity | *argued, not tested* — all versions persist |
| adjoint snapshots | *argued, not tested* — snapshots require explicit copies | *argued, not tested* — versions are the snapshots |

A is sound where structure is constructed; B is sound where structure is read. The boundary
between the two is publication, which is the recommendation. Tests I and J measure the two costs:
I the cost of mutating after publication, J the cost of never mutating.

## 8. How values are kept outside the kernel

The kernel declares no `carries_value`, no `payload`, no `value()`. Its code lines contain none
of `real`, `field`, `coordinate`, `symbol`, `payload`, `arena`, `handle`, `node`, `storage`. A
graph carries `branch(2)` and an identity.

Attributes are bound in `graph_views.f90`, in a partial map keyed by identity:

```
attribute_map : identity -> (number, symbol, index)
```

```fortran
call m % bind(two,   number = 2.0_dp)
call m % bind(three, number = 3.0_dp)
call m % bind(plus,  symbol = '+')
```

`(2 + 3) * 4 = 20` evaluates from that map (test D); the operands are `(NULL,NULL)` graphs whose
branch states encode nothing of 2 or 3. A further column adds no character to the graph type.

`id()` returns the token by value rather than requiring callers to key on
`type(graph), pointer`. A token is serializable and independent of the referent's lifetime, which
the checkpoint/restart row requires; a pointer key would bind every attribute map to graph
lifetime.

## 9. How CSR compiles without changing the ontology

`compile_csr` traverses the graph and returns

```fortran
type :: csr
   integer, allocatable :: rowptr(:), colidx(:)
end type csr
```

which does not extend `graph`, declares no branch, and is not named by the kernel: the module
dependency is `graph_views -> fractal_graph` and not the reverse. Test 11 compiles the
finite-volume face sequence to `rowptr = [1,2,4,5]`, `colidx = [2,1,3,2]`, then verifies that
every tuple of the relation view appears in the compiled rows. Semantic representation and
compiled representation are related by a traversal, not by a type.

One sequence carries four views and no additional type: `relation_view`, `residual`
(`R(Q) = 0` for `Q` constant), `compile_csr`, and — the same pair-and-sequence encoding —
`evaluate`. §12 holds: relational, computational, ordinary and mesh are views, not types.

## 10. Public kernel API

| symbol | why it is not derived |
|---|---|
| `graph`, `graph_branch` | the ontology and its recurring component |
| `graph % branch(2)` (public) | the self-similar structure; traversal is direct |
| `graph_branch % status` (public) | the state, explicit, never inferred from association |
| `graph_branch % known` (public) | the reference to `graph` |
| `GRAPH_NULL`, `GRAPH_UNKNOWN`, `GRAPH_KNOWN` | the three admissible branch states |
| `null_branch()`, `unknown_branch()`, `known_branch(g)` | the only admitted branch constructions; `known_branch` is where §4's validity condition is enforced |
| `declare()` | identity is minted, not chosen |
| `id()` | exports identity so attributes can be bound outside (§8) |
| `same_as(other)` | *derivable* from `id`; retained because §5 and §13 state the comparison in it |

Eight module-level names and three type-bound procedures. Evaluation, relations, residuals, CSR,
numbers, symbols and indices are outside. A graph label is outside too: a label is an attribute,
so `declare` takes no argument and the kernel declares no `label`.

Two refusals are kernel laws: `twice` (*graph identity is assigned once*) and `undeclared`
(*KNOWN requires a graph with assigned identity*, enforcing §4). Six are view laws, two of them
because §4 is easy to violate one level up: a view answering "absent" for a graph outside the
map's domain, or for an UNKNOWN member, would merge UNKNOWN into NULL. `noindex` and
`unknownmember` refuse instead.

**One admitted weakness.** §2 requires `g % branch(i) % known % branch(j)` to traverse, so
`status` and `known` are public, so `g % branch(1) % status = GRAPH_KNOWN` can be assigned
directly, producing a KNOWN branch with a disassociated reference. The three branch constructors
are the only admitted constructions and cannot produce that state. The alternative — private
components behind `status(side)` and `known(side)` — is the interface §2 prohibits.

## §16. Visual acceptance

```fortran
type :: graph
   type(graph_branch)   :: branch(2)
   type(token), private :: identity
end type graph
```

and `graph_branch` declares `type(graph), pointer :: known`. Therefore graph -> branch -> graph
is visible in the declarations themselves. `(NULL,NULL)` is the default initialization of the
type, so §3's "an atom is a graph" requires no code.
