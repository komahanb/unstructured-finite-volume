# Topology ownership and traversal

For a stored directed graph, the endpoint arrays, compressed incidence
indices, counts, tags, carrier identities and partition relation describe
one topology. Its constructor validates the supplied extents, tail
members and partition metadata before constructing any derived index.
Global index maps are injective, and a carrier identity cannot acquire
different counts from the part and whole descriptions.
These components are private. Replacing a complete graph value is allowed;
changing one of its arrays independently is not.

The ordinary directed interpretation retains its boundary convention:
an out-of-range head denotes an edge without a head. Every tail belongs
to the vertex carrier. This convention does not impose a degree limit on
the general relation representation.

## Read access

Scalar graph queries return values. The existing neighbourhood queries
return owned arrays. A whole incoming traversal can use
`g % read_incoming(reader)`, which calls
`reader(offsets, indices, sources)` once. The edge indices for row `v`
are `indices(offsets(v):offsets(v+1)-1)`; `sources(e)` is the tail of
edge `e` in the graph's current orientation.

The callback arguments are integer arrays with `intent(in)` and without
`target` or `pointer`. A conforming callback cannot modify them or retain
a pointer to them. The graph selects orientation before the callback,
without copying arrays or allocating a neighbourhood. The stencil uses
this interface for the same ordered row sums formerly computed through
public graph components. A default empty graph supplies offsets `[1]`
and empty index and source arrays.

For a binary relation `P` contained in `A x B`, `image_view(a)` and
`preimage_view(b)` return an `integer_fibre`. Its private reference
provides:

- `num_members()`: the fibre cardinality, in constant time.
- `member(i)`: the member at local position `i`, in constant time, with
  an explicit index check.
- `read(reader)`: one callback with an `intent(in)` integer array and
  no allocation. Its argument has neither `target` nor `pointer`.
- `read(reader, context)`: the same protected array read with explicit
  per-call state. The callback receives `members` and an `intent(inout)`
  polymorphic `context` owned by the caller.
- `values()`: an independent, allocated copy of the members.

Default fibres and fibres of absent members are empty. CSR construction
and tuple ordering are unchanged; source and target indices remain
consistent because the returned fibre cannot write either index.
The allocating relation queries `image` and `preimage` retain their
interfaces. Existing relation algorithms traverse each fibre through a
module procedure with explicit local context. This avoids a checked
scalar call per edge, captured-procedure overhead and module-global
mutable traversal state.
Test consumers also use the new readers; the former writable pointer
result is removed.

## Lifetime and transpose

A stored graph and a CSR relation own their arrays. An `integer_fibre`
borrows its array: its source must have `target`, remain allocated, and
not be reassigned for the duration of every use of the fibre. Copying a
fibre copies its reference and does not extend that lifetime. Constructors
of other binary relation representations can use `integer_fibre(members)`
under the same requirement. Use `values()` when a result must outlive its
source. Privacy enforces access; it does not manage borrowed lifetimes.

`transpose_of(relation)` is a constant-cost borrowed view. Its base must
outlive the view, and transposed fibres obey the same lifetime rule.
`stored_directed_graph % transpose()` remains an owning copy with reversed
orientation and the same carrier identities. It has storage cost
proportional to the graph. Sharing that storage requires a separate
ownership change; this interface does not claim constant-cost graph
transposition. Graphs with headless edges still refuse directed transpose.

## Counted ownership of immutable storage

`util_counted_storage` is the primitive for jointly owned immutable
storage. A cell (`counted_storage`, extended by its owner with the
stored contents and a `clear` procedure) records a version and the set
of live binding serials. A `counted_reference` is bound to a cell by
`acquire(template)` or by defined assignment, each binding with its
own serial; finalization or assignment over it removes that serial.
The last removal clears the cell's contents; the cell object is
retained and reused for the next acquisition of the same dynamic type,
with a new version, so a reference whose binding is stale addresses a
cell that exists, reports `live() = .false.`, returns no storage and
releases nothing. `num_owners()` is the number of live bindings.

The contract was measured for GNU Fortran 15.2 at `-O0` and `-O3`
(`artifacts/remaining-work-2026-09-11/r04/mechanisms/`, the same
program is the suite's `check_counted_storage`):

| Copy mechanism | Result |
|---|---|
| intrinsic assignment: scalar, container, nested container, allocated allocatable scalar | owner; released with the value |
| element and whole-array assignment; array constructor into an allocated array of references | owner |
| `intent(out)` dummy, block scope end, `move_alloc`, `value` dummy | released or transferred correctly |
| function result assigned to a scalar, a container or a component; `allocate(source=function())` | owner |
| `allocate(mold=)` followed by a dispatching copy of the same dynamic type | owner |
| structure constructor with a live reference component | the constructed object owns; the source keeps reading while the copy lives, and is released with it |
| `allocate(source=variable)`, scalar, container, `class(base)`, `class(*)` | not an owner: a bitwise copy with its twin's serial; whichever twin is finalized first releases the binding, the other is not live afterwards |
| polymorphic intrinsic assignment `class(base), allocatable :: b; b = a` | the same bitwise copy, whether or not `b` was allocated |
| reallocating `arr = [arr, x]` of a type containing the reference | gfortran 15.2 runtime bounds failure |
| assignment to an unallocated allocatable scalar of the type | gfortran 15.2 segmentation fault: allocate first |
| whole-array assignment, array constructor or array function result reaching the reference through two component levels | gfortran 15.2 internal compiler error; element assignment compiles and binds |

gfortran assigns a containing object by finalizing the destination
component in place and then calling the component's defined assignment
on a temporary that still contains the destination's former bytes, and
returns function results through bitwise-moved temporaries of which
only the last is finalized. Registration keyed by the reference's
address was measured and rejected: it leaks on every container
assignment and function result. Registration by binding serial, with
idempotent release, is correct for every mechanism that invokes
defined assignment and never reads released storage under the others.

The stored graph does not use the primitive. The driver copies
every rule and datum by `allocate(source=)`, minimizers copy
themselves and their actions the same way, and a rule is a polymorphic
operation that contains stencils and residual operators with graphs:
each such copy would be a bitwise duplicate of the graph's reference,
and its finalization would release the stored twin. A shared-cell
graph was implemented and passed the ownership laws in isolation, then
refused the first `constrain` of every demonstration through exactly
that path; the implementation is retained as
`artifacts/remaining-work-2026-09-11/r04/shared-cell-graph.patch`.
Moving graphs onto shared cells requires that no object containing a
graph is copied by `allocate(source=)`, structure constructor or
polymorphic assignment: the seventeen operation types and eight
minimizer types need a dispatching copy, and driver data must be
excluded or dispatched. That is the execution copy-semantics boundary
of R05 and R06. The module state (version and serial counters, the
recycled cells) is not synchronised across threads.

## Verification

`test/graph-topology-ownership/run.sh` checks the access boundary with
compiler refusals, exercises sparse fibres and transpose/incidence laws,
and counts allocations during traversal. The suite is part of
`./verify.sh`. Existing partition tests continue to refuse transfer
descriptions with incorrect identities, counts or mapped indices.

## Measured traversal cost

Seven alternating runs against commit `8d59154`, with processor affinity
fixed for each comparison, gave these median execution times. The builds
used GNU Fortran 15 with optimisation and bounds checks. Every numerical
or structural result matched its baseline.

| Consumer | Baseline seconds | Protected access seconds |
|---|---:|---:|
| Stencil, forward | 0.745259 | 0.728233 |
| Stencil, transpose | 0.764293 | 0.746573 |
| Sparse unsuccessful reachability | 0.309707 | 0.317108 |
| Dense unsuccessful reachability | 0.341709 | 0.369802 |
| Topological ordering | 0.445018 | 0.454466 |

The stencil comparison applies one million coefficients on 200,000 rows
120 times in each orientation. Its 400,000 final values are byte-identical
to the baseline and equal independent edge sums. The graph algorithm
comparisons use 500 unsuccessful queries on 6,001 vertices and 36,000
edges, 100 unsuccessful queries on 513 vertices and 261,632 edges, and
200 orderings on a 2,000-vertex, 23,922-edge DAG. The last vertex is isolated
in each reachability case; the DAG's required order is `1,...,2000`.
These are focused local measurements, not a whole-application speed
guarantee.

Scalar `member(i)` validates each access and has more call overhead than
a raw array subscript. Bulk readers avoid that per-member call; the graph
algorithms use the explicit-context form. The standalone benchmark reports
scalar, stateless callback, explicit-context callback and complete incoming
traversal separately. All four traversal paths must perform zero measured
allocations; owning copies intentionally allocate.

For 20,000 degree-six fibres traversed 1,000 times, the isolated benchmark
measured scalar access at 2.93 times baseline cost, stateless fibre reads
at 1.59 times, explicit-context fibre reads at 1.48 times, and complete
incoming traversal at 1.00 times. Thus privacy has a measurable per-fibre
cost even though the consumer measurements above remain within 8.3% of
baseline. This distinction is part of the access contract's measured
limits, not a claim that every read is equally fast.
