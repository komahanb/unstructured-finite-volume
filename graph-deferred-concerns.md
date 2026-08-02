# deferred concerns from the graph rewrite

Things noticed while building the new graph classes that are real but
were not worth stopping for. Nothing here is a bug in what has been
written; each is either a cost we chose to pay for now, or a place
where the locked contract makes a thing awkward. Recorded so they get
reviewed on purpose rather than discovered.

Written as the concretions landed, in the order they came up.

## performance

**get_data returns a copy.** The contract declares the result
`class(graph_data), allocatable, intent(out)`, so `get_data` has to
allocate a fresh object and source-copy the whole field, values
included. On a million cells that is eight megabytes for cell volumes
alone. It cannot be fixed inside the frozen signature.

The convention that keeps it harmless is written into the header of
`class_graph.f90`: an operation fetches what it needs once, at the top
of `apply`, never inside a loop over faces. That is a convention and
nothing enforces it. If profiling later shows this dominating, the
correct fix is a signature change, not a workaround.

**reverse_lookup is a linear scan.** `part_vertex_id(full_id, part_id)`
walks the part's `vfull` array looking for a whole-graph id. That is
O(cells in the part) per call. Assembly does not use it - assembly goes
the cheap way, local to full - so nothing hot touches it today. If
something later needs full-to-local in a loop, the part graph should
carry the inverse map rather than searching for it each time.

**The named sets allocate on every call.** `boundary_edges` and its
siblings build a fresh support each time they are called. That is by
design: they are meant to be called once when an operation begins, and
the split between these and the bare-id walking queries exists exactly
so the per-cell path allocates nothing. Worth a note only because the
cost is invisible at the call site.

**Breadth-first partitioning uses a plain queue and no weights.** It
grows each part outward from the first unclaimed cell until it has its
share. It does not balance by cell volume or by work, and it does not
try to minimise the cut beyond following the connections. It is
simple and deterministic, which is what the identity checks need; it
is not a serious partitioner. The old `partition_rcb` weighted by
coordinates, and that path has not been ported yet.

## where the contract makes something awkward

**One part cannot rebuild a whole.** `assemble_graph(part_graph,
full_graph)` is handed a single part, so with more than one part no
single call can restore the whole graph. What it does instead is fill
in that part's own share and leave the rest at zero, so adding the
answers across the parts rebuilds the whole. That is correct and it is
tested, but it means the law

    assemble( partition( G ) )  ==  G

is exact only when the graph was cut into one piece. For more pieces
the tested statement is the weaker, true one: the owned sets do not
overlap, their union is everything, and the pieces sum to the whole.

Worth revisiting if the assembler ever needs to return a genuinely
whole graph from a genuinely cut one. That would need a signature that
takes all the parts, which the contract does not have.

**A part graph resolves queries only for its own part.** `part_vertex_id`
returns zero for a part this graph is not. That falls out of
`full_vertex_id` taking no part argument - a part graph must therefore
record which part it is, and resolve queries only for that one. So a
replicated partitioner has to give each image its own part graph
rather than one object answering for every part, which is a change
from how the old code worked. Nothing depends on the old behaviour
yet.

**Complex reductions are sums only.** Minimum, maximum and norm are
real, because complex numbers do not order. Average is real too. A
complex-step objective is a weighted sum, so the case that matters is
covered, but an averaged complex objective would need the work state
to carry a complex tally.

**Character fields cannot be reduced.** A field may hold names, and a
functional may hold a word, but no reduction rule produces one.
Concatenation and lexicographic extremum are the only candidates and
neither has a use here. The state would have to carry a deferred
length, which makes combine reconcile two lengths and finalize pick
one.

## still missing, by plan rather than oversight

**Nothing declares what data an operation needs.** The input argument
to `apply` is positional and undocumented, and computed fields do not
live on the graph. That is fine while operations are wired together in
code. It is not fine for a configuration-driven framework, which is
the stated direction, because a config file has no way to connect an
operation to its inputs by name. The right shape will be clearer once
several real operations exist.

**Geometric partitioning is not ported.** `partition_rcb` cut
on cell coordinates. The new partitioner has linear, breadth-first and
adopted. Geometric needs coordinates off the graph by name, which is
easy, but the recursive bisection has not been rewritten.

**Data exchange between parts has no concrete class yet.** It needs
none - it is an operation like any other, and the design note states this -
but nobody has written it. Until then a partitioned run cannot refresh
its borrowed values.

## what gutting the old graph actually involves

This is the one piece of the plan that is not done, and it is not done
because it is not a sequence of small green steps. Measured rather
than estimated:

Eleven modules use `interface_graph`:

    class_chain        class_fvm_field     class_csr
    class_gauss_seidel class_geometric_multigrid
    class_mesh         class_stored_graph  class_sor
    interface_linear_solver  interface_physics  interface_multigrid

And six types *are* graphs by inheritance:

    mesh          extends graph      73 KB, the whole geometry layer
    csr_matrix    extends digraph    24 KB, held by every solver
    physics       extends graph      so every flux and source is a graph
    chain         extends digraph
    marcher       extends chain      so every solver and integrator too
    stored_graph  extends graph

The architecture note specifies that data, transforms and operations
*hold* a graph rather than *be* one. Every line above is that
relationship the wrong way round. Inverting it touches
`interface_physics`, `interface_flux`, `interface_marcher`,
`class_chain` and everything below them - which is what
`class_assembler` (62 KB) and `class_partitioned_assembler` (41 KB)
are built on, which in turn is what all five solver interfaces hold.

So it is one large change where the suites go red in the middle and
only come back at the end. It also touches the adjoint and
complex-step machinery, which the new contract was widened to support
but which nothing has yet exercised end to end - nothing in the repo
extracts `aimag` anywhere.

Three decisions come before the work, not during it:

**Does the new graph carry the old algorithm surface while migrating?**
Carrying it keeps every suite green throughout and leaves the graph
temporarily fat, against the whole point of the rewrite. Not holding
it is cleaner and goes red for a while. `class_graph_walk.f90` now
gives colouring, visit order, components and depth a home as
operations, so the algorithms at least have somewhere to land either
way.

**Does `physics` invert to hold a graph, or collapse?** The edge terms
already collapsed into one concrete type of a given order, and
`class_graph_derivative.f90` shows that shape working. The same could be
done to `physics` wholesale, or it could keep its hierarchy and simply
stop being a graph.

**Is complex step exercised before or after?** It is provisioned
across 19 files and opt-in behind `COMPLEX=yes`, and the new
reduction's complex road is tested only in the contract suite. Proving
it against real physics before the gutting would mean the rewrite has
a working reference; after means the rewrite is judged on the real
tests but with a longer red stretch.

## build

**A clean build of the library emits 209 warnings that predate this
branch.** `interface_graph.f90` alone accounts for 31. The new files
add two, both unused-dummy warnings on procedures that return a
constant; `interface_physics.f90` contains the same pattern. A
cleanup pass belongs in its own change, not in this rewrite.
