# the graph mathematics machinery, reviewed

A walk through everything the graph engine now proves it can do, one
defect found and fixed during this review, and the straightening-out
list before the machinery carries a flow solver. Written 2026-08-02,
at commit `5ba3bed` plus the fix below.

## what stands

One contract file and thirteen concrete modules. The contract,
`abstract_graph.f90`, declares 91 procedures across 22 types and owns
no storage and no algorithms. The concrete layer behind it:

    structure      class_graph                the stored graph
    sets           class_graph_support        vertex and edge sets
    data           class_graph_field          values on sets, any of
                   class_graph_functional     five kinds, interleaved
                                              components
    one number     class_graph_reduction      sum, integral, average,
                                              min, max, norm, count,
                                              all, any - four steps,
                                              so parts combine exactly
    frame          class_graph_partitioner    whole -> parts
                   class_graph_assembler      parts -> whole
    detail         class_graph_coarsener      fine -> coarse
                   class_graph_refiner        coarse -> fine
    calculus       class_graph_differential_operator
                                              any order, per-edge
                                              coefficients, adjoint
                                              as the reverse walk
    composition    class_graph_balance        vertex source plus
                                              edge terms reduced
                                              through incidence
    algorithms     class_graph_walk           colouring, visit order,
                                              components, depth

The theory guide, `doc/graph-differential-operators.pdf`, states the
conventions and proves the exactness claims; the contract suite
carries the same numbers as 156 checks; all seventeen suites pass.

## what the tests actually pin

The properties worth trusting because a check fails without them:

    the ordering rule        a value sits at (entry-1)*ncomp + c;
                             nothing else in the library can catch a
                             violation, so it has its own check
    the three identity laws  assemble(partition(G)) round-trips; every
                             cell owned exactly once under three cut
                             rules; a sum computed piecewise equals
                             the sum computed whole
    the reduction steps      parts averaging 2 and 7 combine to 4,
                             either arrival order
    conservation             a balance on a closed ring sums to zero
                             under three states and two terms; opening
                             the ring breaks the cancellation, so the
                             check has a failure mode
    the derivative numbers   a line slopes, a parabola gives 2, the
                             fourth power gives 24 at order 4, a
                             weighted incidence step gives the hand-computed 8
    components               a line and a parabola ride one field and
                             one application answers both slots
    adjoints                 forward and reverse walks pair to machine
                             precision, scalar and two-component; the
                             weighted order-2 operator is its own
                             adjoint
    curl                     a difference field cancels around the
                             square's border graph; a circulating
                             field does not
    complex step             the imaginary part of a functional
                             survives a reduction at 4e-20

## found in this review

**A per-edge coefficient was spent twice on one path.** The vertex
operator fed an edge field scales the samples by the coefficient
before its first incidence step. At order above one it then continues into the
deeper chain - and passed its coefficient array along, where the
innermost step applied it again. The comment said "already spent";
the code disagreed. Fixed by passing a never-allocated array, and
pinned by a hand-computed ring case that only passes with a single
application. No released result was wrong: no test or caller had
exercised that path.

## to straighten before CFD, in order

**1. The mesh migration.** Everything else on this list is small
beside it. The operators take their geometry as arrays - face area
and conductivity in the coefficients, centre distance in the
spacings, cell volume in the measures - and today those arrays are
filled by hand. The mesh (`class_mesh.f90`, 73 KB, which IS a graph
by inheritance the wrong way round) holds all of it. Until the
migration lands, a flow model can run on the engine only by doing its
own geometry bookkeeping. This is the recorded gutting task, and it
is the gate.

**2. The balance holds no vertex-side terms.** It carries edge terms
and one scalar source. Two consequences for a flow solver: the
pressure force - an incidence step with normal-weighted coefficients, which is a
vertex operator - cannot join a balance today; and a source that
varies cell to cell cannot either, because the source is one number.
The fix is small and additive: a `vertex_terms` array of the vertex
operator type beside the existing `edge_terms`, and a source field
option beside the scalar.

**3. No exchange between parts.** A partitioned run cannot refresh
its borrowed values; the design says exchange is an ordinary
operation, and nobody has written it. Serial CFD works without it;
parallel CFD does not exist until it lands.

**4. No inner product of two fields.** The reductions take one field.
Solvers stop on dot products and energy norms of pairs - the
architecture note lists the inner product among the reductions, and
it is the one rule not implemented. Small: one more rule with a
second-field argument, or a two-field variant of `reduce`.

**5. Copies per application.** Each `apply` builds a fresh field,
copies values in through `set_real_vector`, and copies again through
`allocate(source=...)` - roughly three traversals of a full field per
operator call, and a supplied buffer saves none of it because the
implementation deallocates and reallocates anyway. The contract
permits reuse; the implementation does not yet exploit it. Harmless
at test size, real at a million cells. Same family as the recorded
`get_data` copy.

## fences to state, not fix

Limits that are correct to have, so long as they are said out loud:

- **Boundary closure of composed orders.** Every average and
  difference step uses the stored boundary value where a head is
  missing - the boundary value of q, even when the intermediate field
  is no longer q. Orders one and two are right; order three and above
  are naive within two rings of a boundary. The suite tests deep
  interiors, which is where the claims hold.
- **The adjoint is of the linear part.** Boundary values make an
  operator affine; the reverse walk transposes the linear part, which
  is the standard meaning and the one a sensitivity solver wants.
- **Per-edge coefficients are single numbers.** The dissipation term
  of a wave-speed flux couples components through a small matrix per
  edge. That formula is a concrete physics class on the edge leaf,
  which can compute anything; the engine's built-in operators offer
  scalars.
- **One metric for all components.** Coefficients, spacings, measures
  and boundary values are shared across a field's components; a
  component needing its own gets its own operator instance.
- **Odd orders above one** mix one-sided and centred steps and match
  a stencil, not a calculus formula.

The previously recorded concerns - the `get_data` copy, the linear
reverse lookup, the unweighted breadth-first cut, one part per
assembly call, complex reductions being sums, character fields not
reducing, and operations not declaring their inputs by name - stand
unchanged in `graph-deferred-concerns.md`.

## the shape of a flow solver on this engine

For the record, the delegation as it now stands. A model brings:

    convection      order 1, coefficients = mass flux per face,
                    refreshed each evaluation
    viscosity       order 2, coefficients = viscosity times area
                    over distance
    pressure force  interpolation reduced with normal-weighted
                    coefficients, one call per direction
    continuity      the same reduction read as a constraint
    wave-speed
    fluxes,
    limiters        concrete classes on the edge leaf
    boundaries      per-edge boundary values
    time            the existing integrators

and writes no loop over cells or faces. Items 1 and 2 above are what
stand between this table and running it.
