# the pipeline: implementation, and what it absorbs

The pipeline is the reviewer-named home for a simple fact the engine
keeps rediscovering: every evaluation is stages joined by "computed
from", and every derivative, checkpoint and parallel schedule is a
way of walking those stages. `derivative-walks-plan.md` motivates it;
this plan says what gets built and what existing machinery turns out
to have been a pipeline all along. Written 2026-08-02.

## the object

One concrete module, `class_graph_pipeline`, no contract change. A
pipeline IS an operation - it takes a graph and data and returns
data - so it extends the existing vertex-field leaf, and the contract
never learns a new word.

    type, extends(graph_vertex_field_operation) :: graph_pipeline

       type(stored_graph) :: stages         ! what fed what: stage
                                            ! vertices, dependency
                                            ! edges

       type(edge_differential_operator),   allocatable :: edge_stages(:)
       type(vertex_differential_operator), allocatable :: vertex_stages(:)
       integer, allocatable :: itinerary(:) ! stage order, which list,
                                            ! which entry

    end type

The stage record is itself a stored graph, so the engine's own
machinery serves the pipeline: the dependency order is the visit
order, the reversal is the adjoint's itinerary, colouring finds the
stages that can run together, and partitioning a pipeline is
partitioning a graph.

Stages come in exactly the two built-in operator types, for the
reason the derivative plan gives: derivative walks only ever traverse
the linear shadow of a computation, and the linear shadow lives in
those two types. A nonlinear formula joins a pipeline through its
linearization - the coefficient arrays it freezes at the current
state - which is the obligation it already carries.

## the four walks

    value walk      visit stages in itinerary order; each stage's
                    apply feeds the next. The plain evaluation.

    tangent walk    the same itinerary, the same operators - a linear
                    stage is its own tangent - fed the perturbation
                    instead of the state.

    adjoint walk    the itinerary reversed, every stage with its
                    adjoint flag raised. One walk, whole gradient.

    segment walk    any of the above between two named stages.
                    Transposes factor stage by stage, so a segment is
                    a complete answer, not an approximation: the
                    sensitivity with respect to that segment's input.
                    Checkpointing and mid-chain sensitivities are
                    segment walks; nothing in the type may assume a
                    walk runs first stage to last.

All four are loops over the typed stage arrays. Work buffers are
allocated once per walk, not per stage - which is also where the
recorded copies-per-apply concern gets its fix for free: inside a
walk, stages hand each other work arrays, and fields appear only at
the pipeline's boundary.

## what the pipeline absorbs, audited

Each existing structure, its pipeline reading, and the action taken.

**The balance - absorbed first.** A balance is a two-stage pipeline:
edge terms, then the incidence step, plus a vertex source. It stays
as the convenient front door - models declare terms, not itineraries -
but its apply becomes "build the two-stage pipeline, walk it", and
the review's item about vertex-side terms lands here naturally: a
balance with a pressure-force term is just a pipeline with one more
vertex stage, and the balance type grows the `vertex_terms` array the
machinery review already called for. One implementation of walking,
not two.

**The operator chains - formalized, not absorbed.** The S/G/D
composition inside a differential operator is a pipeline in the
mathematical sense, and stays a fused kernel loop in the engineering
sense: no fields, no allocation, one component pass. The pipeline
composes OPERATORS; the operator composes STEPS. Two granularities,
each stated in the other's banner, both walked by the same
mathematics.

**The time march - the nesting target.** The integrator already
forms a chain and walks it forward, and its adjoint walks it
backward; that chain is a pipeline whose stages are instants and
whose dependency edges are the stencil, with the stencil coefficients
as edge coefficients. The old `chain` type (which `marcher` extends,
the wrong way round) is what the pipeline replaces in the new world:
when the gutting reaches the marchers, they hold a pipeline instead
of being a chain. Not built now; the shape is fixed now so the
migration lands on it.

**Solver iteration - one sweep is a pipeline.** A linear solver's
pass - precondition, operate, correct - is a pipeline walked once;
the iteration loop and its tolerance stay with the marcher, which
walks the pipeline repeatedly. The pipeline never loops; loops own
pipelines. That division keeps termination logic out of structure.

**The one-line law - the formalization prize.** Partition, operate,
assemble is itself three stages. A partitioned pipeline is the law
as a runnable object: partition stages at the head, operator stages
on the parts, assembly at the tail - and the exchange operation,
still unwritten, finally gets its natural slot: an ordinary stage
sitting before any stage that reads borrowed values. When exchange
is written, it is written as a pipeline stage, not as solver
plumbing.

**The reduction tail.** A pipeline ending in one number - an
objective - ends with a reduction. The reduction's four steps
already handle parts; the pipeline records which stage feeds the
reduction and the adjoint walk starts by broadcasting the seed back
through it. This is the reduction seat the derivative plan listed,
placed.

**Not absorbed, stated to prevent scope creep:** the reductions'
combine tree across parts (already correct, four steps, leave it);
the walks module (colouring and friends are single operations, not
compositions); the transforms (partitioner, assembler, coarsener,
refiner act on graphs, not through stage chains - they appear IN
pipelines, they are not made OF them).

## construction

Models never write itineraries. Three builders, in order of arrival:

    from a balance        the two-stage form, plus vertex terms
    from an operator list ordered stages, dependencies inferred from
                          which side each operator reads and writes
    from an integrator    the time chain, when the marcher migration
                          lands

## the checks

    equivalence     a balance evaluated directly and through its
                    pipeline agree to machine precision
    pairing         the pipeline's tangent and adjoint walks are
                    matrix transposes - the unit-field test, lifted
                    from operators to pipelines
    gradient        the adjoint-walk gradient of a small objective
                    equals the complex-step gradient entry for entry
    segments        an end-to-end adjoint equals the composition of
                    its two segment walks split at any stage
    conservation    a pipeline ending in the incidence step conserves
                    on a closed ring, whatever stages precede it
    cost            one adjoint walk per gradient, counted

## order of work

1. The type, the value walk, and the balance equivalence check.
2. Tangent and adjoint walks with the pairing and gradient checks;
   segment walks with the split check.
3. The balance grows vertex terms and a source field (the machinery
   review's item 2), implemented as pipeline stages.
4. Hessian-vector products by the complex road over the adjoint walk
   (derivative plan, step 3).
5. The time nesting and the exchange stage - both gated on the mesh
   migration and the marcher gutting, both with their shapes fixed
   by this plan.

Steps 1 through 4 need no contract change and no mesh, and run on
the chains and rings the suite already uses.
