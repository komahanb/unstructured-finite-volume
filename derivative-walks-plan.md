# derivatives of any order, as walks over walks

A design for computing gradients, Hessians, and higher derivatives of
anything the graph engine evaluates, using the two directions a walk
can go. Written 2026-08-02, as a proposal. Nothing here is built yet.

## the idea in one picture

Every evaluation the engine performs is already a walk over a graph.
A residual walks edges and folds; a time march walks a chain of
instants; a functional walks a support and reduces. Draw what fed
what and the evaluation IS a graph:

    q ---> edge samples ---> folded balance ---> J
              (operator)         (fold)        (reduction)

Differentiation is then not a new kind of mathematics. It is a
direction of travel on that graph, and there are exactly two:

    the tangent walk    same graph, same direction, carrying
                        "how much does this stage move when the
                        input moves" forward

    the adjoint walk    the same graph reversed, carrying
                        "how much does the answer care about this
                        stage" backward

One reverse walk yields the whole gradient, however many inputs
there are. One forward walk yields one directional derivative. That
asymmetry is the entire economics of sensitivity analysis, and it
falls straight out of the picture: the reverse walk starts at the one
answer, the forward walk starts at one input.

This is the time-integrator construction, generalized. The integrator
already says it in its own comment: it forms the chain at the
scheme's stencil depth and walks it forward, and a discrete adjoint
walks the same chain backward. The proposal below makes that sentence
the design for everything, not just for time.

## the crux: derivatives only ever walk the linear shadow

Here is the observation that makes the design small.

An evaluation may contain anything - wave-speed formulas, limiters,
coefficients recomputed from the state. But its DERIVATIVE walks
never touch those formulas directly. They walk the linearization: the
frozen-coefficient linear operator each nonlinear piece produces at
the current state. That tier already exists - it is the
lagged-coefficient path the suite exercises - and its two concrete
types are the edge and vertex differential operators.

So the derivative walks only ever traverse operations of exactly two
concrete types. Two consequences:

    no heterogeneity wall    a plain array holds the stages, because
                             they are one type per side - the same
                             resolution the balance and the fluxes
                             already use

    no contract change       tangent and adjoint of the built-in
                             operators exist today: the tangent of a
                             linear operator is itself, and the
                             adjoint is the flag, proven equal to the
                             matrix transpose at every order tested

A nonlinear piece participates in exactly one way: it owes its
linearization. Given the state, it fills coefficient arrays - the
mass flux, the wave speeds, the limiter slopes - and hands back
built-in operators. That is the whole obligation, and it is the
obligation the freeze-and-thaw pattern in the old assembler already
imposed. Physics delegates; the walks are engine.

## the evaluation graph, as an object

The one new concrete thing: a type that records what fed what.

    type :: evaluation_graph
       type(stored_graph)                              :: stages
       type(edge_differential_operator),   allocatable :: edge_stages(:)
       type(vertex_differential_operator), allocatable :: vertex_stages(:)
       ...which stage is which, in walk order
    end type

Its vertices are the stages of an evaluation - state in, samples,
folds, functional out. Its edges say "computed from". A model builds
it from what it already declares: the balance's term lists, the
integrator's stencil, the reduction at the end. And because it is a
stored_graph underneath, everything the engine does to graphs applies
to it: it partitions, it walks, its dependency order is the visit
order, and its reversal is the adjoint's itinerary.

The two walks over it:

    tangent(eg, x, dx)  ->  dJ      visit stages in order, each
                                    applying its operator to the
                                    incoming tangent

    adjoint(eg, x, dJ)  ->  dx      visit stages in reverse, each
                                    applying its operator with the
                                    adjoint flag raised

Both walks are loops over typed arrays, like every walk in the
engine. Neither allocates per stage.

## any order, by walks over walks

A walk is itself an evaluation, so it has its own evaluation graph,
and the two walks apply to it again. Order grows by composition,
exactly as operator order grew by repeating steps:

    gradient          one adjoint walk                        A
    Hessian times v   a tangent walk over the adjoint walk    T A
    full Hessian      one T A per direction, n of them
    third order       another letter                          T T A

The Hessian needs one thing the first-order walks do not: the
derivative OF the linearization - how the frozen coefficients move
when the state moves. Two ways to get it, and the design supports
both:

    the honest seat     a physics class that can differentiate its
                        own coefficient fill provides it, and the
                        tangent walk uses it

    the complex road    run the ADJOINT walk in complex arithmetic
                        with a tangent perturbation riding the
                        imaginary part. The imaginary part of the
                        gradient is the Hessian-vector product,
                        exact to machine precision, with only
                        first-order seats written by hand.

The complex road is why the fields were taught to carry complex
values and why the functional keeps its imaginary part - the
machinery for it exists end to end, checked at 4e-20. A physics
class that provides nothing beyond its linearization still gets
exact Hessians.

## time is the same structure, nested

The unsteady case adds nothing new. The march is an evaluation graph
whose stages are instants and whose edges are the stencil - the
time-integrator chain, verbatim. Each instant's stage operator is
itself an evaluation graph over space. The unsteady adjoint is then
the reverse walk of the march, where visiting a stage means reverse
walking its spatial graph:

    forward:   instant 1 ---> instant 2 ---> ... ---> instant N ---> J
    adjoint:   instant N <--- ... <--- instant 2 <--- instant 1

    each box, opened, is a spatial evaluation graph
    walked in the matching direction

Graphs whose stages hold graphs - the same move as the border graph
for curl, repeated one level up. The stencil coefficients the
integrator already exposes are the per-edge coefficients of the time
chain, and the backward sweep's future-step coupling is the reversed
chain's incidence. Nothing about time is special; it is the outermost
graph.

## what each piece owes, in full

    a linear operator      nothing. it is its own tangent, and its
                           adjoint is the flag.

    a nonlinear formula    its linearization: fill coefficient
                           arrays at the given state, return built-in
                           operators. (already the lagged tier.)

    a reduction            its tangent is the same reduction on the
                           tangent field; its adjoint broadcasts the
                           seed to the support. two one-line seats.

    a time scheme          its stencil, which it already exposes.

    the engine             the evaluation graph type, the two walks,
                           and the composition of letters.

## the checks that would pin it

Every claim above has a machine-checkable form, in the style the
suite already speaks:

    pairing        tangent and adjoint walks of one evaluation graph
                   pair as transposes - the matrix test, one level up
    cross-check    the adjoint-walk gradient equals the complex-step
                   gradient, entry for entry
    symmetry       Hessian-vector products satisfy
                   sum((Hv) w) = sum(v (Hw))
    cost           the gradient costs one reverse walk regardless of
                   the number of inputs - counted, not assumed
    nesting        the unsteady adjoint on a two-instant march equals
                   the hand-unrolled derivative

## order of work

1. The evaluation graph type and the two first-order walks, with the
   pairing and cross-check tests. Small: the walks are loops, the
   operators exist.
2. The reduction seats (two one-line procedures on the reduction).
3. Hessian-vector by the complex road over the adjoint walk, with
   the symmetry test.
4. The time nesting, once the mesh migration lands and a real
   unsteady residual runs on the engine.

Steps 1 through 3 need no contract change and no mesh; they can be
proven on the same chains and rings the contract suite already uses.
