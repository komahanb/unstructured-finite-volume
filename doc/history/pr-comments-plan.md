# plan for the review of class_graph_differential_operator

Fourteen new comments, one file-level note, and three threads parked
earlier. Sorted into what gets addressed and what deliberately does
not, with the mathematical answer worked out here where the comment
asked a mathematical question. Written 2026-08-02.

## items that need addressing

### 1. "fold" is not a mathematical word - rename it

Two comments, one angry emoji, and the point stands. "Fold" is a
programming word I imported; the operation it names - the signed sum
over the edges at a vertex, divided by the measure - has a familiar
mathematical name, and the named layer already uses it: divergence.

Action: the D step becomes **the incidence step** - resolved purer than
this plan first proposed. Divergence is calculus's word and stays on
the named constructor; incidence is the graph's own word, already in
the architecture note ("incidence reduces z_e onto vertices"), and
names the step in the engine's register. `fold_step` becomes
`incidence_step`, with every comment, label and document following.

### 2. the division by two is systematic, and the comment will say why

The question answers itself in the second comment: the average is an
integral in disguise. Walk along one edge from tail to head; the
integral of a straight-line q over that edge, divided by the edge's
length, is exactly (q_i + q_j)/2. So the S step is **the edge
average** - the mean value of q along the edge - and the one-sided
variant is the same integral evaluated by the value at one end. Both
are the two simplest ways of averaging along an edge; neither number
is arbitrary.

Action: the S banner states this, S is called the edge average step,
and the guide gains the same sentence.

### 3. c is used before it is defined

The header's order table writes c q before saying what c is.

Action: one line above the table - "with c the coefficient the
operator carries" - and the table reads clean.

### 4. tabulate both sides, and state the recurrences explicitly

The vertex side has a table; the edge side is prose. And the comment
asks for the recurrence itself, which is the cleanest statement in
the whole design:

    edge(0)   = S                the edge average
    vertex(0) = c q              the term itself

    edge(n)   = G of vertex(n-1)     n >= 1
    vertex(n) = D of edge(n-1)       n >= 1

with the closed forms falling out by parity:

    vertex(2k)   = (D G)^k           even orders
    vertex(2k+1) = (D G)^k D S       odd orders

Action: both tables and both recurrences go into the header,
mirrored in the guide.

### 5. "ring" is not in the taxonomy - define it, with a drawing

Fair. A ring is the set of vertices a walk of exactly r edges can
reach and a shorter walk cannot:

         2   2   2
       2   1   1   2          ring 1 of v: the neighbours
       2   1 v 1   2          ring 2 of v: the neighbours'
       2   1   1   2                       neighbours, minus
         2   2   2                         what ring 1 took

The engine already computes rings: they are the level sets of the
depth walk. Action: the definition and drawing go in where "ring"
first appears, and the word joins the working vocabulary.

### 6. the uniform chain is a calibration bed, not an assumption

The comment worries that uniformity is baked in. It is not: every
parameter is per-entity - spacings per edge, measures per vertex -
and the operators are defined on any graph. The uniform chain appears
only where exactness is CLAIMED, because that is where the discrete
formulas coincide with the calculus ones, which is what makes the
claims checkable by hand.

Action: the exactness paragraph says this in one sentence: defined
everywhere, calibrated on the uniform chain.

### 7. odd and even orders have different adjoints - say so

The parity structure is clean and deserves its paragraph:

    even orders    vertex(2k) = (D G)^k is its own adjoint whenever
                   the same weights sit on its steps - the transpose
                   of a power is the power of the transpose, and
                   each D G is self-adjoint

    odd orders     the adjoint reverses the sampled end: the
                   transposed chain ends with the average step's
                   transpose, which returns each edge value to the
                   end it was read from. A one-sided operator's
                   adjoint is one-sided the other way - the
                   downstream walk of an upstream sample

Action: this paragraph joins the adjoint banner, and the guide's
adjoint section gains the parity statement.

### 8. the ladder of orders, with what each one models

The comment asks the analogy be extended as far as it goes. The
ladder, each rung named by the physics that lives on it - names that
belong to model classes, listed here only as motivation:

    order 0    the term itself         storage, reaction, mass
    order 1    first derivative        transport along a flow
    order 2    second derivative       diffusion, conduction,
                                       viscosity, pressure fields
    order 3    third derivative        dispersion - waves whose speed
                                       depends on their length
    order 4    fourth derivative       bending - beams, plates, and
                                       the interface physics of phase
                                       separation
    order 6    sixth derivative        pattern-forming films and
                                       crystals

Action: the ladder goes into the closing banner of the operator file.

### 9. must an adjoint walk run end to end? no - and this is worth a section

The deepest comment of the batch. The answer: the adjoint of a
composition factors stage by stage,

    (C B A) transposed = A* B* C*

so a reverse walk may stop at any stage - yielding the sensitivity
of the answer WITH RESPECT TO that intermediate stage - or start
from one, or run any contiguous segment. Nothing in the mathematics
requires end-to-end; end-to-end is just the segment that happens to
reach the input. The segmented walks are not a curiosity; they are:

    checkpointing        long time marches store a few instants and
                         re-walk segments between them, trading
                         compute for memory
    stage sensitivities  the derivative of the answer with respect
                         to an edge sample or an intermediate field,
                         read off mid-walk
    preconditioning      partial transposes as building blocks

The current adjoint flag runs the whole chain because an operator IS
one chain; no limiting assumption is introduced there. The place
segments belong is the pipeline (item 10), whose walks are over
stages by construction - a segment is just a sub-walk, and the plan
document now says so.

Action: a paragraph in the adjoint banner stating the factorization,
and a section in `derivative-walks-plan.md` making segment walks an
explicit requirement of the pipeline design.

### 10. class_graph_pipeline - yes, and it already has a design

The file-level comment proposes it, and it is precisely the
"evaluation graph" of `derivative-walks-plan.md`: stages, what fed
what, walked forward for values and tangents, backward for adjoints,
in segments for checkpoints. The proposed name is better than mine -
pipeline says what it is in one word.

Action: the plan document adopts the name `graph_pipeline` in module
`class_graph_pipeline`, and the pipeline becomes the home of the
segment walks of item 9. Building it remains the plan's step one,
after these comment fixes land.

## items that do not need addressing, and why

### the scalar-plus-array parameter pairs (line 160)

The comment says it can live with the duality for now and asks that
purer alternatives be sought. Agreed on both halves. Recorded in
`graph-deferred-concerns.md` with the candidate worth exploring
when the time comes: a constant field - a graph data type whose
storage is one number but whose contract answers per entity - would
collapse the pair into one member with no storage cost, at the price
of a level of indirection in the hot loops. That trade is measurable,
so it should be measured, and the mesh migration will provide the
mesh to measure it on. No action now.

### the three parked threads on the contract

- **why should a graph know where it came from** (the part relation):
  parked by the earlier ruling - typed procedures beat stringly data
  for compile-time safety. Stands.
- **the reverse lookup is an efficiency concern**: agreed and
  recorded; nothing hot touches it.
- **defined_on_graph as an admission about robustness**: the comment
  wording was fixed and the coarsener no longer refuses; whether the
  gate should exist at all is a contract question deferred until a
  caller actually suffers it.

### compression, "with an obsessive compulsive graph disorder"

The instinct is right and the bar is known: consolidation is judged
by net lines deleted, and a pass over the forward/reversed kernel
pairs and the two applies is plausible. But the file is still
moving - the divergence rename and the comment work above touch it
throughout - and compressing a moving file wastes the compression.
Scheduled for after this batch lands and the file goes quiet.

## order of work

1. The divergence rename, everywhere at once - code, tests, guide,
   documents - so the word "fold" is gone in one commit.
2. The comment mathematics: items 2 through 8, one commit, comments
   only, verified by a clean build and the unchanged 161 checks.
3. The plan-document updates: segment walks and the pipeline name.
4. Replies on the pull request, resolving what is fixed, leaving
   open what is parked.
