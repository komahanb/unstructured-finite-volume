# dereplication

## philosophy

a vector is values on the graph's vertices. today every image holds all n
values of every vector and co_sum rebuilds all of them after every product -
the images are photocopies doing identical work. that was the right first
rung: it let the solver stay untouched while the system side learned to
partition, and it left behind proofs we now stand on.

the graph says what an image should actually hold: the values at vertices it
can see. its owned set, plus its halo - the ghosts, which the graph already
lists, and which the parallel suite now proves are exactly the vertices its
owned rows reach. nothing more is ever read.

    today                            after
    ┌───────────────┐                ┌──────────┐   halo   ┌──────────┐
    │ image 1: all n│                │ own | gh │ <------> │ own | gh │
    │ image 2: all n│      ──▶       └──────────┘          └──────────┘
    │ image 3: all n│                each image: its part of the graph,
    └───────────────┘                one row of ghosts around it

communication stops being "sum everything everywhere" and becomes what the
graph drew all along: one value per cut edge. the words moved per product ARE
the edge cut - which rcb already minimizes. partition quality and
communication volume become the same number.

the solvers never find out. cg is inner products and matvecs; the system owns
both. change what a vector means underneath and every algorithm above stands
unchanged - the same trick the frozen linearization played on the solver, now
played on storage.

## what stands ready (proven, on the shelf)

- the exchange sets: ghosts(k) = exactly the reach of the owned rows
  (parallel suite, 0 violations, 1/2/4 images)
- the address book: the partition is deterministic and replicated - every
  image holds every part's own_list and gh_list, so each image can compute
  where any ghost sits in its owner's frame with no communication at all
- the trio: gather / scatter / dot on the graph; matvec_rows on the matrix;
  the per-image block preconditioner

## the local frame

an image's world, in order: its owned dofs, then its ghosts.

    local vector:   [ owned dofs        | ghost dofs      ]
                      1 .. nown           nown+1 .. nloc

    the frame is graph data: dofs_of(owned(me)) ++ dofs_of(ghosts(me)),
    and the global-to-local map is its inverse. both live on the graph.

the local operator: the sub-digraph induced by the owned rows, far ends
renumbered into the frame. rows 1..nown, columns anywhere in 1..nloc.
(principal_submatrix generalizes: rows from one list, far ends through the
frame's map.)

invariants the frame keeps:

- elementwise arithmetic (the axpys inside cg) preserves ghost freshness:
  ghosts are copies, copies of both operands stay copies of the result,
  scalars from reductions are identical everywhere
- the one place freshness dies is a product's output: w has only owned rows.
  so the one seat that restores it is the next product's entry: exchange the
  input's halo, then dot the owned rows. one seat, no discipline spread
  through the code.

## the exchange

one coarray buffer per image, holding its owned values. a ghost's owner and
its position in the owner's frame are computed once from the replicated
lists. the exchange is each image reading its ghosts from the owners'
buffers - one-sided, neighbour to neighbour, no collective.

    image me:  buffer[me](1:nown) = x(1:nown)     post
               sync
               x(nown+j) = buffer[owner(j)](slot(j))   pull, j = 1..ngh

## stages

1. the frame - graph machinery: the local numbering of a part (owned then
   ghosts) and its inverse map; the matrix answers the induced local block.
   serial checks pin both entry by entry.

2. the vectors move in - the partitioned system reports the local length;
   source and residual assemble owned rows only; the product becomes
   exchange + owned-row dots; inner_product becomes the owned prefix dot
   plus one scalar reduction. the only co_sum left carries one scalar.

3. the proofs - the parallel suite solves the same problems and gathers the
   owned answers to compare against serial to machine precision; a memory
   assertion (nloc = nown + ngh, not n); a communication assertion (words
   per exchange = the ghost count = the cut's dofs).

4. on demand - block amg already lives on the owned block; writers gather at
   the door on image 1; transient marches and adjoints join when a transient
   parallel problem exists to demand them.

## books

the partitioned assembler sheds its full-length paths; the frozen-refusal
door stays until stage 4 lifts it. no graph machinery changes meaning -
stage 1 only adds the frame queries next to the lists they read. vectors
shrink from n to n/p + cut. the co_sum of whole vectors disappears.

## check

gates, in order: the frame pinned serially; distributed cg == serial cg to
machine precision on square-20/40 and box-36 at 1/2/4 images; block-amg
still cuts iterations; the halo-reach check keeps its zero; memory and
word-count assertions hold.
