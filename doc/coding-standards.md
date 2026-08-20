# coding standards

one rule sits above the rest: write for the next human. the code, the
comments, and the commit all get read by someone who was not there when you
wrote them. explain the way feynman explained — simple words, a concrete
picture, nothing invented.

## how we use graphs

see the problem as a graph first. a mesh is a graph (cells are vertices,
shared faces are edges). a sparse matrix is a graph with a weight on every
edge. the time steps of a march are a chain. the partition, the halo, the
coloring, the frame — all of it is vertices, edges, and who is next to whom.
when a new problem arrives, the first question is: what are the vertices,
what are the edges?

machinery lives on the graph, written once. every other class either IS a
graph (it extends one — the mesh, the matrix, the chain) or CARRIES one and
asks it questions (the assembler, the smoothers, the time integrator):

```
   solvers        ask: residual? product? dot?
      │           (never look inside)
   assemblers     answer on the whole field, or on one
      │           processor's slice + fetched neighbours
   graphs         know who is next to whom, who owns what,
                  and do all the sorting, ordering, coloring
```

the division of labor, precisely: STRUCTURE belongs to the graph — who
neighbours whom, who owns what, where each ghost value lives, in what order
to visit. COMMUNICATION belongs to the system — the sending, the receiving,
the syncs. a lookup table is structure even when only the parallel code
reads it; the wire that moves the values is not.

pass the graph, not its belongings. if a routine takes an index list, a
coordinate array, or a dof list that some graph already owns, hand it the
graph and let it ask. handing a class the graph's arrays one by one is how
duplicate bookkeeping is born.

two kernels do all the heavy lifting, and nobody re-rolls them:
`counting_sort` (count, prefix-sum, scatter — every adjacency, table, and
grouping is built through it) and `power_iteration` (iterated products and
dots — every spread measurement). if you are writing a count-then-fill loop
or a multiply-measure-normalize loop, stop: the kernel exists.

the smell that finds hidden graphs: any do-loop whose bounds encode reach,
cutoff, or window arithmetic (`k-p:k`, `if (j > p_k) cycle`, a hand-kept
cursor) is a graph query in hiding. replace the arithmetic with
`neighbours`, `in_neighbours`, an ordering, or a table the graph builds.

absorbing work into the graph is judged by what it deletes from the
callers. if the graph gains a method and the callers do not shrink, it was
a rename, not machinery.

## spacing

declarations stand in columns. inside one block, every `::` lines up:

```fortran
   ! wrong - the :: marks wander
   integer, allocatable :: ghost_dofs(:), inverse(:)
   integer :: j, g, v, p

   ! right - one column
   integer, allocatable :: ghost_dofs(:), inverse(:)
   integer              :: j, g, v, p
```

group the declarations: dummy arguments first (in call order), one blank
line, then locals. dummies carry their intent, aligned like everything else:

```fortran
   class(graph)        , intent(in)  :: this
   integer             , intent(in)  :: k
   integer, allocatable, intent(out) :: owner(:), slot(:)

   integer, allocatable :: ghost_dofs(:), inverse(:)
   integer              :: j, g, v, p
```

one statement per line. never stack with semicolons:

```fortran
   ! wrong
   allocate(a(n)); a = 0; allocate(b(n)); b = 0

   ! right
   allocate(a(n))
   allocate(b(n))
   a = 0.0_dp
   b = 0.0_dp
```

space around `%`: `this % grid % num_cells`, never `this%grid%num_cells`.
blank lines separate thoughts — before and after every loop nest, block
construct, and comment banner. a wall of statements with no air in it is
as unreadable as a wall of prose.

## the comment banner

every type definition and every procedure gets a banner, and the banner
carries a drawing when there is any structure to draw — boxes, arrows,
before/after, the ring of neighbours. the test of a good banner: a reader
understands what the code below does WITHOUT reading it.

```fortran
   !===================================================================!
   ! The halo exchange: one slab out, one frame in.
   !
   !    ┌─ image 1 ────┐          ┌─ image 2 ────┐
   !    │ own:   a b c │ ──a───▶  │ ghost: a     │    each arrow is one
   !    │ ghost: x     │ ◀───x──  │ own:   x y z │    cut edge's value -
   !    └──────────────┘          └──────────────┘    the traffic IS the cut
   !===================================================================!
```

what belongs in a banner: the picture, the reason the routine exists, and
any promise it makes (what it refuses, what it leaves alone, what must be
true before calling it). what does not: a restatement of the code line by
line, or a history of how it got here.

inline comments are for the one fact the code cannot say itself — a
constraint, a subtle invariant, a deliberate refusal. never narrate the
next line.

## names

use the formal name when the thing has one: `counting_sort`, `transpose`,
`frame`, `neighbours`, `degree`. a real name can be looked up; an invented
one cannot. if you find yourself coining a word for an identifier, stop and
find what the textbook calls it.

names say what a thing IS, not how it is spelled in some other layer's
vocabulary. a graph method speaks graph (`vertices`, `edges`, `owned`),
never the vocabulary of its callers (`cells`, `matrix rows`).

no abbreviation soup. `num_vertices`, not `nv` in a field; short names
(`i`, `v`, `e`, `p`) live only inside a loop where the banner already
said what they run over.

identifiers stay plain and standard; the COMMENTS may be vivid. the
drawing goes in the banner, not in the variable name.

## the naming law (2026-08-19)

the whole repository follows one scheme. a name that breaks a rule
below is a defect, not a style choice.

1. module = `<prime>_<subject>`, prime one of graph, token, relation,
   field, operation, transform, map, view — plus `util_` for the
   non-mathematical helpers (string, file, verbosity). file = module,
   one to one. TYPES read english order (`binary_relation`,
   `stored_field`); MODULES read namespace order (`relation_binary`,
   `field_stored`). corollary, statically checkable: no type name
   begins with a prime word plus underscore, so a type can never
   collide with a module name.
2. cardinality is `num_<plural>` — bindings and components alike.
   `size`, `size_of`, `x_size`, `_count`, and `n_` prefixes are out.
3. a predicate is a verb phrase asserting the property it tests
   (`has(x)`, `matches(t)`, `labelled(g)`); set membership is `has`.
4. a reader is the bare noun of what it returns; `get_` is banned.
   a writer is `set_<noun>`.
5. a function is named by its value (a noun); a subroutine by its
   act (an imperative verb, e.g. `derive_faces`).
6. constructors: the generic interface bears the type name; its
   implementations are `create` or `create_<discriminant>`. free
   builders are `x_of(y)` ("the x of y") or `x_from_y`.
7. british spelling with oxford -ize: `colouring`, `neighbours`,
   `centres`, `fibre`, `materialized`.
8. a public constant is `<SUBJECT>_<VALUE>`; the subject names the
   vocabulary it belongs to (`SWEEP_FORWARD`, `WALK_COLOURING`).
9. an abstract interface is `<type>_<binding>_interface`.
10. one type, one name: import-renaming aliases (`x => y`) are out,
    and no module re-exports a name it does not define.
11. an abbreviation is allowed only when the abbreviation IS the
    formal name (`csr`, `gmres`, `bdf`, `lhs`/`rhs`).
12. every module carries `private` plus an explicit public list.

## the writing voice

comments and commit messages are written like a person explaining to a
person: short sentences, concrete pictures, the reason before the
mechanism. read it aloud — if it sounds like a machine wrote it, rewrite it.

say what happens: "the ghosts are pulled fresh from their owners every
product, so the only invariant a vector keeps is a valid owned slab."
not: "ghost coherency is maintained via on-demand synchronization."

banned, from hard experience: invented compound jargon, nicknames coined
mid-project, abbreviation chains, and stacked metaphors. one metaphor,
carried carefully, beats three mixed ones.

## the commit

a commit message is a graph of the change, not an essay. the code itself
is a graph — classes, methods, tables, kernels are the NODES; who calls
whom, who builds what, who asks whom are the EDGES. so draw exactly that:
the before-graph with its flaw visible (a missing edge, a duplicated
node), the after-graph with the change visible (the node added, the edge
wired, the duplicate crossed out ✗). words appear only as labels on nodes
and verbs on arrows.

a real one:

```
   the exchange tables move to the graph

   BEFORE ──────────────────────────────────────────────

      (graph)                    (assembler)
      partition lists,             │ hand loop ──▶ (face table)
      never asked                  │ hand loop ──▶ (node table)
         ▲                         │ hand loop ──▶ (reverse table)
         ╳ no edge                 └─ re-rolled sort ✗
                                      (the kernel exists!)

   AFTER ────────────────────────────────────────────────

      (graph) ──ghost_owners──▶ (face table)
         │  └──ghost_owners──▶ (node table)      three new edges,
         │ ────ghost_copies──▶ (reverse table)   one new node each
         │            │
         │            └──groups by──▶ (counting_sort)
         │                             the one kernel
      (assembler) ──reads tables, sends values──▶ (the wire)

   CHECKED ──────────────────────────────────────────────
      tables pinned entry-by-entry on the chain;
      parallel 1/2/4 unchanged to machine precision
```

the reader should be able to answer, from the picture alone: which node
was added? which missing edge got implemented? which duplicate died?
if the message cannot be drawn this way, the commit is doing too many
things — split it.

rules: lowercase; plain english on every label (a name the textbook does
not define gets its one-clause explanation right on the node); no tool
signatures; CHECKED carries the numbers, never just the word "tested".

## before committing

- the `::` columns line up in every block you touched
- every new type and procedure has a banner, with a drawing if there
  is structure to draw
- no semicolon-stacked lines, no dead code left commented out
- names are formal or plain, never coined
- read the comments aloud once - would a person say that?
- the commit message shows the change as a picture and ends with what
  was checked
