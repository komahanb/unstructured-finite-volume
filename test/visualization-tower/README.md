# The Visualization Tower

An orbital application client of the relational nucleus, and the
seventh in the series. Its purpose is not to print a graph.

Its architectural question is:

> **Can mathematical structure describe itself visually, independently
> of numerical execution?**

The persistent specimen is a chain of three mathematical operators,

```
    X0 --A1--> X1 --A2--> X2 --A3--> X3
```

which through all of Gate A have **no coefficients at all**. What
exists instead is each operator's *structural dependency*,

```
    D1 : X0 -> X1        D2 : X1 -> X2        D3 : X2 -> X3
```

derived from primitive incidence, composed along the chain,
transposed, and rendered.

**Level 5 is where the numbers arrive** — `w_k : E_k -> R`, one
coefficient per dependency occurrence — and its finding is that
nothing structural moves when they do:

> **structure != value**, and `w(e) = 0` does not make `e` or its
> dependency disappear.

---

## Status

**Levels 0–6 are built and green.** Review Gate A is behind them.

**Review Gate B is NOT reached.** Gate B comes after Levels 5, 6 and
7; Level 7 is UNBUILT, and so are 8 and 9. The tower is not sealed and
does not claim to be.

Gate A completed at `44ae3da`/`29c0ccd`; Level 5 at `b134a1f`; Level 6
on top of that.

---

## The ten-level Rosetta table

| Level | Mathematics | Framework object | Input domain | Output domain | Forward structure | Reverse structure | Representation | Exact test | Status |
|---|---|---|---|---|---|---|---|---|---|
| **L0** | seven declared carriers | `counted_set` | — | — | `X0 X1 X2 X3` | — | none | \|X0\|=4, \|X1\|=\|X2\|=3, \|X3\|=2, \|E1\|=5, \|E2\|=4, \|E3\|=3; all 42 ordered pairs distinct | PASS |
| **L1** | typed occurrence incidence | `csr_relation` | `E_k` | `X_(k-1)`, `X_k` | `T_k <= E_k x X_(k-1)`, `H_k <= E_k x X_k` | — | none | twelve occurrences, each with exactly one tail and one head | PASS |
| **L2** | relation algebra | `compose_binary`, `transpose_of` | `X0` | `X3` | `D_k = H_k o T_k^T`; `D2:1`; `D3:1` | `D_k^T`; `D1^T o D2^T o D3^T` | none | exact extensions; 7 walks → 6 tuples; `(D3 o D2 o D1)^T = D1^T o D2^T o D3^T` | PASS |
| **L3** | one relational ownership environment | `relational_graph` | — | — | owns 7 carriers + 6 primitive relations | derived on demand | none | signature closure over all twelve slots; whole chain re-derived from graph-owned relations alone | PASS |
| **L4** | structural interpretation | test-local `structural_renderer_fixture` | any binary relation | text | chain line + sparsity | chain line + sparsity | A, B, C, D, E | five generated representations, cell by cell against `relation % has` | PASS |
| — | ===== **REVIEW GATE A** ===== | | | | | | | | |
| **L5** | fields `w_k : E_k -> R` | `class_graph_field` / `field` | occurrence carriers `E1, E2, E3` | scalar values | structure **unchanged** | structure unchanged; **no numerical reverse** | structural sparsity **+ coefficient view** | `level-5-field-calculus/test.f90` | PASS |
| **L6** | production dependency projection | `discretization_operator % dependencies()` → `class(graph)` | one ordinary vertex carrier | the same carrier | Boolean coordinate pattern equals `D2 : X1 -> X2` | structure unchanged; no numerical reverse | signature **+** sparsity, shown side by side | `level-6-discretization/test.f90` | PASS |
| **L7** | minimization | — | — | — | — | — | — | — | UNBUILT |
| **L8** | constitution | — | — | — | — | — | — | — | UNBUILT |
| **L9** | statement | — | — | — | — | — | — | — | UNBUILT |

No level is marked N/A, and every unbuilt level is marked UNBUILT.
Levels are implementation units; gates are horizontal separators where
review happens, and are never directories. **Level 5 being green does
not make Gate B reached** — Gate B reviews Levels 5, 6 and 7 together,
and two of those do not exist.

---

## The specimen

Seven structurally distinct carriers. Raw integer labels overlap on
purpose — `a`, `p`, `u`, `m`, `e11`, `e21` and `e31` are all the
integer `1`, and all seven are different objects.

```
    X0 = { a  b  c  d }              E1 = { e11 e12 e13 e14 e15 }
    X1 = { p  q  r }                 E2 = { e21 e22 e23 e24 }
    X2 = { u  v  w }                 E3 = { e31 e32 e33 }
    X3 = { m  n }
```

Twelve dependency occurrences, each with one tail and one head:

```
    A1                A2                A3
    -------------     -------------     -------------
    e11 : a -> p      e21 : p -> u      e31 : u -> m
    e12 : b -> p      e22 : q -> u      e32 : v -> n
    e13 : b -> q      e23 : q -> v      e33 : w -> n
    e14 : c -> q      e24 : r -> w
    e15 : d -> r
```

---

## Five different things, and none of them is another

This is the section to read before any of the pictures below.

### 1. The operator `A_i`

A mathematical map `X_(i-1) -> X_i`. **At Gate A it has not been
instantiated numerically.** It has no entries, it is never evaluated,
and no vector is ever pushed through it. It is the thing the tower is
*about*, and it is the only one of these five that this tower does not
yet possess.

### 2. The dependency relation `D_i`

A typed binary relation `D_i <= X_(i-1) x X_i` saying *which* members
could influence which. It exists, it is derived, and it is the
tower's whole content at Gate A.

### 3. The numerical coefficient

A number attached to a dependency occurrence — `w_k : E_k -> R`.
**Absent through all of Gate A; introduced at Level 5.** It lives on
`E_k` and on nothing else, and attaching it changes no relation, no
tuple and no structural picture.

### 4. The sparsity representation

A rectangular arrangement of filled and empty cells with rows and
columns in a chosen order. Derived from `D_i`, and not the same
object: `D_i` has no rows.

### 5. The rendered picture

Characters on a page. Derived from the representation. It is
downstream of everything above it and feeds back into none of it.

**Therefore the tower's central Gate-A hypothesis, now tested from
both sides:**

> The structural picture can precede the numerical coefficients
> entirely — and once they arrive, it does not change.

Levels 0–4 exercise items 2, 4 and 5 while items 1 and 3 do not exist.
Level 5 adds item 3 and shows items 2 and 4 unmoved.

That is **enforced rather than asserted**. Refusing
`class_graph_field` and `graph_field_calculus` would not have been
enough — a level could simply have written `real(dp) :: w = 2.0_dp`
and helped itself. So `check_imports.sh` carries a *numberless law*:
no `real` or `complex` declaration, and no literal carrying a decimal
point or an exponent, in any source **below Level 5**. Comments are
stripped first (the prose says "real" and means the English word);
integers are untouched, because a carrier's size is 4 and that is not
a coefficient.

**The law is a ceiling that lifts at exactly one rung, not a rule that
was deleted when it became inconvenient.** Level 5 and its two
coefficient fixtures are entitled to numbers; Levels 0–4 and the
structural renderer are not, and still are not. The selftest checks
both directions, and the gate refuses a planted
`real :: sneaky = 2.0` in Level 4 while accepting the real fields in
Level 5.

---

## Forward and reverse, side by side

```
    FORWARD

    X0 --D1--> X1 --D2--> X2 --D3--> X3


    STRUCTURAL REVERSE

    X3 --D3^T--> X2 --D2^T--> X1 --D1^T--> X0
```

and the law the tower exists to state:

```
    D3:1  =  D3 o D2 o D1

    D3:1^T  =  D1^T o D2^T o D3^T
```

verified extensionally — domains by `same_as`, cardinality, and
membership both ways round — never by comparing tuple enumeration
order, because a relation is a set.

### This is structural dual orientation. It is not yet adjoint numerical evaluation.

The only adjoint-shaped statement Gate A makes is

```
    supp(A^T)  =  supp(A)^T
```

and even that is a statement about *supports*, made in a tower where
no `A` has been given entries. There is no `A^T v`, no `lambda`, no
dual pairing, no gradient, no sensitivity, no reverse accumulation and
no adjoint solver. Those questions belong to the Adjoint Tower, which
already owns them.

What `D^T` says, and all it says, is: **which members of `X0` could
have reached this member of `X3`.**

---

## Level 5 — what changes when values arrive

Nothing structural. That is the finding.

```
    w1 : E1 -> R  =  [ 2, -1,  0,  3,  4 ]
    w2 : E2 -> R  =  [ 1,  5, -2,  2 ]
    w3 : E3 -> R  =  [ 3, -1,  4 ]
```

### The coefficients live on `E`, and on nothing else

`e12` and `e13` both read `b`; they carry `-1` and `0`. A field on
`X0` could not hold both, and a field on `X1` could not either. `E_k`
is where the two ends meet, and it is the only one of the three
carriers that can seat the number.

This is why Level 0 gave the occurrences first-class identity five
levels before anything needed a value: **a coefficient needs somewhere
to live that survives being zero.**

The domain is checked by **identity**, never by size. The specimen
supplies the hostile case for free — `|X0| = |E2| = 4` — so a
four-valued field on `X0` exists, is perfectly valid, and is *refused*
as `E2`'s coefficients.

### `0` is not `.`

Two representations of the same relation:

```
    D1 STRUCTURE                 D1 VALUES

            a b c d                         a    b    c    d
    p       # # . .              p          2   -1    .    .
    q       . # # .              q          .    0    3    .
    r       . . . #              r          .    .    .    4
```

```
    #  = structural dependency
    0  = structural dependency carrying the value zero
    .  = no structural dependency
```

and therefore, within this representation:

> **`0` != `.`**

`w1(e13) = 0` where `e13 : b -> q`, and `b -> q` is still in `D1`. The
cell at row `q`, column `b` prints `0`. The cell at row `q`, column
`a` prints `.` — because there is no dependency there, so there is
nothing to carry a value.

A representation that wrote `0` for both would have thrown the
distinction away.

### Two independent queries, joined only for the picture

The pipeline forks and rejoins exactly once:

```
    primitive incidence                occurrence carrier E
            |                                  |
            v                                  v
    derived dependency structure          field w : E -> R
            |                                  |
            v                                  v
    structural picture                  coefficient values
            |                                  |
            +---------------+------------------+
                            |
                            v
                   valued visualization
```

**Never** redrawn as `coefficient field -> dependency structure`. That
would invert the architecture. The forbidden inference is

```fortran
if (coefficient /= 0) then draw('#')
```

and it appears nowhere. A cell is present because
`relation % has(x, y)` says so, and for no other reason; the value is
consulted only to decide *how the number looks*, never whether the
cell exists.

### The independence probe

A second field on the **same** `E1` — `w1_alt = [9, 8, 7, 6, 5]`,
sharing no value with the first and containing no zero:

```
    D1 VALUES (w1)               D1 VALUES (w1_alt)

               a    b    c    d            a   b   c   d
    p          2   -1    .    .   p        9   8   .   .
    q          .    0    3    .   q        .   7   6   .
    r          .    .    .    4   r        .   .   .   5
```

Identical structural picture; different coefficient picture; the dots
in exactly the same places. So

> **`D1` does not determine `w1`, and `w1` does not determine `D1`.**

### What Level 5 deliberately does not do

- **No composed coefficients.** `D2:1` and `D3:1` stay structural.
  Giving them values would mean choosing an algebra for numerical
  composition — sums over intermediate members — which is operator
  mathematics, not field calculus.
- **No `w^T`.** Transposing a coefficient is a numerical act. Gate A's
  structural reverse asked for none and still does not.
- **No numerical operator is executed.** No `A` is applied to
  anything, and no `A^T v` exists.
- **Nothing in the Level-4 renderer changed.** The valued renderer
  *imports* it; Level 4 remains able to render mathematics that has no
  numbers at all, and the numberless law still forbids it from holding
  one.

---

## Level 6 — what production actually exposes

The first level to name production machinery, and the **first
executable consumer `dependencies()` has ever had.**

### The production census, at `b134a1f`

```
discretization_operator                    src/graph_calculus.f90:219
│   deferred :: dependencies
│   subroutine (this, pattern)
│       class(graph), allocatable, intent(out) :: pattern
│
├── stencil_operator                       src/class_graph_stencil.f90:47
│   ├── implementation   stencil_dependencies              :168
│   ├── meaning          a copy of its own stored pattern, built at
│   │                    construction as
│   │                    stored_graph(nv, tails=columns, heads=rows),
│   │                    nv = size(constant)
│   ├── returns          class(graph) → concretely stored_graph
│   └── repo callers     0
│
└── step_operator                          src/class_graph_step.f90:46
    ├── implementation   step_dependencies                 :127
    ├── meaning          a freshly built linear chain of reach+1
    │                    instants, tails=[1..reach], heads=[2..reach+1]
    ├── returns          class(graph) → concretely stored_graph
    └── repo callers     0

repo-wide executable callers of % dependencies():  ZERO
```

Three separate facts, and they must not be run together:

| | |
|---|---|
| the family contract exists | `graph_calculus.f90:219`, with a signature at `:406–410` |
| two implementations exist | `stencil_dependencies`, `step_dependencies` |
| **no consumer exercises it** | zero call sites, anywhere in the repository |

The contract's own prose says *"the minimizers one level up
interrogate the pattern — the diagonal, the colouring, the
triangularity, the Galerkin road — so it is exposed by law, never by
inspection."* No minimizer does. That sentence is recorded here as

```
    PROSE INTENTION — NO EXECUTABLE CALLER FOUND
```

and is **not** repeated as an established repository fact. It is a
**latent contract**, not dead code: two faithful implementations
waiting for a consumer. Level 6 is the first.

### Measurement one — the stencil witness

A production `stencil_operator` is built carrying D2's Boolean
occupancy in production's own coordinates, with Level 5's `w2` as its
weights. Then `dependencies()` is called.

```
    RELATIONAL D2                    STENCIL dependencies()
    signature: X1 -> X2              signature: vertices -> vertices

            p q r                            1 2 3
    u       # # .                    1       # # .
    v       . # .                    2       . # .
    w       . . #                    3       . . #

    coordinate pattern equal ........................ YES
    typed source identity equal ..................... NO
    typed target identity equal ..................... NO
```

`D2` stands on **two** declared carriers and production's answer on
**one**, which is neither of them. All three sets hold three members.

> **STENCIL-B** — coordinate-equivalent to `D2`, and it loses the
> distinct typed source and target carriers.

### Measurement two — the rectangular witness

`D1 : X0 -> X1` runs `4 -> 3`. The production constructor takes a
single vertex count, so `|V|` would have to be 4 and 3 at once. Given
4 it holds every one of D1's five arrows — and then reports a row that
D1's codomain does not have:

```
    RECTANGULAR D1                   STENCIL dependencies()
    signature: X0 -> X1              signature: vertices -> vertices
    shape: 3 rows x 4 columns

            a b c d                          1 2 3 4
    p       # # . .                  1       # # . .
    q       . # # .                  2       . # # .
    r       . . . #                  3       . . . #
                                     4       . . . .    <-- phantom

    production contract preserves this typed signature . NO
```

The arrows survived. The **signature** did not. No union carrier was
manufactured, no domain was padded, and nothing was indexed out of
range — the fourth row is simply what one carrier serving both axes
looks like.

> **RECT-B** — the ordinary-graph answer cannot intrinsically
> represent distinct 4-member source and 3-member target carriers.

### Measurement three — the step witness

The same family verb, asked of the other concrete citizen:

```
    BDF2 dependencies()
    signature: vertices -> vertices

            1 2 3
    1       . . .
    2       # . .
    3       . # .

    wrapped state pattern equal ..................... NO
    same BDF2 motif under second wrapped sparsity ... YES
```

Three vertices, like the stencil's answer — and a different edge
structure entirely. The stencil says *state 1 depends on state 1*; the
step says *no instant is its own predecessor*. Wrapping a second
action whose state sparsity is genuinely different (a diagonal) yields
the **identical** motif, which is executable evidence that the step's
`dependencies()` describes the scheme's **time axis** and not the
wrapped action's algebra.

> **FAMILY-B** — the two concretes use one procedure for two different
> structural axes: `stencil` for state sparsity, `step` for temporal
> reach.

**This is not a defect in either.** Each faithfully implements the
narrower mathematics it actually represents, and no executable
contract, test or caller promises otherwise. A defect would require a
promise that behaviour violates; the census found no promise with a
consumer behind it.

### The visual equality theorem

```
    V(R1) = V(R2)   does NOT imply   R1 = R2
```

when the representation `V` forgets carrier identity. `D2 : X1 -> X2`
and `pattern : vertices -> vertices` render an identical grid and are
different objects.

**Therefore a richer visualization preserves the signature.** Every
Level-6 picture prints `signature:` above its grid, because the grid
alone is exactly the part that cannot tell the two apart. This is the
tower demonstrating, on itself, which information a representation
retains and which it discards.

`(signature, sparsity)` is **not** turned into a production type here.

### Nothing was applied

`stencil % apply` and `step % apply` are called **zero** times, and
that is mechanical rather than promised: `check_imports.sh` refuses a
`% apply(` in any tower source. The weights, constants and step size
exist only because the production constructors require complete
objects.

No numerical composition of `A3 A2 A1`, no transpose stencil, no
`A^T v`. `D2:1` and `D3:1` remain relation-algebra results.

---

## Mathematics and Fortran spelling

The repository's composition reads its arguments in the order the data
flows, and its result in the order mathematics writes it:

```
    compose_binary(R_AB, R_BC)  =  R_BC o R_AB
```

so the **first argument is applied first**. Written out:

| Mathematics | Fortran |
|---|---|
| `D_k = H_k o T_k^T` | `compose_binary(T_k^T, H_k)` |
| `D2:1 = D2 o D1` | `compose_binary(D1, D2)` |
| `D3:1 = D3 o D2:1` | `compose_binary(D2_1, D3)` |
| `D_rev = D1^T o D2^T o D3^T` | `compose_binary(compose_binary(D3^T, D2^T), D1^T)` |

The two directions of reading are genuinely opposite. That is resolved
in exactly one place — `common/visualization_algebra_fixture.f90` —
and every level above it says `derive_composition(first, second)` and
means the arrow.

---

## What Level 4 generates

Every line below was produced from the twelve occurrences at run time.
None of it is stored, and the renderer has never seen the expected
output.

```
    FORWARD CHAIN
    X0 --D1--> X1 --D2--> X2 --D3--> X3

    D1                 D2                 D3
            a b c d            p q r              u v w
    p       # # . .    u       # # .      m       # . .
    q       . # # .    v       . # .      n       . # #
    r       . . . #    w       . . #

    D2:1               D3:1
            a b c d            a b c d
    u       # # # .    m       # # # .
    v       . # # .    n       . # # #
    w       . . . #

    STRUCTURAL REVERSE
    X3 --D3^T--> X2 --D2^T--> X1 --D1^T--> X0

    D3:1^T             D1^T o D2^T o D3^T
            m n                m n
    a       # .        a       # .
    b       # #        b       # #
    c       # #        c       # #
    d       . #        d       . #
```

The last two grids are drawn from two independently composed
relations, and are required to agree cell for cell while differing in
their name line. That is Level 2's transpose/composition law seen
through a representation instead of through `same_extension`.

The borderless grid style follows the brief's §30 exactly. The
bordered style shown in §5 and §16 (`+------` / `|`) is the same
content with rules drawn in; both sections state that the artistic
formatting is not prescribed and that the semantic content is the
oracle.

### The witness collapse, visible in the picture

Seven walks run `X0 -> X1 -> X2`; `D2:1` holds six tuples. Both
`b->p->u` and `b->q->u` witness `b->u`, and the cell at row `u`,
column `b` is filled **once**. A relation is a set: two witnesses of
one fact are one fact.

### Order comes from declaration, always

Every axis is walked `member(1), member(2), ..., member(n)`. Never
sorted, never assumed to be `1..n`, never read off the tuple table —
which is a set and has no order to read.

Level 4 checks this with a carrier built to be hostile: a
`subset_set` declaring its members `{ 30, 10, 20 }` in that order
renders

```
    P
            30 10 20
    m       #  .  .
    n       .  .  #
```

and not `10 20 30`. No other carrier in the specimen could catch a
renderer that sorted, because every other one counts from one and
would look correct either way.

---

## Rendering is an interpretation, not a source of truth

The required flow, and the only one the tower permits:

```
    primitive incidence
            |
            v
    relation algebra
            |
            v
    dependency relation
            |
            v
    structural representation
```

The forbidden flow — write the expected picture, then invent the
relation that reproduces it — is closed off structurally: the oracle
pictures live in `level-4-graph-calculus/test.f90`, and not one line
of them is reachable from `common/structural_renderer_fixture.f90`.

The renderer knows nothing about the specimen. It has never heard of
`a`, `p`, `u` or `m`; it does not know that `D1` has four columns, or
that its first column is full. Every one of those facts is obtained,
at the moment of drawing, by asking `domain(1)`, `domain(2)`,
`size()`, `member(i)`, `has([col, row])`, `name()` and
`label_for(carrier, member)` — and nothing else.

### Interpreting needed strictly less of the nucleus than deriving did

Worth stating separately, because it makes the Gate-A answer stronger
than "the nucleus was sufficient":

| Job | What it needed |
|---|---|
| **deriving** the structure | `binary_relation` (`source`, `target`, `transpose_of`), `compose_binary`, `relational_graph` |
| **interpreting** it | the **root `relation`** contract (`arity`, `domain(k)`, `has`, `name`) and `member_set` (`size`, `member`, `local_index`, `name`) |

`structural_renderer_fixture` names no binary relation at all. It never
calls `source()`, `target()`, `image_view()` or `transpose_of()` — it
is handed transposed and composed relations and cannot tell them from
any other. Rendering a rectangular typed dependency turned out to need
*less* of the nucleus than building one.

---

## What was *not* built

- **No production visualization abstraction.** No
  `graph_visualization`, `class_graph_visualization`,
  `graph_observation`, `graph_interpretation`, `graph_renderer` or
  `graph_visualizer`; no `visualize()`, `print()` or `dependencies()`
  added to `graph_operation`.
- **No union carrier.** There is no
  `V = X0 u X1 u X2 u X3`. Level 3 asserts the graph does not hold
  one, and Level 4 shows what adopting one would cost.
- **No production change at all.** `git diff 44ae3da -- src/`,
  `git diff 29c0ccd -- src/` and `git diff b134a1f -- src/` are all
  empty. Level 5 needed no correction — the field machinery already
  seats a coefficient on an occurrence carrier and already answers
  domain questions by identity. Level 6 needed none either: it is a
  measurement, and what it measured was a specialization, not a
  defect.
- **`dependencies()` was not generalized**, its return type was not
  changed, and nothing was moved onto `graph_operation`.

The tower succeeded, and success here means: **the nucleus was already
sufficient.** That is a reason to build nothing, not a reason to build
a hierarchy.

---

## The Gate-A architectural experiment

**Question one.** Can the rectangular typed relation `D1 : X0 -> X1`
be rendered directly from the relational nucleus?

**Yes**, and the whole of Level 4 is the demonstration. Recorded:
*ordinary graph not required for structural visualization at this
radius.*

**Question two.** Would forcing `D1` into `ordinary_graph_view`
collapse `X0` and `X1`, or otherwise change its mathematical meaning?

The two profile contracts were inspected rather than provoked:

| Profile | Schema it demands | Why the specimen cannot satisfy it |
|---|---|---|
| `ordinary_graph_view` | `T <= E x V` and `H <= E x V`, one `V` | `T1` lands in `X0`, `H1` in `X1`; refuses with *"the head relation must share the tail's domains"* |
| `directed_adjacency_view` | `A <= V x V`, one `V` | `D1`'s source is `X0`, its target `X1`; refuses with *"a directed adjacency runs over one domain"* |

Both readings demand a **single** vertex carrier. Satisfying either
would require manufacturing `V = X0 u X1 u ...`, and it is worth being
exact about what that costs and what it does not.

**Direction survives.** An ordinary directed graph preserves
ordered-pair direction perfectly well: `directed_adjacency_view`
documents that "the tuple order of a same-domain binary relation IS
the direction", and `(a, p)` in a union carrier would still be
distinguishable from `(p, a)`. Nothing here claims otherwise.

**The typed signature does not.** What a union loses is the *intrinsic
distinct source and codomain carriers* of the rectangular relation:

```
    D : X_i -> X_j        against        D^T : X_j -> X_i
```

Under a union both become relations over one domain `V`, and the two
declared carriers that made them different KINDS of object are gone —
recoverable only from an offset convention the mathematics never
stated. The first thing `same_extension` tests is domain identity, and
after a collapse that test has nothing left to compare.

So the honest statement is:

> An ordinary graph can preserve ordered-pair direction, but it does
> not intrinsically preserve the distinct typed source/codomain
> carriers of the rectangular relation.

This is not a defect in the profile. It is a specialization the
profile documents, and it is the right specialization for the mesh
path it was written for.

**Gate-A answer: C.** The relational nucleus is sufficient; the
ordinary-graph profile is intentionally too specialized for
rectangular typed dependency.

---

## Import stratification

`check_imports.sh`, fail-closed, keyed per level and per shared
fixture, with `--selftest`.

| Level | Ceiling |
|---|---|
| L0 | `graph_carrier` |
| L1 | `+ graph_relation`, `graph_binary_relation` |
| L2 | `+ graph_relation_algebra` |
| L3 | `+ graph_structure` |
| L4 | `+ structural_renderer_fixture` |
| — | ===== **REVIEW GATE A** ===== |
| L5 | `+ class_graph_field`, `graph_field_calculus`, `visualization_values_fixture`, `valued_renderer_fixture` |
| L6 | `+ graph_grammar`, `class_graph_stencil`, `class_graph_step`, `production_discretization_fixture`, `production_pattern_renderer_fixture` |

The shared fixtures are keyed by file and classified by the first
level that earns them: `visualization_assert` below everything,
`visualization_carriers_fixture` at L0,
`visualization_relations_fixture` at L1,
`visualization_algebra_fixture` at L2,
`structural_renderer_fixture` **at L4** — so no level below Level 4
can quietly become a picture — `visualization_values_fixture` plus
`valued_renderer_fixture` **at L5**, and
`production_discretization_fixture` plus
`production_pattern_renderer_fixture` **at L6**.

Two families are refused **universally**, at every level built so far.
Values used to be a third, and Level 5 is where it stopped being one:

0. **Values** — `class_graph_field`, `graph_field_calculus`, *and any
   real number written by hand*, refused at L0–L4 and **earned at
   L5**. See the numberless law above; the module refusal alone would
   have left the back door open, and the ceiling lifting at exactly
   one rung is what keeps Level 4's claim mechanical after Gate A.
2. **Operators, refused at L0–L5 and PARTLY earned at L6.**
   `graph_grammar`, `class_graph_stencil` and `class_graph_step` lift
   at Level 6 — and only those three, because only those three are
   used. `graph_calculus`, `graph_fitting`,
   `class_graph_linearization`, `graph_minimization`,
   `class_graph_gmres`, `class_graph_newton`, `class_graph_marcher`
   and the derivative/adjoint fixtures stay refused **everywhere**,
   Level 6 included.

   The type names `discretization_operator`, `stencil_operator` and
   `step_operator` are additionally refused by a direct source scan,
   in case a later level finds another road to them — and that scan is
   level-sensitive in the same way: **forbidden L0–L5, allowed L6.**
   A planted `use class_graph_stencil` in a Level-5 source is refused
   twice over, by the module scan and by the name scan.

2b. **Application, refused everywhere.** A `% apply(` in any tower
   source is refused at every level. This tower interrogates
   structure; the moment a level evaluates a production operator it
   has started doing arithmetic, and Level 6's "zero applies" would
   become a promise instead of a fact.
3. **The ordinary graph** — `graph_profile`, `graph_algorithms`,
   `class_graph`, `class_stored_graph`. Level 5 included: values
   arrived, machinery did not.

**That third refusal is load-bearing evidence, not hygiene.** Gate A
concludes that a rectangular typed dependency needs no ordinary-graph
reading. Because no level of this tower may name `graph_profile`, the
pictures Level 4 produced cannot have leaned on one. The gate is what
makes "no ordinary graph was required" checkable instead of promised.

---

## The frontier — documented, not implemented

Level 5 asked *"does a numerical zero erase structural dependency?"*
and answered **no**. Level 6 asked whether a production
discretization operator exposes the same skeleton, and answered
**STENCIL-B / RECT-B / FAMILY-B**; see the Level-6 section above.
What follows is still unbuilt.

### Level 6 — done, and what it left open

Where production first enters. The expected question:

> How does the structural relation `D` correspond to
> `discretization_operator % dependencies()`?

The answer was the second of the anticipated findings:
`dependencies()` is ordinary-graph/square-domain specialized and
cannot intrinsically express `X -> Y`. No RED occurred, because no
executable contract promised otherwise.

**What Level 6 did NOT answer**, and deliberately: whether
`dependencies()` belongs on `graph_operation` at all. This tower
inspected `discretization_operator` and its two concretes and nothing
else. The root question —

> Does every `graph_operation` possess one meaningful dependency
> structure?

— is untouched, and the ledger records root pressure as **OBSERVE**,
never REFACTOR.

**A deeper possibility the two witnesses raise.** The stencil exposes
a state axis and the step a temporal one. The pressure this tower
found may therefore not be *"we need one universal
`dependencies()`"* but *"an operation may admit several legitimate
structural interpretations"* — `D_state`, `D_time`, `D_space`,
`D_block`. The evidence is recorded in VIZ-20 and VIZ-24. **The
abstraction is not designed here.**

### Level 7 — the question Level 6 sharpened

Visualization does not obviously require minimization, and Level 7
was always a candidate for N/A. Level 6 gave the question a sharper
form: the census found **zero** existing callers of `dependencies()`,
and the contract's prose names *minimizers* as the intended consumer.
So Level 7 should first test whether minimization actually consumes
discretization structure at all.

The level stays UNBUILT rather than marked N/A. It must be inhabited
or explicitly refused **after** evidence, never to fill a row.

### Level 8 — constitution

> Can an actual operator chain and its structural representation
> coexist compositionally?

```
    A3 o A2 o A1
         |
         +-- numerical execution
         |
         +-- structural interpretation
                  |
                  v
                picture
```

A shape like `operation (x) representation` may emerge. It is not
implemented, not named, and not designed here — Gate B must earn it.

### Level 9 — statement

Potentially `visualize_structure(A3 o A2 o A1)`, returning forward
chain, reverse chain, individual sparsity and composed sparsity while
leaving the mathematical operator unchanged. The exact interface is
deliberately unresolved.

---

## Layout

```
test/visualization-tower/
├── README.md                    this file
├── NUCLEUS-OBSERVATIONS.md      the observation ledger
├── run.sh                       the ladder, stopping at the frontier
├── check_imports.sh             fail-closed, per level, --selftest
├── common/
│   ├── visualization_assert.f90            below everything
│   ├── visualization_carriers_fixture.f90  earned at L0
│   ├── visualization_relations_fixture.f90 earned at L1
│   ├── visualization_algebra_fixture.f90   earned at L2
│   ├── structural_renderer_fixture.f90     earned at L4
│   ├── visualization_values_fixture.f90    earned at L5
│   ├── valued_renderer_fixture.f90         earned at L5
│   ├── production_discretization_fixture.f90    earned at L6
│   └── production_pattern_renderer_fixture.f90  earned at L6
├── level-0-carrier/
├── level-1-relation/
├── level-2-relation-algebra/
├── level-3-graph/
├── level-4-graph-calculus/
├── level-5-field-calculus/
└── level-6-discretization/
```

There is no `gate-a/` directory and no `gate-b/` directory, and there
will not be either. Gates are where you stop and report; levels are
where you build.
