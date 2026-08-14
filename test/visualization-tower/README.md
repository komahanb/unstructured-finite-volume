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

which at Gate A have **no coefficients at all**. Not one real number
appears anywhere in this tower. What exists instead is each
operator's *structural dependency*,

```
    D1 : X0 -> X1        D2 : X1 -> X2        D3 : X2 -> X3
```

derived from primitive incidence, composed along the chain,
transposed, and rendered.

---

## Status

**REVIEW GATE A — reached.** Levels 0–4 are built and green. Levels
5–9 are UNBUILT. The tower is not sealed and does not claim to be.

Starting HEAD: `44ae3da`.

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
| **L5** | values on dependency occurrences | — | — | — | — | — | — | — | UNBUILT |
| **L6** | production operator structure | — | — | — | — | — | — | — | UNBUILT |
| **L7** | minimization | — | — | — | — | — | — | — | UNBUILT |
| **L8** | constitution | — | — | — | — | — | — | — | UNBUILT |
| **L9** | statement | — | — | — | — | — | — | — | UNBUILT |

No level below Gate A is marked N/A, and no level above it is marked
anything but UNBUILT. Levels are implementation units; gates are
horizontal separators where review happens, and are never
directories.

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
**Absent at Gate A**, and the subject of the future Level 5.

### 4. The sparsity representation

A rectangular arrangement of filled and empty cells with rows and
columns in a chosen order. Derived from `D_i`, and not the same
object: `D_i` has no rows.

### 5. The rendered picture

Characters on a page. Derived from the representation. It is
downstream of everything above it and feeds back into none of it.

**Therefore the tower's central Gate-A hypothesis:**

> The structural picture can precede the numerical coefficients
> entirely.

Levels 0–4 exercise items 2, 4 and 5 while items 1 and 3 do not exist.
That is the claim, and it is executable: `visualization_assert` holds
no real constant of any kind, and the import gate refuses
`class_graph_field` and `graph_field_calculus` at every level.

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
- **No production change at all.** `git diff 44ae3da -- src/` is
  empty.

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
would require manufacturing `V = X0 u X1 u ...` — and then `D3:1` and
`D3:1^T` would both be relations over one domain, so the orientation
Level 2 proved would stop being expressible at all. The first thing
`same_extension` tests is domain identity.

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

The shared fixtures are keyed by file and classified by the first
level that earns them: `visualization_assert` below everything,
`visualization_carriers_fixture` at L0,
`visualization_relations_fixture` at L1,
`visualization_algebra_fixture` at L2, and
`structural_renderer_fixture` **at L4** — so no level below Level 4
can quietly become a picture.

Three families are refused **universally**, at every level:

1. **Values** — `class_graph_field`, `graph_field_calculus`. Gate A
   holds no number, and the gate is what makes that mechanical.
2. **Operators** — `graph_grammar`, `class_graph_stencil`,
   `class_graph_step`, `graph_calculus`, `graph_fitting`,
   `class_graph_linearization`, `graph_minimization`,
   `class_graph_gmres`, `class_graph_newton`, `class_graph_marcher`,
   and the derivative/adjoint fixtures. The type names the brief
   forbids by name — `discretization_operator` and `stencil_operator`
   — are additionally refused by a direct scan, in case a later level
   finds another road to them.
3. **The ordinary graph** — `graph_profile`, `graph_algorithms`,
   `class_graph`, `class_stored_graph`.

**That third refusal is load-bearing evidence, not hygiene.** Gate A
concludes that a rectangular typed dependency needs no ordinary-graph
reading. Because no level of this tower may name `graph_profile`, the
pictures Level 4 produced cannot have leaned on one. The gate is what
makes "no ordinary graph was required" checkable instead of promised.

---

## The frontier — documented, not implemented

### Level 5 — values on dependency occurrences

The expected experiment attaches `w_k : E_k -> R` with at least one
zero coefficient, e.g.

```
    w1 = [ 2, -1,  0,  3,  4 ]
    w2 = [ 1,  5, -2,  2 ]
    w3 = [ 3, -1,  4 ]
```

and asks:

> **Does a numerical zero erase structural dependency?**

Expected answer: **no**. A structurally existing occurrence with
coefficient zero remains structurally present. This tests
`structure != value`, and it is why `E_k` was given first-class
identity at Level 0 — a coefficient needs somewhere to live that
survives being zero.

### Level 6 — production operator structure

Where production first enters. The expected question:

> How does the structural relation `D` correspond to
> `discretization_operator % dependencies()`?

`discretization_operator`, `stencil_operator` and `step_operator` are
to be audited **only at Level 6**, and Level 6 must RED against the
persistent typed chain before changing anything. Possible findings
include *"`dependencies()` already expresses the correct skeleton"*,
*"`dependencies()` is ordinary-graph/square-domain specialized and
cannot faithfully express `X -> Y`"*, or something else. Nothing is
pre-decided, and the existing `dependencies()` is **not** generalized
at Gate A.

### Level 7 — likely N/A, not yet refused

Visualization does not obviously require minimization. But a later
layout algorithm might genuinely introduce optimization, so the level
is left UNBUILT rather than marked N/A. It must be inhabited or
explicitly refused **after** evidence, never to fill a row.

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
├── run.sh                       the ladder, stopping at Gate A
├── check_imports.sh             fail-closed, per level, --selftest
├── common/
│   ├── visualization_assert.f90            below everything
│   ├── visualization_carriers_fixture.f90  earned at L0
│   ├── visualization_relations_fixture.f90 earned at L1
│   ├── visualization_algebra_fixture.f90   earned at L2
│   └── structural_renderer_fixture.f90     earned at L4
├── level-0-carrier/
├── level-1-relation/
├── level-2-relation-algebra/
├── level-3-graph/
└── level-4-graph-calculus/
```

There is no `gate-a/` directory, and there will not be one. Gates are
where you stop and report; levels are where you build.
