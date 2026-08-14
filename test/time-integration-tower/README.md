# Time Integration Tower

The sixth orbital client of the core-math nucleus, after Calculator,
Learning, Derivative Action, Adjoint Sensitivity and Partitioned
Implicit PDE — and the first to ask what **time** is, structurally,
before anything moves through it.

> **This tower is built vertically, and it is not a pretext for a
> refactor.** Levels 0–4 ask one question only: *what is time, before a
> state value, an ODE law, a scheme or a solver exists?*

## Status — REVIEW GATE A

```text
TIME INTEGRATION TOWER

    L0 carrier ......................... PASS
    L1 relation ........................ PASS
    L2 relation algebra ................ PASS
    L3 relational graph ................ PASS
    L4 graph calculus .................. PASS

    ===== REVIEW GATE A =====

    L5 field calculus .................. UNBUILT
    frontier stops here
```

The tower is **not sealed**. Levels 5–9 are unbuilt, and nothing in
this document claims a truth they have not established.

## Levels are the architecture; gates are only review checkpoints

Ten levels, 0 through 9. Each level is a directory, a test, an import
ceiling and a distinct mathematical responsibility. Review happens at
three checkpoints — after Level 4, after Level 7, after Level 9 — and
that is *all* a gate is. Gates own no mathematics and appear in no
directory name.

---

# A. The structural specimen

Gate A carries **no numbers**. Its specimen is three carriers and the
incidence between two of them:

```text
STATE DOMAIN Q                      TIME DOMAIN T

    x     y                             t0 → t1 → t2 → t3 → t4
                                          e1   e2   e3   e4

    two state coordinates           five instants, four steps
```

```text
Q = { x, y }                    state coordinates      |Q| = 2
T = { t0, t1, t2, t3, t4 }      time instants          |T| = 5
E = { e1, e2, e3, e4 }          time steps             |E| = 4
```

**These are not one domain.** A state coordinate is not a time instant;
a time instant is not a step. The tower's first act is to make that
structural rather than notational.

## The eventual specimen — declared, not implemented

The dynamical law this tower will eventually carry:

\[
\frac{dx}{dt} = -x, \qquad \frac{dy}{dt} = x - y,
\qquad\text{i.e.}\qquad
S(q) = \begin{bmatrix} x \\ y - x \end{bmatrix},
\quad \dot q = -S(q),
\]

with \(q_0 = [2,0]^{T}\) and nominal step \(h = 1/2\).

**None of this exists at Gate A.** No value, no step size, no residual,
no scheme. It is written here so the reader knows the specimen is one
specimen carried through, and so that Levels 5+ inherit it rather than
invent it.

---

# B. The ten-level Rosetta table

| Level | Mathematical object | Framework object | Domains | Relations | Interpretation | Exact test | Truth established | Production consequence |
|---|---|---|---|---|---|---|---|---|
| **0** | \(Q, T, E\) | `counted_set` | three carriers | — | none — sets only | `level-0-carrier/` | three independent identities; member 1 is \(x\), \(t_0\) **and** \(e_1\), and only identity separates them | none |
| **1** | \(\mathrm{Tail},\mathrm{Head}\subseteq E\times T\) | `csr_relation` | \(E\times T\) | primitive temporal incidence | time acquires *direction* | `level-1-relation/` | direction is relation structure; \(Q\) participates in **no** relation, and that absence is the point | none |
| **2** | \(A_1=\mathrm{Head}\circ\mathrm{Tail}^{T}\); \(A_2=A_1\circ A_1\) | `transpose_of`, `compose_binary` | \(T\to T\) | **derived** | one-step and two-step temporal *reach* | `level-2-relation-algebra/` | reach is generated, not stored; \(A_1\neq A_2\) | none |
| **3** | \(G_{\text{time}}=(\{Q,T,E\},\{\mathrm{Tail},\mathrm{Head},A_1,A_2\})\) | `relational_graph`, `held_set`, `held_relation` | three owned carriers | four owned relations | one relational structure, not one domain | `level-3-graph/` | signature closure holds, **and** \(Q\) is lawfully owned while naming no relation | none |
| **4** | directed views over \(A_1\), \(A_2\) | `directed_adjacency_view`, `graph_profile`, `graph_algorithms` | \(T\) (both views) | borrowed from the graph | **causality** | `level-4-graph-calculus/` | forward causal order \([t_0..t_4]\); two views, one carrier | none |
| **5** | — | — | — | — | — | — | **UNBUILT** | — |
| **6** | — | — | — | — | — | — | **UNBUILT** | — |
| **7** | — | — | — | — | — | — | **UNBUILT** | — |
| **8** | — | — | — | — | — | — | **UNBUILT** | — |
| **9** | — | — | — | — | — | — | **UNBUILT** | — |

The road as far as it has been built:

```text
L0  three carriers, three identities
L1  primitive temporal incidence — time gets a direction
L2  derived one-step and two-step reach
L3  one relational structure owning both axes without conflating them
L4  causal traversal: sources, sinks, reachability, topological order
--- REVIEW GATE A ---
L5  field calculus       UNBUILT
L6  discretization       UNBUILT
L7  minimization         UNBUILT
L8  constitution         UNBUILT
L9  statement            UNBUILT
```

---

# C. The persistent object, and the axis that is *not* time

```text
STATE DOMAIN Q                 TIME DOMAIN T
                               
   x      y                      t0 --e1--> t1 --e2--> t2 --e3--> t3 --e4--> t4
                               
   no values yet                 Tail ⊆ E×T    e_i → t_{i-1}
   no field yet                  Head ⊆ E×T    e_i → t_i
   no q(t) yet                   A1 = Head ∘ Tail^T : T → T
                                 A2 = A1 ∘ A1      : T → T
```

Both axes live inside one relational structure, and neither becomes the
other:

```text
                    G_time
                      │
      ┌───────────────┼───────────────┐
      │               │               │
      Q               T               E
   (owned,       (Tail, Head,      (Tail, Head
    named by      A1, A2 all         run from
    nothing)      run over T)        here)
```

**Q is owned, and no relation mentions it.** That is lawful: a
relational graph is a collection of member sets and typed relations over
them, and nothing requires it to be connected. Inventing a relation
merely to attach \(Q\) to the time chain would be manufacturing
structure to satisfy an aesthetic.

What the future will mean, stated now so it is not misread later:

```text
      q(t_i)  is a value on Q, indexed by an instant of T

NOT   "Q becomes T"
NOT   "state coordinates become graph vertices of the time chain"
```

## The collision, and why it is harmless

All three carriers enumerate from one, so the integer **1** is a member
of \(Q\), of \(T\) and of \(E\) — meaning \(x\), \(t_0\) and \(e_1\)
respectively. Raw integer equality is not domain identity. This is the
same discipline the partitioned tower proved on \(V, E, K\), meeting a
second, independent specimen here.

A sharper form of the same hazard lives at Level 1. Over the specimen's
numbering,

```text
Tail = { [1,1] [2,2] [3,3] [4,4] }        signature  E × T
A1   = { [1,2] [2,3] [3,4] [4,5] }        signature  T × T
```

which are, tuple for tuple, the same integers a six-vertex chain would
produce. The signature is the entire difference.

---

# D. Levels 0–4 in detail

## Level 0 — carrier

\(Q\), \(T\), \(E\), and nothing else. Cardinalities, mutual
non-identity, `member`/`local_index` round trips, and boundary refusals.
No relation, no graph, no value, no step size.

## Level 1 — relation

```text
Tail ⊆ E×T     e1→t0   e2→t1   e3→t2   e4→t3
Head ⊆ E×T     e1→t1   e2→t2   e3→t3   e4→t4
```

Signatures pinned by carrier identity. Extensions pinned exactly. Every
step has exactly one tail and exactly one head — *time's direction is
relation structure, not the order of a loop index*.

\(Q\) participates in no relation. At this level time has structure,
state coordinates exist, and **nothing says a state coordinate is a time
instant**.

## Level 2 — relation algebra

\[
A_1 = \mathrm{Head}\circ\mathrm{Tail}^{T} : T \to T,
\qquad
A_2 = A_1 \circ A_1 : T \to T
\]

read through the repository's convention
`compose_binary(R_AB, R_BC) = R_BC ∘ R_AB`.

```text
A1 = { t0→t1, t1→t2, t2→t3, t3→t4 }        one-step reach
A2 = { t0→t2, t1→t3, t2→t4 }               two-step reach
```

Both run \(T\to T\); both are derived, never written down as data; and
\(A_1 \neq A_2\) both as identities and extensionally.

### The semantic boundary this level defends

> **\(A_2\) is NOT BDF2.**

\(A_2\) says only that an instant two steps later is structurally
reachable. Later, at Level 6, a scheme such as BDF2 *may consume*
present, one-step history and two-step history — but that consumption is
an interpretation laid on the structure, not the structure itself.

```text
temporal reach   ≠   temporal discretization scheme
```

No union, no transitive closure: neither is earned by any caller yet.

## Level 3 — relational graph

\[
G_{\text{time}} = \big(\{Q,T,E\},\ \{\mathrm{Tail},\mathrm{Head},A_1,A_2\}\big)
\]

Exactly three member sets, exactly four relations. Every slot of every
owned relation resolves to a carrier \(G_{\text{time}}\) owns — the
signature validity law. `Tail`/`Head` remain \(E\times T\) and
\(A_1\)/\(A_2\) remain \(T\times T\) *after* ownership.

And the level's own truth: \(Q\) is owned while naming no relation. The
graph validates relations against its sets; it never demands that every
set be named by a relation. A relational structure need not be
connected.

## Level 4 — graph calculus

\(A_1\), interpreted as a `directed_adjacency_view` borrowed from
graph-owned storage:

```text
sources(A1)            = { t0 }
sinks(A1)              = { t4 }
reachable(t0, t4)      = true
reachable(t4, t0)      = false
topological_order(A1)  = [t0, t1, t2, t3, t4]        FORWARD CAUSAL ORDER
```

The order comes from the carrier's declaration order via `local_index`,
never from arithmetic on member values.

**Reverse causal order** — \([t_4,t_3,t_2,t_1,t_0]\) — is *stated* as
the opposite traversal of the established forward order. No new
algorithm is written for it, no adjoint is implemented, and no
derivative machinery is imported. That reverse order exists structurally
*before* an adjoint exists is the observation; building one is not this
level's business.

\(A_2\), interpreted the same way over the **same carrier**:

```text
successors_A2(t0) = { t2 }     reachable_A2(t0, t4) = true    (t0→t2→t4)
successors_A2(t1) = { t3 }     reachable_A2(t0, t3) = false   ← two-step
successors_A2(t2) = { t4 }                                      reach lands
                                                                only on even
                                                                offsets
```

That last refusal is the executable form of *reach is not a scheme*: a
relation that were BDF2's dependency would have to see \(t_3\) from
\(t_0\), and \(A_2\) does not.

The Rosetta truth a later scheme can consume:

```text
A1 and A2 are DIFFERENT STRUCTURAL VIEWS over the SAME T carrier
```

> **REVIEW GATE A** — after Level 4. The frontier stops here.

---

# E. What Gate B will ask

Stated as a question, because it has not been answered:

> When \(q\) becomes a field on \(Q\), can the existing temporal
> discretization and minimization stack preserve \(Q\) as the state
> domain, independently of whatever graph supplies structural context?

Neutral record of what will be tested there:

```text
the current step / march code (class_graph_step, class_graph_marcher)
will be exercised at Gate B; it is NOT imported at Gate A
```

Gate A imports none of it deliberately. Levels 0–4 establish time's
structure *independently of the current implementation of time
marching*; otherwise this client would merely redescribe production and
could not be evidence about it.

One further fact, recorded neutrally and **not** verified here:
`test/graph-marching/test.f90` opens by stating *"TIME IS A GRAPH: the
marcher's instants stand as a chain — one vertex per instant, one edge
per step, walked in order."* That is production's own description of
itself, read but not imported. It means the Gate-B comparison is at
least well posed — there is a chain-of-instants notion on the other side
to compare against. Whether it agrees with \(A_1\), and whether it
preserves \(Q\), is exactly what Gate B is for.

Nothing here says the production marcher is wrong, that **seam A2** of
the reverse architecture review must be closed, or that
`graph_operation` must change. Those are questions for levels that have
not been built.

> Read carefully: **seam A2** is the review's *"operations should carry
> their own domain"*, and \(A_2\) is this tower's *two-step reach
> relation*. They share a label and nothing else. This document uses
> "seam A2" for the first and \(A_2\) for the second, everywhere.

---

# F. What this tower proves / does not prove

## Proven, at Gate A

```text
Q, T and E are independent carriers; raw integer equality is not
    domain identity
temporal direction is relation structure, not loop index order
one-step and two-step reach are GENERATED from primitive incidence,
    never stored independently
a relational graph may lawfully own a state carrier that no relation
    names
causality is a graph-calculus INTERPRETATION of A1, not a property
    of A1 itself
reverse causal order exists structurally before any adjoint exists
```

## Not proven — anywhere yet

```text
anything numerical           any value of q
Forward Euler                BDF2 or any scheme
step size h                  a residual or action S(q)
a time-marching loop         the production marcher's contract
whether an operation can carry its own domain   ← SEAM A2,
                                                  NOT exercised here
```

**Seam A2** of
[the reverse architecture review](../../doc/REVERSE-ARCHITECTURE-REVIEW.md)
— *operations should carry their domain rather than ask a graph for one*
— is **not** advanced by Gate A. Nothing here has attached an operation
to anything.

---

# G. Code map

```text
test/time-integration-tower/
├── README.md                    this document — the Rosetta stone
├── NUCLEUS-OBSERVATIONS.md      the evidence ledger (TI-*), by level
├── run.sh                       level-by-level runner; stops at Gate A
├── check_imports.sh             fail-closed allowlists, PER LEVEL,
│                                + its own --selftest
├── common/
│   ├── time_assert.f90                    (below everything)
│   ├── time_carriers_fixture.f90          earned at Level 0
│   ├── time_relations_fixture.f90         earned at Level 1
│   └── time_algebra_fixture.f90           earned at Level 2
├── level-0-carrier/             test.f90
├── level-1-relation/            test.f90
├── level-2-relation-algebra/    test.f90
├── level-3-graph/               test.f90
└── level-4-graph-calculus/      test.f90

    level-5-field-calculus/      NOT YET
    level-6-discretization/      NOT YET
    level-7-minimization/        NOT YET
    level-8-constitution/        NOT YET
    level-9-statement/           NOT YET
```

The fixture ladder is the tower's own stratification applied to itself:

```text
Level 0    time_carriers_fixture      declares Q, T, E
Level 1    time_relations_fixture     states Tail, Head over them
Level 2    time_algebra_fixture       composes what follows
```

The relation fixture does not *import* the carrier fixture — its
constructors receive \(Q,T,E\) as arguments, because a Level-1 file may
state facts over sets but may not name a set into existence. The ladder
is enforced by the import gate's per-file allowlists, by each level's
Makefile, and by `check_imports.sh --selftest`, which asserts that a
Level-0 source saying `use time_relations_fixture` is refused.

There is no `check_marker.sh`: Gate A produces no numerical result, and
a result contract for a tower with no result would be a claim it has not
earned.
