# Reverse Architecture Review

**Mode:** reverse — evidence flowing from the sealed orbital towers back
toward the nucleus, per
[FRACTAL-GRAPH-ARCHITECTURE.md](../FRACTAL-GRAPH-ARCHITECTURE.md) §15.

**Status at review:** `tower-graph-as-sets-relations @ abb9167`.
Four towers sealed. **No production change is proposed by this document
for immediate action.**

> The governing rule, restated because this document exists to obey it:
>
> \[
> \boxed{\text{evidence, not elegance, changes the nucleus}}
> \]
>
> A repeated inconvenience is not automatically the wrong abstraction,
> and a high-radius ownership need must not be destroyed because a
> radius-0 operation ignores something numerically.

---

# 1. Method

## 1.1 Evidence thresholds

The repository's own scale (fractal document §16):

```text
1 tower                    observation
2 unrelated towers         recurring seam
3+ independent towers      strong architectural evidence
high-radius failure        high-value evidence
several high-radius        nuclear refactor candidate
```

**Counting rule applied here.** A tower counts *once*. Multiple levels,
gates or radii within one tower do **not** multiply into independent
towers. Where a seam appears at several radii of a single tower, that
strengthens *confidence in the observation*, not the tower count.

**Second counting rule, added by this review.** Evidence that a contract
is *unused by a caller* is not evidence that the contract is *wrong*.
The two must be separated explicitly, because the strongest-looking seam
in the ledgers fails exactly on this distinction (§3.5).

## 1.2 Sources reviewed

| Source | Observations | Note |
|---|---|---|
| `test/derivative-action-tower/NUCLEUS-OBSERVATIONS.md` | 24 (DA-0 … DA-9C) | full ledger |
| `test/adjoint-tower/NUCLEUS-OBSERVATIONS.md` | 16 (AD-0 … AD-15) | full ledger |
| `test/calculator-tower/` | — | sealed before the fractal document existed; evidence read from its README, tests and commit record |
| `test/learning-tower/` | — | as above |
| `src/` production contracts | — | traced directly for this review |

The calculator and learning towers predate the ledger discipline. Their
evidence is real but was never written in ledger form; this review reads
it from their code and READMEs and says so wherever it is used.

## 1.3 The four sealed towers

| Tower | Result | Levels inhabited | Production changes |
|---|---|---|---|
| Calculator | \(20\) | 0–9 | NONE |
| Learning | \(w=3\) | 0–9 | NONE |
| Derivative Action | \([3,3]\) on \(X\) | 0–6, 8, 9 (**7 N/A**) | NONE |
| Adjoint Sensitivity | \(df/dp=7\) | 0–9 | NONE |

Four consecutive clients closed with zero production changes. That is
itself the review's most important background fact: **the nucleus is not
currently obstructing its clients.** Nothing below should be read as a
complaint.

---

# 2. The recurring seams, named

Three candidates were nominated by the ledgers. A fourth, minor, is a
contract-wording defect.

```text
A   the graph_operation host                 (DA-8F, DA-9A, AD-9, AD-13)
B   bidirectional linearization              (DA-8B, AD-10, AD-12, AD-15)
C   model-owned domains vs domain handles    (AD-13, and §4 of the brief)
D   marker contract wording                  (defect, not architecture)
```

---

# 3. CANDIDATE A — the `graph_operation` host

## 3.1 The contract as it stands

```fortran
! graph_grammar / graph_calculus
subroutine apply (this, input_graph, input_data, output)
subroutine domain(this, input_graph, domain)
   class(graph), intent(in) :: input_graph
```

Every operation family inherits it:

```text
graph_operation
├── differential_operator          class_graph_differential_operator.f90
├── balance                        class_graph_balance.f90
├── fit                            graph_fitting.f90
├── walk                           class_graph_walk.f90
├── graph_reduction  (abstract)    graph_calculus.f90
├── graph_broadcast  (abstract)    graph_calculus.f90
├── discretization_operator (abs)  graph_calculus.f90
├── linearization_operator  (abs)  graph_calculus.f90
│     └── difference_linearization class_graph_linearization.f90
└── minimizer        (abstract)    graph_minimization.f90
```

## 3.2 What each implementation actually reads

Measured, not assumed — by counting genuine topology reads
(`num_edges`, `edge_tail`, `edge_head`, `edge_has_head`, `incidence`,
`breadth_first`, `components`, `edge_set`) against mere domain reads
(`vertex_set`, `num_vertices`, `all_vertices`):

| Implementation | topology reads | domain reads | **Class** |
|---|---:|---:|---|
| `differential_operator` | 28 | 6 | **1 — topology** |
| `walk` | 6 | 5 | **1 — topology** |
| `balance` | 4 | 2 | **1 — topology** |
| `graph_reduction` / `graph_broadcast` | 2 | 1 | **2 — domain only** |
| `fit` | 0 | 3 | **2 — domain only** |
| `difference_linearization` | 0 | 4 | **2 — domain only** |
| `minimizer` | 0 | 0 | **3 — ignores it** |

This is a genuine three-way factorization, and it is the review's
central finding:

```text
CLASS 1   needs the graph AS A GRAPH
          it traverses edges, incidence, adjacency
          differential_operator, balance, walk

CLASS 2   needs a DOMAIN, and asks the graph for one
          all_vertices() / vertex_set() and nothing else
          fit, reduction, broadcast, difference_linearization

CLASS 3   needs nothing; carries the argument regardless
          minimizer  ->  `associate (u1 => input_graph); end associate`
```

## 3.3 What the towers observed

| Tower | Radius | Exact caller | Observation |
|---|---|---|---|
| Calculator | L7, L9 | `solver % attach(oracle, host, u)` with a 7-vertex host | host unrelated to \(U\), never read |
| Learning | L7, L9 | same shape, 7-vertex host | recorded in README as "compatibility scenery" |
| Derivative Action | Gate B | its evaluators take **no** graph at all | DA-8F: graph not an operand at evaluator radius |
| Derivative Action | Gate C | statement owns structure | DA-9A: graph **is** the ownership environment |
| Adjoint | Gate B | `attach(eq, host, Q)` / `(eq, host, Y)`, 5-vertex host | AD-9: host provably neither \(Q\) nor \(Y\) nor their size |
| Adjoint | Gate C | model graph **and** host coexist | AD-13: two graph roles at once |

**Tower count for "the minimizer does not read its host": three**
(calculator, learning, adjoint) → **strong evidence** by the
repository's threshold.

**Tower count for "a radius-0 numerical operation needs a domain, not a
graph": two** (derivative action, adjoint) → **recurring seam**, with
three independent *production* corroborations (`fit`, `reduction`,
`difference_linearization` are all Class 2).

## 3.4 The graph's four roles, separated

The brief asked for this separation and the towers supply it:

| Role | Where it is real | Where it is absent |
|---|---|---|
| **numerical operand** | Class 1: `differential_operator` reads incidence to compute | Classes 2 and 3 |
| **compatibility argument** | the minimizer's `on`, from the perspective of the *minimizer itself* | — |
| **contextual environment** | multigrid's `this % on`; a level's real graph | test-local towers |
| **structural owner** (`relational_graph`) | adjoint Gate C model; learning L9; derivative Gate C | never the same object as the legacy host |

Note the last row: `relational_graph` and the legacy `class(graph)` host
are **different types holding different things**, and the adjoint
statement holds both simultaneously (AD-13). No proposal below conflates
them.

## 3.5 The decisive counter-evidence

The obvious refactor — *remove `class(graph)` from `graph_operation`,
since three towers show it unread* — is **REJECTED**, on production
evidence the towers could not see.

The minimizer ignores the graph **itself**, but it is the **conduit** by
which its attached action receives one. Traced to real callers:

```fortran
! test/graph-minimization/test.f90:345
call cg % attach(vertex_differential_operator(order=2,          &
     &           coefficients=c, spacings=deltas), m, m % vertex_set())
!                                                  ^^^ the real mesh

! src/class_graph_multigrid.f90:123
call this % smoother % attach(this % action, this % on, this % unknown_domain)
!                                            ^^^^^^^^^ the level's own graph

! src/graph_fitting.f90:233
call solver % attach(dual, dual % pattern, dual % pattern % vertex_set())
!                          ^^^^^^^^^^^^^^ the fit's pattern graph
```

The attached action there is a **Class-1** operation. The solver must
carry a graph precisely so it can hand the mesh to an operation that
genuinely traverses it.

All four towers attached **Class-2/3** actions — test-local operations
that compute from their input field alone. That is why the host looked
gratuitous. It is gratuitous *for their actions*, and load-bearing for
the CFD path.

\[
\boxed{
\begin{gathered}
\text{"three towers never read the host" is TRUE}\\
\text{and does NOT imply "the host is unnecessary"}
\end{gathered}
}
\]

This is the fractal document's **Case III** (\(\neg L \wedge G\)) in the
flesh, and it is exactly what §27.1 warns against destroying.

## 3.6 What the evidence *does* support

The seam is not the graph argument. It is narrower and sharper:

> **Class-2 operations obtain their domain *from* the graph instead of
> carrying it.**

`fit`, `reduction`, `broadcast` and `difference_linearization` all
answer `domain()` by calling `input_graph % all_vertices()` or
`% vertex_set()`. Consequences, both observed:

- they cannot serve a client whose domain is not a graph's vertex set —
  which is exactly why `difference_linearization` cannot express
  \(R_q : Q \to Y\) (AD-12);
- they force a graph to exist wherever they are used, even when the
  mathematics has only domains.

**And the nucleus has already solved this once.** The `minimizer` carries
`unknown_domain` and `residual_domain` explicitly — earned during the
learning tower's Phase 5B — and *because* it does, the adjoint tower
could exchange the orientation \(Q\to Y\) / \(Y\to Q\) with no solver
change (AD-7). The precedent exists, is in production, and works.

**Proposed mathematical contract** (recorded, not implemented):

```text
an operation is a map between DOMAINS:

    domain_in  : member_set          (what it consumes)
    domain_out : member_set          (what it answers on)
    apply(input_data, output)

whether it additionally needs a GRAPH is a separate question,
answered by the operation, not by the interface
```

Class 2 would then carry the domain it is given; Class 1 would keep the
graph it genuinely traverses; Class 3 would keep a graph only as an
explicit *conduit* for a possibly-Class-1 action — which is honest
context, not debt.

**Verdict: BUILD ANOTHER TOWER FIRST.** See §7.

---

# 4. CANDIDATE B — bidirectional linearization

## 4.1 The recurring mathematical contract

Written by hand, from scratch, in two independent towers:

```text
input domain    U
output domain   Y

forward         L  : U  -> Y
reverse         L* : Y* -> U*

one underlying law generating both
```

with, in every instance: no assumption that \(U=Y\); no assumption that
either is a graph's vertex set; no assembled matrix.

| Tower | Radius | Use |
|---|---|---|
| Derivative Action | Gate B | evaluation only — \(Jv\) and \(J^{T}\bar z\) over an execution order |
| Adjoint | Gate B | inside a **solve** — \(R_q^{T}\lambda=f_q^{T}\) |
| Adjoint | Gate C | composed into a complete statement; \(f_q^{T}\) itself generated by the reverse reading |

**Tower count: two** → **recurring seam**, stressed at three radii.
Radii do not multiply the count (§1.1).

## 4.2 Why `difference_linearization` is narrower

Precisely, from `src/class_graph_linearization.f90`:

```fortran
state = field('state', input_graph % vertex_set(), ncomp=ncomp)   ! ~line 156
state = field('state', input_graph % vertex_set(), ncomp=ncomp)   ! ~line 162
state = field('J v'  , input_graph % vertex_set(), ncomp=ncomp)   ! ~line 170
```

Three narrowings, each independent:

1. **Same-domain.** Input and output are the *same* set, so it can only
   express \(L : U \to U\). \(R_q : Q \to Y\) with \(Q \neq Y\) is
   inexpressible.
2. **Graph-as-domain.** That set is a graph's vertex set, so the domain
   cannot be an arbitrary `member_set` — it is Class 2 (§3.2).
3. **Forward-only.** It exposes `apply` (a \(Jv\) product by finite
   difference) and no reverse action at all. The transpose is not
   narrow here; it is *absent*.

It is a legitimate specialized citizen — a finite-difference \(Jv\) on a
mesh — and this review proposes no change to it.

## 4.3 Is the recurrence a nucleus contract or application law?

Distinguishing carefully, because the answer decides everything:

| Duplicated across towers | Kind | Healthy? |
|---|---|---|
| coefficients \(A,b,c,d\); \(u=xy,\ z=u+y\) | **application law** | yes — specimens differ, of course they differ |
| "one law, read forward and reverse" | **protocol shape** | recurring |
| explicit input/output domains | **protocol shape** | recurring |
| support-driven traversal of incidences | **protocol shape** | recurring |
| `+=` accumulation on the reverse side | **protocol shape** | recurring |

So the recurrence is *not* merely application law. But two towers is a
recurring seam, not strong evidence — and both towers are in the **same
orbital family** (derivative). The fractal document's threshold asks for
*independent* towers; two members of one family are weaker than two
unrelated clients.

**Not added, deliberately:** `apply_transpose`, `adjoint_operation`,
`bidirectional_operator`. Naming a thing is not the same as needing it.

**Verdict: BUILD ANOTHER TOWER FIRST.** See §7.

---

# 5. CANDIDATE C — model-owned relations, externally-held domains

## 5.1 The observation

The sealed adjoint statement destroys its relation selectors and
continues through model-owned relations — but it goes on using the
external `p_dom`, `q_dom`, `y_dom`, `z_dom` handles for every field and
every `local_index`.

Is that an incomplete ownership boundary?

## 5.2 Classification: **healthy identity handle**

`member_set` identity is a **value-semantic token**: `held_set` copies
the set into the graph, and the copy is `same_as` the original because
identity travels with the value (`graph_identity`'s opaque
image+serial). `model % holds_set(q_dom)` answers true for the external
handle, and answers false for the host's vertex set.

So a domain does *not* need to be fetched from its owner to be the same
domain. A relation does — which is why the relation lifetime test is
meaningful and a domain lifetime test would be theatre.

```text
relations  are OWNED objects        -> fetch from the owner (proved)
domains    are IDENTIFIED values    -> a handle IS the domain
```

The asymmetry is not an oversight; it follows from the identity design
the towers repeatedly leaned on.

## 5.3 What a larger-radius client would want

Recorded as an open question, not a defect. A partitioned or
multiphysics statement will hold *many* domains and will want to
**discover** them (“give me this model's residual domain”) rather than
carry handles for all of them. That is a **lookup/context API** question
— `relational_graph` currently exposes `member_set_at` and
`holds_set`, which is enough to build such a lookup locally, exactly as
the towers build every other convenience locally.

**Verdict: REJECT the change; KEEP the current design.** Nothing is
broken. Revisit only when a client actually needs domain *discovery*
rather than domain *possession*.

---

# 6. CANDIDATE D — the marker contract wording (fixed)

Not architecture; a truthfulness defect in two sealed towers' result
checkers, and repaired as part of this review.

Both `check_marker.sh` scripts claimed to validate a **finite** real
token. The grammar

```text
^[+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eEdDqQ][+-]?[0-9]+)?$
```

is **lexical**: it excludes `NaN` and `Inf` spellings because they
contain letters the grammar does not admit, but it accepts `1e999999`,
which parses to infinity. "Finite" was an overclaim about a syntactic
check.

**Fixed** by making the wording honest — the checkers now say they
validate a *real numeric literal*, syntactically, and state explicitly
that non-finite spellings are excluded by that syntax while overflowing
literals are not caught. No behaviour changed; both towers remain green.

---

# 7. REVERSE EVIDENCE MATRIX

| # | Seam | Towers | Independent? | Radii | Exact caller / contract | Local nec. | Global nec. | Graph role | Current workaround | Evidence grade |
|---|---|---|---|---|---|---|---|---|---|---|
| A1 | minimizer does not read its graph host | Calculator, Learning, Adjoint | **3, yes** | 1 | `graph_minimization.f90` `solver_apply` (`associate`-and-discard) | no | **YES** — it is a conduit to Class-1 actions | compatibility / conduit | supply an unrelated host | **strong evidence for the FACT; the fact does not imply the refactor** |
| A2 | Class-2 operations take their domain from the graph | Derivative Action, Adjoint (+3 production sites) | 2 towers, same family | 0–1 | `fit`, `graph_reduction`, `graph_broadcast`, `difference_linearization` → `vertex_set()` / `all_vertices()` | yes | unknown | graph-as-domain | test-local operations carrying explicit domains | **recurring seam** |
| A3 | `relational_graph` is load-bearing as owner | Learning, Derivative Action, Adjoint | **3, yes** | 1 | L9 statements; `relation_at` + `same_as` + selector destruction | yes | unknown | structural owner | none needed — it works | **strong evidence that the design is RIGHT** |
| B | one law, forward and reverse, between different domains | Derivative Action, Adjoint | 2, one family | 0–1 | test-local fixtures; `difference_linearization` is narrower on three axes | yes | unknown | none | test-local adapters | **recurring seam** |
| C | domains held externally while relations are owned | Adjoint | 1 | 1 | Gate C statement | n/a | n/a | n/a | none — identity handles are sound | **observation; design is correct** |
| D | "finite real" wording in marker checkers | Derivative Action, Adjoint | 2 | n/a | `check_marker.sh` | n/a | n/a | n/a | — | **defect, fixed** |

## ACT NOW

```text
D   marker contract wording          — done in this review, no architecture
```

Nothing else clears the bar. In particular **A1 does not**, despite
being the only seam with three independent towers: the towers establish
that the minimizer never reads its host, and production establishes that
the host is how a `differential_operator` receives its mesh. Strong
evidence for a true fact that does not license the refactor it appears
to license.

## BUILD ANOTHER TOWER FIRST

```text
A2  operations should CARRY their domain rather than ask a graph for one
B   a bidirectional linearization contract with explicit input/output
    domains and a reverse action
```

Both are recurring seams from a single orbital family. Both would be
settled by **one** experiment, described in §8.

## REJECT / KEEP CURRENT DESIGN

```text
A1  removing class(graph) from graph_operation
    REJECTED — Case III; it is the conduit for Class-1 actions

A3  the relational_graph ownership contract
    KEEP — three towers show it working; no change proposed

C   fetching domains from their owning graph
    KEEP — member_set identity is value-semantic; a handle IS the domain
```

## Requirements not met by any ACT-NOW candidate

For completeness, the bar an architectural change must clear here, and
which none of A2/B currently meets:

| Requirement | A2 | B |
|---|---|---|
| ≥ strong evidence (3 independent towers) | ✗ (2, one family) | ✗ (2, one family) |
| precise mathematical replacement contract | ✓ (§3.6) | ✓ (§4.1) |
| known affected production files | ✓ `graph_fitting`, `graph_calculus`, `class_graph_linearization` | ✓ `class_graph_linearization`, `graph_calculus` |
| migration strategy | partial — Class 1 untouched, Class 2 gains a stored domain | not designed |
| acceptance towers that protect it | ✓ all four sealed towers + 17 suites | ✓ the two derivative-family towers |
| performance implications | unmeasured; `graph-benchmark` exists | unmeasured |

---

# 8. The next experiment

One tower would settle both A2 and B, and it must be **outside the
derivative family** to count as independent evidence.

```text
proposed: a PARTITIONED or COUPLED tower

why it discriminates:
  · its domains are built by partitioning, not by a graph's vertex set
      -> forces A2 to be true or false, not merely inconvenient
  · it has real topology at the coupling boundary
      -> exercises CLASS 1 in the same program as CLASS 2,
         which no tower has yet done
  · it will attach a Class-1 action to a minimizer
      -> the first tower to see the host as a CONDUIT rather than
         as scenery, testing §3.5 directly
  · if it also needs a linearization between two different domains,
      B reaches three independent towers and clears the bar
```

Until then the nucleus stands unchanged, and rightly so: four clients
have inhabited it without asking it to move.

---

# 9. Summary for the architect

```text
strongest seam by tower count      A1 (3 towers) — and REJECTED on
                                   production evidence
strongest seam still open          A2 / B (2 towers each, one family)
most confirmed design decision     A3 relational_graph ownership, and
                                   the minimizer's explicit
                                   unknown/residual domains
changed by this review             documentation only (candidate D)
recommended next action            one non-derivative-family tower,
                                   not a nucleus mutation
```

The single most useful thing this review found is not a refactor. It is
that the seam with the most tower-votes was the one the towers were
least equipped to judge — because all four attached operations that had
no topology to traverse. That is the value of reverse mode: it is where
a locally-obvious conclusion meets the callers that would have paid for
it.
