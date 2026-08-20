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

    L5 field calculus .................. PASS
    L6 discretization .................. PASS
    L7 minimization .................... PASS

    ===== REVIEW GATE B =====

    L8 constitution .................... PASS
    L9 statement ....................... PASS

    ===== REVIEW GATE C =====

    TIME_INTEGRATION_RESULT = 2.9166666666666674E-01 5.7638888888888895E-01

    TOWER SEALED.
```

All ten levels are built. \(q(t_4) = [7/24,\, 83/144]\), on \(Q\).

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
| **3** | \(G_{\text{time}}=(\{Q,T,E\},\{\mathrm{Tail},\mathrm{Head},A_1,A_2\})\) | `graph_fractal`, `view_relational` | three bound carriers | four bound relations | one relational structure, not one domain | `level-3-graph/` | signature closure holds, **and** \(Q\) is lawfully owned while naming no relation | none |
| **4** | directed views over \(A_1\), \(A_2\) | `directed_adjacency_view`, `graph_profile`, `relation_algorithms` | \(T\) (both views) | borrowed from the graph | **causality** | `level-4-graph-calculus/` | forward causal order \([t_0..t_4]\); two views, one carrier | none |
| **5** | \(q_0:Q\to\mathbb R\); \(\mathrm{time}:T\to\mathbb R\); \(h:E\to\mathbb R\) | `field` | \(Q\), \(T\), \(E\) — three distinct | consumed, not made | values, no scheme | `level-5-field-calculus/` | **values live on domains, not on graphs** — \(q_0\) needs no graph, and \(\mathrm{time}(\mathrm{head}(e))-\mathrm{time}(\mathrm{tail}(e))=h(e)\) | none |
| **6** | \(S:Q\to Q\); FE; BE residual; BDF2 residual | `step_operator`, `backward_euler`, `bdf` | \(Q\), beside a 5-vertex host \(H_t\) | \(A_1,A_2\) supply the history roles | the discrete law | `level-6-discretization/` | **temporal discretization preserves \(Q\)** — RED first, then a narrow correction | `class_graph_step.f90` — domain and width delegated to the action |
| **7** | solve temporal residual \(=0\) | `gmres` | unknown \(Q\), host \(H_t\) | — | the implicit solve | `level-7-minimization/` | **an explicit unknown domain \(Q\) survives minimization** while the host has five unrelated vertices | none — the minimizer was already right |
| **8** | four-step marches: FE, BE, BDF2 | `marcher`, `newton`, `gmres` | \(Q\); control chain; \(H_{\text{context}}\) | control chain realizes \(A_1\) extensionally | multi-step constitution | `level-8-constitution/` | **the marcher's state domain is independent of its host**; three five-element objects stay apart | `class_graph_marcher.f90`, `class_graph_linearization.f90` — both RED-earned |
| **9** | the complete IVP: \(q(0)=[2,0]\), \(\dot q=-S(q)\), \(h=1/2\), BDF2 → \(q(t_4)\) | the whole constitution | \(Q\), in and out | — | the statement | `level-9-statement/` | \(q(t_4)=[7/24,\,83/144]\), **a field on \(Q\)** | none |

The road as far as it has been built:

```text
L0  three carriers, three identities
L1  primitive temporal incidence — time gets a direction
L2  derived one-step and two-step reach
L3  one relational structure owning both axes without conflating them
L4  causal traversal: sources, sinks, reachability, topological order
--- REVIEW GATE A ---
L5  values on Q, T and E — three domains, no graph
L6  the scheme: S on Q, and a discretization that keeps Q
L7  the implicit solve, on Q, beside a five-vertex host
--- REVIEW GATE B ---
L8  the multi-step march, constituted from every earned piece
L9  the complete initial-value problem, answered on Q
--- REVIEW GATE C ---   TOWER SEALED
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
`compose_binary(P_AB, P_BC) = P_BC ∘ P_AB`.

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

> **REVIEW GATE A** — after Level 4.

---

# D2. Levels 5–7 — values, scheme, solve

## The road, and the arrow that is never drawn

```text
        TIME STRUCTURE                      STATE DOMAIN
        T, E, A1, A2                            Q
                                                |
        (structure: which instants,             |  q_n : Q -> reals
         in which order, reachable              |
         from which)                            v
                                    +-----------------------+
                                    | temporal step operator|
                                    |  a0 q + a1 qold       |
                                    |  + a2 qolder + h S(q) |
                                    +-----------------------+
                                          ^          |
                    H_t --------------------          |
                    (compatibility host,              |
                     5 vertices, carried              v
                     as context, unread)        residual on Q
                                                      |
                                                      v
                                              GMRES, unknown = Q
                                                      |
                                                      v
                                                q_{n+1} : Q -> reals
```

**No arrow in that picture says \(Q = V(H_t)\), and none may.** The host
enters from the side because the `graph_operation` contract requires
one; it leaves nothing behind. \(Q\) has two members and \(V(H_t)\) has
five, throughout.

## Level 5 — field calculus

Three fields on three domains:

```text
q0   : Q -> reals  [2, 0]
time : T -> reals  [0, 1/2, 1, 3/2, 2]
h    : E -> reals  [1/2, 1/2, 1/2, 1/2]
```

pinned by domain **identity**, and pairwise distinct. The central truth:

\[
\boxed{\text{STATE VALUES LIVE ON } Q \text{, INDEPENDENTLY OF ANY GRAPH.}}
\]

No graph is constructed at Level 5. None is needed — production field
calculus already says a field is a function over one `member_set`, and
\(Q\) is one. That capability is not discovered here; it is *exercised*
by a client whose state domain is emphatically nobody's vertex set.

Values are not structure. The instant \(t_2\) is a **member** of \(T\);
the real \(1.0\) is a **value** at \(t_2\). Their consistency is proved,
not assumed, through the relations Level 1 earned:

\[
\mathrm{time}(\mathrm{head}(e)) - \mathrm{time}(\mathrm{tail}(e)) = h(e)
\qquad\text{for every step } e.
\]

## Level 6 — discretization, and the seam

The first rung to touch production, and the first place **seam A2** is
genuinely exercised. The action \(S([x,y]) = [x, y-x]\) carries its own
domain and stores no graph. The compatibility host \(H_t\) is a
five-vertex chain:

```text
|V(H_t)| = 5        |Q| = 2        V(H_t) is NOT Q, and NOT T
```

The mismatch is load-bearing: at equal cardinality a substitution of one
carrier for the other would produce plausible numbers and the seam would
hide.

```text
S(q0)    = [2, -2]                          direct action, no fix needed
q_FE,1   = q0 - h S(q0) = [1, 1]            by hand, not by marcher
q_BE,1   = [4/3, 4/9]                       residual zero, by substitution
q_BDF2,2 = [5/6, 47/72]                     residual zero, by substitution
```

**RED came first.** Asked `step % domain(H_t, d)`, the production
reviewed at Gate A answered the host's five vertices, and `step % apply`
then died inside `set_real_vector`. The failure is recorded verbatim in
[`NUCLEUS-OBSERVATIONS.md`](NUCLEUS-OBSERVATIONS.md) TI-8, written
before production was touched.

The correction is narrow and lives in `src/class_graph_step.f90` alone:

\[
\boxed{\text{TEMPORAL DISCRETIZATION PRESERVES THE DOMAIN OF THE ACTION IT DISCRETIZES.}}
\]

A step is an operation *built from* another operation, so its residual
is a statement about the same unknown. The domain question is delegated
to the action; the component width is read from the input field; the
answer lands on the action's domain. For every graph-based action —
which is all of them on the ordinary-graph road, since each answers
`input_graph % all_vertices(...)` — this returns exactly what asking the
graph returned, and `test/graph-marching` passes unchanged.

### Structure and scheme, finally joined

Level 2 derived \(A_1\) and \(A_2\) and *refused* to call either a
scheme. Level 6 supplies the other half. At instant \(t_2\):

```text
A1-predecessor of t2  =  t1        one-step history role
A2-predecessor of t2  =  t0        two-step history role

bdf-2:   a0 = 3/2 at t2     a1 = -2 at t1     a2 = 1/2 at t0
```

\[
\boxed{A_1, A_2 \text{ supply STRUCTURAL REACH; the scheme supplies the NUMBERS.}}
\]

Neither contains the other, which is why Level 2 was right to refuse the
name.

## Level 7 — minimization

Production GMRES, unknown domain \(Q\), host \(H_t\) with its five
unrelated vertices. The affine constant is measured, not assumed, and
the implicit step becomes a linear system:

```text
backward euler   c = R(0) = -q0 = [-2, 0]
                 rhs = [2, 0]        ->  q1 = [4/3, 4/9]

bdf-2            c = -2q1 + q0/2
                 rhs = [5/3, 8/9]    ->  q2 = [5/6, 47/72]
```

Both right-hand sides are fields **on \(Q\)**; both solutions come back
**on \(Q\)**.

**The minimizer needed no change.** `graph_minimization` already takes
its unknown domain as an explicit argument and asks the *action* for the
residual domain — its own comment says "no hidden fallback to the host's
vertices." That contrast is the evidence: the same seam that was open in
the discretization was already closed in the minimization, by a contract
written the other way round.

> **REVIEW GATE B** — after Level 7.

---

# D3. Levels 8–9 — constitution and statement

## Four different things are called "time", and they are not synonyms

This is one of the tower's main Rosetta contributions, because every
one of these has been called "the time" by somebody:

```text
1.  TIME MEMBERS            T = { t0 t1 t2 t3 t4 }
                            a carrier.  Five identities, no order,
                            no numbers.                        (L0)

2.  TEMPORAL RELATION       A1 = Head ∘ Tail^T,  A2 = A1 ∘ A1
                            which instants follow which.
                            Structure, no coefficients.        (L2)

3.  NUMERICAL COORDINATE    time : T → R = [0, ½, 1, 3/2, 2]
                            a FIELD on T.  Values, not structure,
                            and h : E → R beside it.           (L5)

4.  MARCHING RULE           FE / BE / BDF2, constituted through
                            the marcher.  Coefficients weighting
                            the roles (2) supplies.        (L6, L8)
```

And two more objects that are neither, and must not be confused with
each other:

```text
    TIME GRAPH              the marcher's CONTROL CHAIN, generated
                            by clock % instants(4).  Five vertices.

    OPERATION HOST          H_context, the conduit the
                            graph_operation contract requires.
                            Five vertices — and a different object.
```

At Level 8 all three five-element things exist in one program —
\(T\), the control chain, \(V(H_{\text{context}})\) — and the level
pins that no two are the same carrier. **That they all have five
elements is a coincidence of this specimen**, not a fact about time,
and the assertions refuse to lean on it. Beside them sits the
two-member \(Q\), which is none of them.

## The constitution road

```text
    G_time / A1                      Q
        │                            │
        │ extensional realization    │  q0 : Q → R
        ▼                            ▼
    marcher.instants(4)        triangular_decay  S : Q → Q
        │                            │
        │  four causal steps         ▼
        │                      step scheme : Q → Q
        │                            │
        │                            ▼
        │                         Newton
        │                            │
        │                            ▼
        │            difference linearization : Q → Q
        │                            │
        │                            ▼
        │                          GMRES
        │                            │
        │                            ▼
        └──────────────────────► q_(n+1) : Q
                                     │
                    repeated four times along the chain
                                     │
                                     ▼
                                 q(t4) : Q

    H_context ─────────────────────────────► carried alongside,
    (5 vertices, chain)                      passed to the action,
                                             never the state's seat
```

\[
\boxed{\begin{gathered}
\text{ONE STATE DOMAIN } Q \text{ SURVIVES THE WHOLE ROAD:}\\
\text{field} \to \text{action} \to \text{step} \to \text{Newton}
\to \text{linearization} \to \text{GMRES} \to \text{marcher}\\[2pt]
\text{while the time structure and the compatibility host}\\
\text{remain separate objects throughout.}
\end{gathered}}
\]

## Level 8 — the marches, and two REDs

The control chain first, before any value moves. Production's
`instants(4)` is checked **extensionally** against \(A_1\) — step
\(i\) joins instant \(i\) to instant \(i+1\), exactly the pairs
\(A_1\) holds — and *not* by identity, which is not required and
would be wrong to demand. Two realizations of one structure; agreement
never made two parties one party.

Then the marches:

```text
FE    [2,0] → [1,1] → [1/2,1] → [1/4,3/4] → [1/8,1/2]
BE    [2,0] → [4/3,4/9] → [8/9,16/27] → [16/27,16/27] → [32/81,128/243]
BDF2  [2,0] → [4/3,4/9] → [5/6,47/72] → [1/2,2/3] → [7/24,83/144]
```

every prefix pinned, not only the terminal state. BDF2's first state
*is* the backward-euler one — a two-step scheme cannot reach two steps
back on its first step, which is a structural fact before it is a
numerical one.

**Two REDs, both recorded before production was touched.**

```text
RED 1   read_statement built the state on the HOST's vertex set and
        took its width as size(q)/num_vertices() = 2/5 = 0.
        → class_graph_marcher.f90, narrow fix                 TI-14

RED 2   difference_linearization built its frozen and perturbed
        states on the HOST's vertex set.
        → class_graph_linearization.f90, narrow fix           TI-16
```

The second is the one that matters as evidence. It was reached through

```text
marcher → newton → difference_linearization → gmres
```

— the production path an implicit march *requires* — and **no level of
this tower may import that module**; the import gate refuses it
universally and `--selftest` asserts the refusal. So the Class-2
witness is demonstrably natural rather than manufactured.

The implicit governor is **Newton, not bare GMRES**: the marcher drives
the whole residual to zero via `inner % solve(zeros, q, …)`, and a bare
GMRES `matvec` has already subtracted the affine constant, which would
answer a different question.

### Two specializations, and they are not defects

```text
production regenerates a LINEAR CHAIN from nsteps
    rather than consuming G_time

production carries ONE SCALAR step
    rather than the field h : E → R
```

For this specimen both are **exact specializations** — the time graph
is a simple chain and \(h\) is uniform, and Level 8 checks
\(h(e) = \text{clock}\%\text{step}\) at every step. So no defect is
established. They are recorded as frontier in §F below, for clients
that would supply a nonuniform or nonlinear time structure. One tower
cannot decide that.

## Level 9 — the statement

> Given \(Q=\{x,y\}\), \(q(0)=[2,0]\), \(\dot q = -S(q)\) with
> \(S(q)=[x,\,y-x]\), \(t_0=0\), \(t_4=2\), \(h=1/2\), and BDF2
> with one backward-euler startup step — compute \(q(t_4)\).

Both ends of the statement are **fields on \(Q\)**: the initial state
arrives as one and its vector is fetched once; the answer is written
back into one and its domain checked. The marcher's raw-array core was
left exactly as it is — nothing was refactored to make a public
argument prettier.

The endpoint is *earned*, not assumed: \(\text{time}(t_4)=2\) read off
the Level-5 coordinate field, four steps walked, and the control chain's
terminal instant reached by **following its incidence from \(t_0\)**
rather than by indexing the integer 5.

\[
\boxed{q(t_4) = \left[\tfrac{7}{24},\ \tfrac{83}{144}\right]}
\]

```text
TIME_INTEGRATION_RESULT = 2.9166666666666674E-01 5.7638888888888895E-01
```

Two tokens, not five — the result lives on \(Q\), and a five-token
marker would quietly undo nine levels of argument at the last line of
output. `check_marker.sh` validates shape and syntax only; whether those
numbers *are* \(q(t_4)\) is the Level-9 test's business.

> **REVIEW GATE C** — after Level 9. **Tower sealed.**

---

# E. What Gate B asked, and what it answered

The question, as Gate A posed it:

> When \(q\) becomes a field on \(Q\), can the existing temporal
> discretization and minimization stack preserve \(Q\) as the state
> domain, independently of whatever graph supplies structural context?

**Answered: yes — after one narrow correction, and the two halves of the
stack answered differently.**

```text
graph_minimization      ALREADY preserved Q.  Unknown domain explicit,
                        residual domain asked of the action.  No change.

class_graph_step        DID NOT.  It read the host's vertices for the
                        domain and the host's vertex count for the
                        width.  RED at Level 6; corrected there.
```

The marcher was **not** imported and remains untested. Whether the
machinery that stamps a step along a chain can carry a state domain that
is not its host's vertex set is a question for a level that does not
exist, and the import gate refuses the module at every built level.

One further fact, recorded neutrally and **not** verified here:
`test/graph-marching/test.f90` opens by stating *"TIME IS A GRAPH: the
marcher's instants stand as a chain — one vertex per instant, one edge
per step, walked in order."* That is production's own description of
itself, read but not imported.

Nothing here says the production marcher is wrong, or that
`graph_operation` must change. Those remain questions for levels that
have not been built.

## Seam A2 at seal — the count, and its limit

```text
BEFORE this tower     Derivative Action, Adjoint
                      = 2 towers, one family (derivative)

AFTER this tower      + Time Integration
                      = 3 independent towers, and the FIRST
                        non-derivative-family client
```

**Time contributes ONE tower vote**, not three. The tower produced
three separate REDs — `class_graph_step` at Level 6,
`class_graph_marcher` at Level 8, `class_graph_linearization` at Level 8
— but its ten levels are one client, and counting them separately would
be counting one experiment three times.

What *did* increase is **evidence quality**:

```text
3 independent towers                        (the strong-evidence bar)
+ first non-derivative-family client
+ first full temporal PRODUCTION COMPOSITION
      marcher → newton → difference_linearization → gmres
+ a genuine CLASS-2 WITNESS, reached naturally: no level of this
      tower may import class_graph_linearization, and the import
      gate enforces it
```

What that earned is **three narrow sentences of production**, each
RED-first:

> temporal discretization takes its state domain from the action it
> discretizes; a march takes its state domain from the action it
> marches; a finite difference takes its domain from the operation it
> differences.

It does **not** authorize the broader A2 migration across `fit`,
`graph_reduction`, `graph_broadcast`, or a rectangular linearization.
None of those was exercised by this tower.

**Recommended next action: a dedicated reverse architecture review.**
It has not been performed here, and
`doc/REVERSE-ARCHITECTURE-REVIEW.md` was deliberately left unedited —
that document is a separate artifact with its own process.

\[
\boxed{\begin{gathered}
\text{Strong evidence makes the broader refactor ELIGIBLE for reverse review;}\\
\text{it does not make every possible implementation MANDATORY.}
\end{gathered}}
\]

That question belongs to a reverse review after this tower is sealed —
and this tower is not sealed.

## Seams A1 and B — unmoved

```text
A1   graph host cannot be removed generically
     CLOSED by the Partitioned Implicit PDE Tower, on production
     evidence.  This client locally ignores H_context's topology
     because a triangular 2x2 decay HAS no topology to traverse -
     a property of THIS action, not a counterexample.  The host is
     still passed to every action, at every level, all the way
     through the march.  A1 stays closed.

B    one law, forward and reverse, between different domains
     ZERO new votes, and this needs saying carefully because the
     tower did touch difference_linearization.  It remains a
     SAME-DOMAIN citizen, L : Q -> Q.  Its forward use inside
     Newton is not bidirectional-linearization evidence: no
     tangent, no adjoint, no transpose, no rectangular U -> Y, and
     no reverse action appears anywhere in this tower.  Seam B
     still stands at 2 independent derivative-family towers.

A3   the (S, P) view as structural owner
     Another successful ownership pattern - G_time owns two axes
     and four relations, including a carrier no relation names.
     No production change follows automatically.  KEEP.
```

> Read carefully: **seam A2** is the review's *"operations should carry
> their own domain"*, and \(A_2\) is this tower's *two-step reach
> relation*. They share a label and nothing else. This document uses
> "seam A2" for the first and \(A_2\) for the second, everywhere.

## The specialization frontier

Recorded, deliberately **not** fixed:

| What production does | What the tower declared | Status here |
|---|---|---|
| regenerates a linear chain from `nsteps` | \(G_{\text{time}}\), a relational structure with \(A_1, A_2\) | **exact specialization** — the specimen's time graph *is* a simple chain, checked extensionally at L8 |
| carries one scalar `clock % step` | \(h : E \to \mathbb{R}\), one value per step | **exact specialization** — \(h\) is uniform, and L8 checks \(h(e)=\text{step}\) at every step |

Neither is a defect, and this tower establishes neither as a general
contract. They are the natural experiments for an **Adaptive Time** or
**Composite Time** client, which would supply a nonuniform \(h\) or a
branching time structure and find out whether the specialization holds.
Not implemented here.

---

# F. What this tower proves / does not prove

## Proven, through Gate B

```text
STRUCTURE (L0-L4)
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

VALUES AND MATHEMATICS (L5-L7)
state values live on Q, and a field needs no graph to exist
values are not structure: time(head(e)) - time(tail(e)) = h(e) is
    proved over the already-earned relations
an action may carry its own domain and ignore the host it is handed
temporal discretization PRESERVES the domain of the action it
    discretizes - after a narrow correction earned by RED
A1/A2 supply the history ROLES; bdf-2 supplies the COEFFICIENTS,
    and neither contains the other
production GMRES solves an implicit temporal step on Q while the
    host carries five unrelated vertices

CONSTITUTION (L8-L9)
production's control chain realizes A1 EXTENSIONALLY - two
    realizations of one structure, never one object
the marcher's state domain is independent of its host, and a
    four-step march runs on a two-member Q beside a five-vertex host
FE, BE and BDF2 all march to their exact rational trajectories,
    every prefix pinned, worst error 2.2e-16
a uniform h field and a linear time chain are exact
    SPECIALIZATIONS of what production carries - not defects
the complete initial-value problem is asked and answered on Q:
    q(t4) = [7/24, 83/144]
```

## Not proven — anywhere yet

```text
NONUNIFORM TIME              a variable step size that actually varies
NONLINEAR TIME STRUCTURE     branching, adaptive or composite time
    - production regenerates a linear chain from nsteps and carries
      one scalar step.  Both are EXACT here.  Whether they are
      general contracts is what an Adaptive Time or Composite Time
      client would decide, and this tower cannot

any tangent or adjoint       any transpose, any reverse action
any rectangular U -> Y linearization   ← SEAM B, still zero votes
    - difference_linearization remains SAME-DOMAIN, Q -> Q

whether fit / reduction / broadcast should carry their own domains
    ← the BROADER A2 migration, NOT performed and NOT tested here

a transient adjoint          more than one specimen
    - one 2x2 triangular decay, on one uniform chain
```

**Seam A2** of
[the reverse architecture review](../../doc/REVERSE-ARCHITECTURE-REVIEW.md)
— *operations should carry their domain rather than ask a graph for one*
— reaches **three independent towers** at Level 6, and this is the first
outside the derivative family. What that earned is one narrow sentence
of production, stated above. The migration across the other four call
sites is eligible for reverse review and has not been done here.

---

# G. Code map

```text
test/time-integration-tower/
├── README.md                    this document — the Rosetta stone
├── NUCLEUS-OBSERVATIONS.md      the evidence ledger (TI-*), by level
├── run.sh                       level-by-level runner; seals at Gate C
├── check_imports.sh             fail-closed allowlists, PER LEVEL,
│                                + its own --selftest
├── check_marker.sh              the result contract + its self-test
├── common/
│   ├── time_assert.f90                    (below everything)
│   ├── time_carriers_fixture.f90          earned at Level 0
│   ├── time_relations_fixture.f90         earned at Level 1
│   ├── time_algebra_fixture.f90           earned at Level 2
│   ├── time_fields_fixture.f90            earned at Level 5
│   └── triangular_decay_fixture.f90       earned at Level 6
├── level-0-carrier/             test.f90
├── level-1-relation/            test.f90
├── level-2-relation-algebra/    test.f90
├── level-3-graph/               test.f90
├── level-4-graph-calculus/      test.f90
├── level-5-field-calculus/      test.f90
├── level-6-discretization/      test.f90
├── level-7-minimization/        test.f90
├── level-8-constitution/        test.f90
└── level-9-statement/           test.f90
```

The fixture ladder is the tower's own stratification applied to itself:

```text
Level 0    time_carriers_fixture      declares Q, T, E
Level 1    time_relations_fixture     states Tail, Head over them
Level 2    time_algebra_fixture       composes what follows
Level 5    time_fields_fixture        puts values on all three
Level 6    triangular_decay_fixture   the action S : Q -> Q
```

The import ceiling rises the same way, one rung at a time — the step
operator refused below Level 6, the minimizer and GMRES below Level 7,
the marcher and Newton below Level 8 — so no level can redescribe
machinery it has not yet earned.

`class_graph_linearization` is refused at **every** level, and that
refusal is load-bearing rather than hygienic: Level 8 reached its
Class-2 defect through `marcher → newton → difference_linearization`,
and because no level may name that module, the failure cannot have been
manufactured.

The relation fixture does not *import* the carrier fixture — its
constructors receive \(Q,T,E\) as arguments, because a Level-1 file may
state facts over sets but may not name a set into existence. The ladder
is enforced by the import gate's per-file allowlists, by each level's
Makefile, and by `check_imports.sh --selftest`, which asserts that a
Level-0 source saying `use time_relations_fixture` is refused.

`check_marker.sh` arrived at Level 9, with the statement it validates.
It checks **shape and syntax only** — one marker, two real literals, in
\(Q\)'s declaration order — and says nothing about the values. Two
tokens, not five: the answer lives on \(Q\), and a five-token marker
would undo nine levels of argument at the last line of output.
