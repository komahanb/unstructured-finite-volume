# Adjoint Sensitivity Tower

The second **derivative-family** orbital client of the core-math nucleus,
after Derivative Action. One implicit problem is carried through the
whole tower.

Derivative Action established \(Jv\) and \(J^T\bar y\) from one
computation structure. It deliberately never *solved* anything. The
founding question here is different:

> Can the same nucleus perform a **primal solve**, solve the
> **transposed linearized problem**, and combine the result into a
> **total sensitivity** — without constructing an independent adjoint
> model?

The central law this tower must prove is one computational story:

\[
\boxed{R(q,p)=0},
\qquad
\boxed{R_q^{T}\lambda=f_q^{T}},
\qquad
\boxed{\frac{df}{dp}=f_p-\lambda^{T}R_p}.
\]

---

# 1. The complete mathematical statement

A reader should understand the whole problem before seeing any Fortran.

**Parameter** \(p=2\). **State** \(q=\begin{bmatrix}u\\v\end{bmatrix}\).

**Residual**

\[
R(q,p)=
\begin{bmatrix}
2u+v-4p\\
3u+4v-11p
\end{bmatrix}=0,
\qquad
\textbf{Response}\quad
f(q,p)=u+2v+2p .
\]

Therefore

\[
R_q=A=\begin{bmatrix}2&1\\3&4\end{bmatrix},
\qquad
R_p=\begin{bmatrix}-4\\-11\end{bmatrix},
\qquad
f_q=\begin{bmatrix}1&2\end{bmatrix},
\qquad
f_p=2 .
\]

The solution road, end to end:

\[
p=2
\;\Longrightarrow\;
R(q,p)=0
\;\Longrightarrow\;
\boxed{q=\begin{bmatrix}2\\4\end{bmatrix}}
\;\Longrightarrow\;
\boxed{f=14}
\]

\[
A^{T}\lambda=\begin{bmatrix}1\\2\end{bmatrix}
\;\Longrightarrow\;
\boxed{\lambda=\begin{bmatrix}-0.4\\0.6\end{bmatrix}}
\;\Longrightarrow\;
\frac{df}{dp}=f_p-\lambda^{T}R_p=2-(-5)
=\boxed{7}
\]

The **primary tower result** is \(df/dp=7\). Secondary truths:

```text
p                 = 2
q                 = [2, 4]
f                 = 14
lambda            = [-0.4, 0.6]
dq/dp             = [1, 2]
tangent  df/dp    = 7
adjoint  df/dp    = 7
```

## Why this specimen

It was chosen because it exercises, simultaneously:

```text
different state and residual domains
non-symmetric Jacobian                 A ≠ A^T
primal solve
transpose solve
parameter dependence inside the residual        (R_p)
direct parameter dependence in the response     (f_p)
indirect parameter dependence through the state (f_q q_p)
tangent total sensitivity
adjoint total sensitivity
```

Two choices are **load-bearing**:

- **\(A\) is not symmetric.** An implementation that quietly solves the
  primal system where it should solve the transpose *must fail*. With
  \(A^{T}\) the adjoint is \([-0.4,0.6]\); with \(A\) it would be
  \([-0.6,0.8]\) — different numbers, so the mistake cannot hide.
- **\(f_p=2\neq 0\).** The tower must not collapse total sensitivity to
  \(-\lambda^{T}R_p\) alone. That term is \(5\); the answer is \(7\).

---

# 2. Symbol dictionary

| Mathematics | Framework role | Domain | Orientation |
|---|---|---|---|
| \(p\) | parameter field | \(P\) | — |
| \(q=(u,v)\) | state field | \(Q\) | — |
| \(R=(r_1,r_2)\) | residual field | \(Y\) | — |
| \(f\) | response | \(Z\) | — |
| \(R_q\) | state linearization | — | \(Q\to Y\) |
| \(R_q^{T}\) | reverse / transpose action | — | \(Y\to Q\) |
| \(R_p\) | parameter action | — | \(P\to Y\) |
| \(f_q\) | response state derivative | — | \(Q\to Z\) |
| \(f_p\) | direct response derivative | — | \(P\to Z\) |
| \(\lambda\) | adjoint field | \(Y\) | — |

The adjoint \(\lambda\) lives on \(Y\): it is a covector on the residual
domain, and \(R_q^{T}\) carries it to \(Q\).

## Carriers and roles

```text
V = {p, u, v}          the variable carrier
T = {r1, r2, f}        the target carrier

P = {p}     ↪ V        parameter
Q = {u, v}  ↪ V        state
Y = {r1,r2} ↪ T        residual rows
Z = {f}     ↪ T        response
```

Roles are **domains**, not types: there is no `parameter_field`,
`state_field`, `residual_field`, or `adjoint_field` anywhere.

Note the deliberate index offsets: \(v\) is member 3 of \(V\) but entry
2 of \(Q\); \(f\) is member 3 of \(T\) but entry 1 of \(Z\). A raw
member id is never a storage position — every read goes through
`local_index`.

---

# 3. Domain orientation

```text
        P = {p}                     Q = {u,v}
           │                            │
           │ Rp                         │ Rq
           ▼                            ▼
        Y = {r1,r2}  ◀───────────────  Y = {r1,r2}
           │
           │ Rq^T
           ▼
        Q = {u,v}

        Q --fq--> Z = {f}   <--fp--  P = {p}
```

The two laws that make this tower harder than its predecessors:

\[
\boxed{R_q:Q\to Y}
\qquad\text{and}\qquad
\boxed{R_q^{T}:Y\to Q}
\]

and, stated as sharply as possible:

> **Equal dimensions do not imply equal domains.**
> \(|Q|=|Y|=2\), and \(Q\) is *not* \(Y\). The minimizer may demand
> equal value dimensions; it must never demand
> `Q same_as Y`.

**Relation orientation vs operator direction.** A support relation is
written (row, column):

```text
J_Q  ⊆ Y × Q     is the support of the operator  Rq   : Q → Y
J_Q^T ⊆ Q × Y    is the support of the operator  Rq^T : Y → Q
```

so the relation's first slot is the operator's *codomain*. Reading a
support relation as if it were the operator's direction is the single
easiest way to get this tower wrong.

---

# 4. Tower Rosetta table

| Level | Meaning in *this* problem | Framework abstraction | Source module | Test | New truth | Production |
|---|---|---|---|---|---|---|
| 0 | the four roles exist and differ | `counted_set`, `subset_set` | `graph_carrier` | `level-0-carrier/` | \(Q\neq Y\) though both size 2 | NONE |
| 1 | who may participate in what | `stored_relation` | `graph_relation` | `level-1-relation/` | \(R_{\mathrm{dep}}\subseteq T\times V\), nine facts | NONE |
| 2 | the four derivative supports | `restrict/project/compose`, `inclusion_of` | `graph_relation_algebra`, `graph_binary_relation` | `level-2-relation-algebra/` | \(J_Q,J_P,F_Q,F_P\) **derived** | NONE |
| 3 | one owned structural model | `relational_graph` | `graph_structure` | `level-3-graph/` | six carriers, five relations, closed | NONE |
| 4 | the implicit system is **cyclic** | `directed_adjacency_view`, `topological_order` | `graph_profile`, `graph_algorithms` | `level-4-graph-calculus/` | a topological order **refuses** | NONE |
| 5 | the first numbers | `field` | `class_graph_field` | `level-5-field-calculus/` | \(p=2\), \(q_0=[0,0]\); nothing fabricated | NONE |
| 6 | operator supports and orientation | `transpose_of` | `graph_binary_relation` | `level-6-discretization/` | \(J_Q^{T}\) swaps domain **identities** | NONE |
| 7 | primal and adjoint **solves** | `gmres` / minimizer | `class_graph_gmres` | `level-7-minimization/` | *Gate B — unbuilt* | — |
| 8 | one constitution, both actions | test-local law | — | `level-8-constitution/` | *Gate B — unbuilt* | — |
| 9 | one question: \(df/dp\) | composition | — | `level-9-statement/` | *Gate C — unbuilt* | — |

---

# 5. When structure, numbers and meaning appear

```text
0–4   structure and graph interpretation only — no number exists
5     numerical fields appear: p = 2, q0 = [0,0]
6     derivative/operator SUPPORT is identified — still no coefficient
7     equations are solved, with the laws still supplied opaquely
8     ONE constitution generates primal, tangent and transpose actions
9     ONE statement asks for df/dp
```

Gate A (levels 0–6) contains no \(-0.4\), no \(0.6\), no \(7\), and no
coefficient matrix at all.

---

# 6. Primal, tangent, adjoint

```text
PRIMAL                TANGENT                     ADJOINT
  R(q,p) = 0            Rq q_p = -Rp                Rq^T λ = fq^T
      ↓                     ↓                            ↓
      q                    q_p                           λ
                            ↓                            ↓
                     fq q_p + fp                  fp - λ^T Rp
```

\[
\text{tangent sensitivity}=7
\qquad=\qquad
\text{adjoint sensitivity}=7
\]

The two routes must meet on the same truth **independently** — neither
may call the other.

---

# 7. One source of truth

```text
structural dependency  (R_dep, and J_Q derived from it)
              ↓
constitution / local linearization
              ↓
             Rq
            /  \
          Jv    J^T λ
           |      |
       tangent  adjoint
```

Two duplications are forbidden by construction:

\[
\boxed{\text{adjoint support}=\text{transpose \emph{view}} \text{ of primal support}}
\]

\[
\boxed{\text{one local linearization, read in two directions}}
\]

No independently authored adjoint dependency relation. No separate
reverse coefficient table. No assembled global Jacobian as the
framework abstraction — the 2×2 matrix in this README is an *oracle for
the reader*, not a design.

---

# GATE A — Structural adjoint readiness

Gate A asks:

> Can the nucleus represent the complete structural distinction among
> parameter, state, residual and response **before** any numerical
> adjoint equation exists?

## Level 0 — Domains

Two parent carriers and four role subobjects:

```text
V = {p,u,v}      T = {r1,r2,f}
P = {p} ↪ V      Y = {r1,r2} ↪ T
Q = {u,v} ↪ V    Z = {f} ↪ T
```

**Verified:** cardinalities; both enumeration round trips on every
carrier; outsider rejection; all four embeddings
(`is_subobject_of`); and the identity separations

```text
¬(P same_as Q)     ¬(Y same_as Z)     ¬(V same_as T)
¬(Q same_as Y)     ← both have size 2, and they are still not the same
```

**Negative truth:** no relation, no graph, no field, no operator, no
adjoint, no coefficient.

## Level 1 — One structural dependency source

One base relation says who may participate:

\[
R_{\mathrm{dep}}\subseteq T\times V,
\qquad
|R_{\mathrm{dep}}|=9 .
\]

Every target depends on every variable in this specimen — \(r_1,r_2,f\)
each on \(p,u,v\). Handed ten tuples with one duplicate, the relation
holds nine: a relation is a set.

**Honest note.** This specimen is structurally *dense*: the interesting
structure is the **role partition**, not sparsity. Derivative Action's
message came from a sparsity pattern; this tower's comes from
orientation and identity. Recorded as observation AD-1.

**Negative truth:** \(R_{\mathrm{dep}}\) does not know \(2,1,3,4,-4,-11\)
or \(f_p=2\). Those are numerical law data and belong to Level 8.

## Level 2 — Derived role-specific supports

Nothing is hand-maintained. Using the subobjects' own relational faces
(`inclusion_of`) and the established composition convention
`compose_binary(R_AB, R_BC) = R_BC ∘ R_AB`:

```text
I_Y = inclusion_of(Y) ⊆ Y×T        I_Q = inclusion_of(Q) ⊆ Q×V
I_Z = inclusion_of(Z) ⊆ Z×T        I_P = inclusion_of(P) ⊆ P×V

J_Q = (I_Y ∘ R_dep) ∘ I_Q^T   ⊆ Y×Q
J_P = (I_Y ∘ R_dep) ∘ I_P^T   ⊆ Y×P
F_Q = (I_Z ∘ R_dep) ∘ I_Q^T   ⊆ Z×Q
F_P = (I_Z ∘ R_dep) ∘ I_P^T   ⊆ Z×P
```

with exact extensions

```text
J_Q :  (r1,u) (r1,v) (r2,u) (r2,v)        J_P :  (r1,p) (r2,p)
F_Q :  (f,u)  (f,v)                       F_P :  (f,p)
```

and every signature checked **by identity**, not by size. The two
right-hand inclusions are used transposed — as views, never rebuilt.

## Level 3 — Relational ownership

One `relational_graph` owns the whole structural model: six carriers
(\(V,T,P,Q,Y,Z\)) and five relations
(\(R_{\mathrm{dep}},J_Q,J_P,F_Q,F_P\)).

**Verified:** ownership by identity for every carrier and relation;
signature closure — every slot of every owned relation resolves to an
owned carrier, *including the subobject domains*; and graph identity
against an identically stocked twin.

> The complete adjoint problem has one structural ownership environment
> before it has an adjoint solution.

There is no `adjoint_graph`.

## Level 4 — The implicit system is not a program

This level does **not** fake an execution DAG. The state is implicitly
coupled, and the honest graph-calculus question is whether that coupling
is acyclic. From the same \(J_Q\):

\[
C_Q=J_Q^{T}\circ J_Q\subseteq Q\times Q,
\qquad
C_Q=\{(u,u),(u,v),(v,u),(v,v)\}
\]

— two slots sharing residual rows, so each depends on the other and on
itself. Interpreted as a directed graph, this view is perfectly valid
and it **has a cycle**:

```text
view construction .......... succeeds
reachable(u,v), reachable(v,u) ... both true
topological_order .......... REFUSES, loudly
```

\[
\boxed{\text{a valid directed graph}\neq\text{a DAG}}
\]

The refusal is the level's central truth, checked by message in
`check_refusals.sh`. An implicit system is solved by **minimization**
(Level 7), never by a topological walk. No acyclic fiction was invented
to make the walk succeed.

## Level 5 — The first numbers

```text
p  = [2]     field on P
q0 = [0, 0]  field on Q     — deliberately wrong; the solver's job
```

**Verified:** each field's domain is its declared `member_set` *by
identity*; storage follows domain enumeration and is read through
`local_index`.

**Negative truth:** no residual field, no response field, and above all
**no \(\lambda=0\)** fabricated merely because Level 7 will need an
adjoint. An initial guess belongs to the solve, not to Level 5.

## Level 6 — Operator supports and orientation

The derived relations are read as the discrete supports of

```text
Rq : Q → Y      Rp : P → Y      fq : Q → Z      fp : P → Z
```

and the adjoint support is the **transpose view of the same \(J_Q\)**:

\[
J_Q^{T}\subseteq Q\times Y,
\qquad
\{(u,r_1),(u,r_2),(v,r_1),(v,r_2)\}.
\]

**Verified:** the transpose swaps the actual **domain identities**
(`domain(1) same_as Q`, `domain(2) same_as Y`) — stronger than
transposing a 2×2 array; the view is *not* materialized while \(J_Q\)
is; and no second reverse support relation exists anywhere.

Gate A finishes with **no** primal solution, **no** adjoint solution,
**no** \(\lambda\), **no** coefficients, and **no** \(df/dp\).

---

# GATE B — Solve and constitution *(unbuilt)*

Will ask whether ordinary minimization can govern both \(R(q,p)=0\) and
\(R_q^{T}\lambda=f_q^{T}\), and whether one constituted linearization
generates the forward and reverse actions. Levels 7–8.

# GATE C — The statement *(unbuilt)*

Will ask one question — at \(p=2\), solve and report \(df/dp\) — and
answer it with the single real scalar \(7\). Level 9.

---

# What this tower proves / does not prove

## Proven by Gate A

```text
parameter, state, residual and response are distinguished by DOMAIN
equal cardinality does not collapse domain identity (Q vs Y)
one dependency relation generates all four derivative supports
subobject inclusions carry role structure into the algebra
one relational graph owns the complete structural model
the implicit state coupling is cyclic — and the nucleus says so by
    refusing a topological order rather than inventing one
the adjoint support is a transpose view, never a second relation
all of it with zero production changes
```

## Not proven — deliberately

```text
large PDE adjoints              transient / time-dependent adjoints
checkpointing                   distributed reverse communication
higher-order adjoints           multi-response batching
sparse Jacobian assembly performance                GPU performance
```

and, until Gates B and C are built: any solve, any \(\lambda\), any
sensitivity.

---

# Code map

```text
test/adjoint-tower/
├── README.md                     this document — the Rosetta stone
├── NUCLEUS-OBSERVATIONS.md       the evidence ledger
├── run.sh                        gate-grouped frontier-law runner
├── check_imports.sh              fail-closed per-level allowlists
├── common/
│   └── adjoint_assert.f90        constants + assertion helpers
├── level-0-carrier/              test.f90
├── level-1-relation/             test.f90
├── level-2-relation-algebra/     test.f90
├── level-3-graph/                test.f90
├── level-4-graph-calculus/       test.f90 · refusal.f90 · check_refusals.sh
├── level-5-field-calculus/       test.f90
└── level-6-discretization/       test.f90
```

`common/adjoint_assert.f90` is dependency-free and carries only the
symbolic member names (`VAR_P`, `VAR_U`, `VAR_V`, `TGT_R1`, `TGT_R2`,
`TGT_F`) — no coefficient, no numerical value. It is deliberately
neither `calculator_assert`, `learning_assert`, nor
`derivative_assert`: this is an independent fourth client.

## Gate mechanisms

- **`run.sh`** — the frontier law, grouped by architectural gate: a
  failed level blocks dependent work; the first absent level closes the
  frontier. Gates B and C report `UNBUILT` until they exist.
- **`check_imports.sh`** — fail-closed allowlists per level. Gate
  grouping does not weaken level-by-level dependency checking:
  Gate A may not import minimization or a solver, and levels 0–4 may not
  import the field.

---

# Architectural context

This tower is an orbital client of the core-math nucleus — see
[Fractal Graph Architecture](../../FRACTAL-GRAPH-ARCHITECTURE.md) — the
**second** of the derivative orbital family, at higher orbital energy
than Derivative Action because it must compose structure, constitution
*and* solves. Observations feed
[`NUCLEUS-OBSERVATIONS.md`](NUCLEUS-OBSERVATIONS.md) and are not
automatically promoted into production abstractions.

The sealed calculator, learning and derivative-action towers remain
acceptance oracles: this tower may expose new needs, but it may not
weaken their truths.

---

# Status

```text
adjoint sensitivity tower
├── Gate A · structure
│   ├── 0 carrier ................ PASS
│   ├── 1 relation ............... PASS
│   ├── 2 relation algebra ....... PASS
│   ├── 3 relational graph ....... PASS
│   ├── 4 graph calculus ......... PASS
│   ├── 5 field calculus ......... PASS
│   └── 6 discretization ......... PASS
├── Gate B · solve + constitution . UNBUILT
└── Gate C · statement ............ UNBUILT
```

The final claim this tower is being built to support is **not** "the
framework has an adjoint class". It is:

\[
\boxed{
\begin{gathered}
\text{the same mathematical structure supports primal action,}\\
\text{transpose action, primal solve, adjoint solve and total sensitivity}
\end{gathered}
}
\]

with no second adjoint truth maintained beside the primal one.
