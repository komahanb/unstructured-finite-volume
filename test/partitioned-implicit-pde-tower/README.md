# Partitioned Implicit PDE Tower

The first deliberate **radius-2** orbital client of the core-math
nucleus, after Calculator, Learning, Derivative Action and Adjoint
Sensitivity — and the experiment named by
[the reverse architecture review](../../doc/REVERSE-ARCHITECTURE-REVIEW.md)
as the one that would discriminate.

That review found the strongest-looking seam and rejected it: three
towers reported the minimizer's graph host as scenery, but production
showed the host is the **conduit** by which a topology-consuming
operation receives its mesh. All four towers had attached operations
with no topology to traverse, so none could see it.

> **This tower puts a real topology-consuming operation behind the
> minimizer, on a partitioned graph.** It is deliberately *not* a
> derivative-family client.

---

# A. The mathematical specimen

One problem carried through the whole tower. On the six-vertex chain
\(G\), with \(L\) the production vertex Laplacian:

\[
\boxed{A\,q = b},
\qquad
A = 2I - L .
\]

The exact state and the right-hand side are chosen so that both are
non-uniform:

\[
q^{*}=\begin{bmatrix}1&2&4&7&11&16\end{bmatrix}^{T},
\qquad
L\,q^{*}=\begin{bmatrix}1&1&1&1&1&-5\end{bmatrix}^{T},
\]

\[
b=(2I-L)\,q^{*}
 =\begin{bmatrix}1&3&7&13&21&37\end{bmatrix}^{T}.
\]

These are **test oracles**. \(A\) is never built as a \(6\times6\)
matrix: the execution road is

```text
q ──┬── 2q
    └── L(q)   ← the production Laplacian traverses the graph
         ↓
      2q − L(q)
```

The matrix appears in this README because it is the clearest oracle for
a reader; it is not a design.

## Honest naming

This is a **shifted graph-diffusion (shifted-Laplacian) problem**, not a
complete finite-volume PDE discretization. Its continuous analogue

\[
2q - q'' = b
\]

is intuition only. The framework truth is the discrete graph equation.

## Why this specimen

```text
the topology is genuinely traversed        (a real Class-1 operation)
the cut lies between vertices 3 and 4
owned vertex 3's stencil needs borrowed q(4)
owned vertex 4's stencil needs borrowed q(3)
the exact solution is non-uniform          (no accidental symmetry)
interface values cannot vanish without changing the answer
2I − L is non-singular                     (an implicit solve is posed)
no boundary-condition machinery needs inventing
```

---

# B. Topology

```text
GLOBAL G

  1 --e1--> 2 --e2--> 3 --e3--> 4 --e4--> 5 --e5--> 6

  tails = [1,2,3,4,5]        heads = [2,3,4,5,6]
  6 vertices, 5 edges, uniform spacing = measure = coefficient = 1
```

---

# C. The partition

Two linear parts, `partitioner(PARTITION_LINEAR, nparts=2, part=k)`:

```text
GLOBAL

  1 -- 2 -- 3 -- 4 -- 5 -- 6
              |
             cut

PART 1                          PART 2

  1 -- 2 -- 3 -- (4)           (3) -- 4 -- 5 -- 6
                 borrowed       borrowed
```

**The parentheses mean BORROWED VISIBILITY, not ownership.** A borrowed
vertex is present in the part because a stencil needs to *see* it; it
belongs to the other part.

## Vertex map — all six members

| global | global owner | G1 local | G1 status | G2 local | G2 status |
|---:|---:|---:|---|---:|---|
| 1 | part 1 | 1 | owned | — | absent |
| 2 | part 1 | 2 | owned | — | absent |
| 3 | part 1 | 3 | owned | 4 | **borrowed** |
| 4 | part 2 | 4 | **borrowed** | 1 | owned |
| 5 | part 2 | — | absent | 2 | owned |
| 6 | part 2 | — | absent | 3 | owned |

```text
G1 local order:  [1, 2, 3, 4]      owned owned owned BORROWED
G2 local order:  [4, 5, 6, 3]      owned owned owned BORROWED
```

Note G2's local order is deliberately **not** global order — local
member 1 is global vertex 4. Every value in this tower is therefore read
through `global_vertex_index` and `local_index`, never by position.

## Edge map — all five members

| global | tail→head | global owner | G1 local | G1 status | G2 local | G2 status |
|---:|---|---:|---:|---|---:|---|
| e1 | 1→2 | part 1 | 1 | owned | — | absent |
| e2 | 2→3 | part 1 | 2 | owned | — | absent |
| **e3** | **3→4** | **part 1** | 3 | **owned** | 4 | **borrowed** |
| e4 | 4→5 | part 2 | — | absent | 1 | owned |
| e5 | 5→6 | part 2 | — | absent | 2 | owned |

The **crossing edge e3 is present in both parts** and owned by exactly
one — see §9 of the Gate-A section for how that ownership was *derived*
rather than assumed.

---

# D. Ownership vocabulary

```text
OWNED      this part answers for it; its contribution is assembled
BORROWED   this part can SEE it; its contribution is discarded
OVERLAP    owned ∪ borrowed — everything locally present
ABSENT     not in this part at all
```

## The central partition law

The part carrier is **not** merely the owned subset. A
topology-consuming stencil needs its neighbours:

\[
\boxed{V_{\text{part}} = \text{owned}\;\cup\;\text{borrowed}}
\]

is the numerical **input** carrier. But only owned members contribute
back:

\[
\boxed{
\begin{gathered}
\text{local evaluation domain} = \text{overlap}\\
\text{assembly contribution domain} = \text{owned}
\end{gathered}}
\]

Stated as the law this tower keeps repeating:

```text
borrowed INPUTS  are necessary
borrowed OUTPUTS are disposable copies
```

Conflating the two is the single easiest way to get a partitioned
operator wrong — either by starving a stencil, or by double-counting a
contribution.

---

# GATE A — Partition, ownership, transport

Gate A asks:

> Does the existing graph partition machinery preserve enough structural
> and field truth for a topology-consuming numerical operation to run
> locally?

No operator. No solver. Structure, ownership, visibility, transport and
reconstruction only.

## What is pinned

**Global topology** — six vertices, five edges, each edge's tail and
head, and the stable identities of `vertex_set()` and `edge_set()`.

**Part structure**, for each part: `has_part_relation`, `num_parts = 2`,
`id`, and the full `global_vertex_index` / `vertex_owner_part` maps,
against `owned_vertices`, `borrowed_vertices` and `overlap_vertices`.

**The crossing-edge law, derived not guessed.** The invariant imposed
first is

\[
\boxed{\text{one global entity}\;\longrightarrow\;\text{one assembled
contribution}}
\]

An edge probe field \(z=[10,20,30,40,50]\) is partitioned to both parts,
each part's field is assembled home, and the two contributions summed.
The result must be exactly \(z\) — in particular the crossing edge e3
must contribute \(30\) **once**. Only after that law is pinned does the
tower inspect `edge_owner_part` to *document* which part is canonical.

**Full vertex transport.** \(q^{*}\) partitioned to both parts becomes a
full **overlap** field on each — `q1.domain same_as G1.vertex_set()`,
likewise G2 — with values checked by global member:

```text
G1 globals [1,2,3,4]  →  q1 = [1, 2, 4, 7]
G2 globals [4,5,6,3]  →  q2 = [7, 11, 16, 4]     ← not global order
```

**Round trips.** Assembling each part's field home yields only that
part's *owned* contribution; the two sum to \(q^{*}\) exactly, with no
borrowed value counted twice.

**Proper-subset transport.** A global vertex subset
\(S=\{6,3,4\}\hookrightarrow V_G\), declared in non-global order with
unmistakable values \(600,300,400\), is partitioned, read through
`local_index`, and assembled home to recover \(S\) extensionally. This
carries the learning tower's Phase-5B subdomain law out to radius 2.
Transformed subsets are *not* required to keep \(S\)'s identity token.

## Gate-A negative truths

```text
no laplacian     no PDE residual    no GMRES
no adjoint       no dense matrix    no MPI
no halo-exchange class
```

---

# GATE B — Topology-consuming action

Gate B asks **two** different questions and proves them apart:

> Does a real Class-1 operation, evaluated on overlap-complete part
> fields and assembled through **owned** outputs, reproduce the global
> operator?

> And does the graph a minimizer carries become observably load-bearing
> when the attached action genuinely traverses topology?

## The operation

`common/shifted_laplacian_fixture.f90` — test-local, and it **owns no
graph**. Its graph arrives through `domain()` and `apply()`, and it hands
that graph straight to production:

```text
caller / minimizer
    │  supplies a graph
    ▼
shifted_laplacian        A(q) = 2q − L(q)     ← no topology here
    │
    ▼
production laplacian()   ← traverses the graph it is handed
    │
    ▼
input_graph incidence
```

The adapter never inspects `edge_tail` or `edge_head`, implements no
incidence and reproduces no Laplacian loop. It also refuses a state of
the **right size** on a foreign carrier — six members that are not
`V(G)` — by domain identity (`check_refusals.sh`).

## Road 1 — the global action

```text
q* on V(G)
    ↓ shifted_laplacian(G)
    ↓ production L traverses G
L q* = [1, 1, 1, 1, 1, −5]
A q* = 2q* − Lq* = [1, 3, 7, 13, 21, 37] = b
```

## Road 2 — the partitioned action

```text
q* on V(G)
    ↓ partition_data
q1 on overlap V(G1)          q2 on overlap V(G2)
    ↓ shifted_laplacian(G1)      ↓ shifted_laplacian(G2)
A1 q1                        A2 q2
    ↓ assemble owned only        ↓ assemble owned only
            ╲                   ╱
             sum over parts
                  ↓
        [1, 3, 7, 13, 21, 37]   ← equals the global action, exactly
```

### The local answers, and which of them are authoritative

| part | global member | status | q | L q | A q | authoritative? |
|---|---:|---|---:|---:|---:|---|
| G1 | 1 | owned | 1 | 1 | 1 | **yes** |
| G1 | 2 | owned | 2 | 1 | 3 | **yes** |
| G1 | 3 | owned | 4 | 1 | 7 | **yes** |
| G1 | 4 | *borrowed* | 7 | −3 | **17** | **NO** — global says 13 |
| G2 | 4 | owned | 7 | 1 | 13 | **yes** |
| G2 | 5 | owned | 11 | 1 | 21 | **yes** |
| G2 | 6 | owned | 16 | −5 | 37 | **yes** |
| G2 | 3 | *borrowed* | 4 | 3 | **5** | **NO** — global says 7 |

A part holds enough topology to answer for what it **owns** and no more:
the borrowed vertex sits at the part's edge with half its stencil
missing. Those two rows are why owned assembly is not an optimisation
but a correctness requirement.

\[
oxed{A_G\,q \;=\; \sum_p P_p^{T}\,O_p\,A_{G_p}\,P_p\,q}
\]

where \(P_p\) exposes the overlap state, \(O_p\) keeps only owned
outputs, and \(P_p^{T}\) restores global numbering. These are
*mathematics*, not production types — nothing named `P` or `O` exists.

## Borrowed input is numerically load-bearing

Structural visibility (Gate A) is not the same claim as numerical
dependence. Proved by perturbation, with the borrowed seat located
through `global_vertex_index` rather than assumed:

```text
G1: q(borrowed global 4) += 10  →  owned A at global 3:  7 → −3
G2: q(borrowed global 3) += 10  →  owned A at global 4: 13 →  3
restore the halo → the correct answers return
```

\[
oxed{	ext{borrowed INPUT}\;\longrightarrow\;	ext{owned OUTPUT}}
\]

## Road 3 — the solver-host conduit

The baseline implicit solve puts the Class-1 operation behind production
GMRES on the global graph:

```text
GMRES
  │ carries G                       (it reads no topology itself)
  ▼
attached shifted_laplacian
  ▼
production laplacian
  ▼
G topology
```

`attach(shifted, G, V(G))`; the affine constant is zero (A is linear);
`solver.apply(G, [b], sol)` returns a field on `V(G)` equal to
\(q^{*}=[1,2,4,7,11,16]\).

**And the conduit is proved behaviourally, not by reading code.** The
same `shifted_laplacian` — which stores no graph — is attached to two
solvers over two topologies with the *same* six vertices and *same* five
edges:

```text
G      chain:  1→2→3→4→5→6
G_alt  star :  1→2, 1→3, 1→4, 1→5, 1→6

solver_G   % matvec(q*)  =  b            ✓
solver_alt % matvec(q*)  ≠  b            ← nothing changed but the host
```

If the host were scenery the two would agree. They do not.

> The finding is **not** "GMRES traverses the graph" — it does not. It
> is: **GMRES carries the graph to the attached operation, and that
> operation consumes the topology.** Two distinct roles:
> the minimizer holds the graph as *conduit / context carrier*; the
> differential operation reads it as *numerical topology operand*.

# GATE C — The partitioned implicit statement

Gate C asks whether

> partition + overlap refresh + local topology-consuming actions +
> owned assembly

compose into **one** global matrix-free operation that ordinary
production GMRES can solve with.

## Three things that must never be conflated

\[
\boxed{
\begin{array}{lll}
\textbf{STRUCTURAL PARTITION} & G\to\{G_1,G_2\} & \textbf{once} \\
\textbf{NUMERICAL OVERLAP REFRESH} & q\to\{q_1,q_2\} & \textbf{every matvec} \\
\textbf{OWNED ASSEMBLY} & \{A_1q_1,A_2q_2\}\to Aq & \textbf{every matvec}
\end{array}}
\]

All three are routinely called "partitioning". They are not the same
thing, and the composite's *shape* enforces the distinction: it owns
the structure and owns **no mutable numerical state at all** — no
cached `q`, no cached overlap, no previous result. Nothing in it can go
stale because there is nothing in it to go stale.

## The central road

```text
GLOBAL q
   │
   ├──── refresh G1 overlap ──> q1 ──> A_G1 ──> owned y1 ──┐
   │                                                        │
   └──── refresh G2 overlap ──> q2 ──> A_G2 ──> owned y2 ──┤
                                                            ↓
                                                    assemble + sum
                                                            ↓
                                                        GLOBAL Aq
                                                            ↓
                                                          GMRES
                                                            ↓
                                                           q*
```

\[
\boxed{
\textbf{VISIBILITY}\ \text{governs what a local calculation may READ.}\\
\textbf{OWNERSHIP}\ \text{governs what it may authoritatively WRITE back.}}
\]

## Two kinds of graph-dependence

| | `shifted_laplacian` (Gate B) | `partitioned_shifted_laplacian` (Gate C) |
|---|---|---|
| relationship to a graph | **graph-parameterized** | **decomposition-context-bound** |
| acts on | whatever graph it is handed | only the \(G\) its \(G_1,G_2\) were cut from |
| stores | nothing | \(G\), \(G_1\), \(G_2\), partitioners, assembler |
| a foreign six-vertex host | is a legitimate different problem | is **refused** |

Gate B's star convicted the host by *changing the answer*; here the same
star is refused outright — in `domain()`, so a solver attaching on it
dies before `attach` completes rather than deep inside a matvec with a
chain-derived decomposition. Same cardinality, wrong decomposition.

## What is proved

**Extensional equality on four probes**, by global member —
\(q^{*}\), a mixed-sign vector \([3,-1,4,1,5,-9]\), and the two
**interface basis vectors** \(e_3\) and \(e_4\), which are the ones
that force information across the cut:

\[
A_{\text{partitioned}}(v) = A_{\text{global}}(v).
\]

**No stale overlap.** One composite instance applied five times —
\(q^{*}\), mixed, \(e_3\), \(e_4\), \(q^{*}\) — returns to its
first answer (\(y_1=y_5\)), every intermediate answer matches the
global action *at the time it was asked*, and the interleaved probes
genuinely differ, so the sequence is not a repeated no-op. A halo
cached from a previous matvec would survive into the next and break
this.

**Equivalence at the seat GMRES consumes.** Not only between the
fixtures but at `solver % matvec` for all four probes, with both affine
constants zero.

**The statement.** Both solves run independently from the same
\(b=[1,3,7,13,21,37]\):

\[
q_{\text{partitioned}} = q_{\text{global}} = q^{*} = [1,2,4,7,11,16].
\]

The decomposition changed the road, not the answer.

## The result marker

```text
PARTITIONED_PDE_RESULT =  1.0000000000000002E+00  2.0000000000000009E+00
                          4.0000000000000009E+00  7.0000000000000027E+00
                          1.1000000000000002E+01  1.6000000000000004E+01
```

One marker, six real tokens in global vertex order, at full precision
and **unrounded** — the honest image of the arithmetic. `check_marker.sh`
validates shape and syntax only, never values; whether those numbers
*are* \(q^{*}\) is the Gate-C test's business.

## Why this is NOT a distributed solver

The road still has, and this tower says so plainly:

```text
one process                    one global trial vector
direct global access during partition_data
parts executed sequentially    global assembly in-process
global-array inner products    global-array norms
```

So the following are **not** claimed: distributed GMRES, MPI solver,
parallel halo exchange, distributed vector support. The correct names
are **partitioned matrix-free solve** and *serial semantic model of a
partitioned operator*. What would be needed to make it genuinely
distributed is derived — not implemented — in
[`NUCLEUS-OBSERVATIONS.md`](NUCLEUS-OBSERVATIONS.md) under PIP-8.

---

# J. Graph roles at radius 2

This tower deliberately exercises the **legacy ordinary-graph / HPC
branch** of the nucleus. It introduces no `relational_graph`, because
its mathematics does not need one — and forcing one in merely because
other towers used one would be exactly the aesthetic reasoning the
reverse review forbids.

| Object | Type | Role |
|---|---|---|
| **global G** | `stored_graph` | topology **yes**; global domain source **yes**; GMRES host **yes**; downstream numerical influence **yes** (proved by the star-host control) |
| **parts G1, G2** | `stored_graph` with a part relation | partition frame **yes**; local domain source **yes**; local topology operand **yes**; ownership/global maps **yes** — one object, four roles |
| **borrowed member** | a member of a part carrier | visible to the local operator **yes**; authoritative output contributor **NO** |
| **vertex/edge sets** | `counted_set` / `subset_set` | value domains — what fields live on |

Four distinct senses of "graph" are in play, and the tower keeps them
apart:

```text
graph as TOPOLOGY         what the Laplacian traverses          (Gate B)
graph as CONDUIT          what the minimizer carries so its
                          action can traverse something         (Gate B)
graph as PARTITION FRAME  what holds global↔local maps and
                          ownership                             (Gate A)
member_set as DOMAIN      what a field lives on; never a graph
```

---

# K. Nucleus Rosetta

This tower **consumes** most nucleus levels rather than reconstructing
them — and says so honestly rather than manufacturing per-level code:

| Level | What this tower uses it for | Reconstructed or consumed? |
|---|---|---|
| 0 carrier | global and part vertex/edge domains | consumed |
| 1 relation | incidence, embodied by the ordinary graph topology | consumed |
| 2 relation algebra | none required — no new derivation here | **not used** |
| 3 graph | global and part structural ownership | consumed |
| 4 graph calculus | partition / overlap / local topology interpretation | **exercised at a new radius** |
| 5 field calculus | global, overlap and transported fields | **exercised at a new radius** |
| 6 discretization | \(L\) and \(A=2I-L\) on the graph | Gate B |
| 7 minimization | GMRES | Gate B, C |
| 8 constitution | the shifted-diffusion law \(A(q)=2q-Lq\) | Gate B |
| 9 statement | solve \(Aq=b\); **implementation** production GMRES; **matvec realization** serial partitioned composite; **structural context** \(G\to G_1,G_2\) once; **state context** \(q\to q_1,q_2\) every matvec; **result** field \(q^{*}\) on \(V(G)\) | Gate C |

---

# L. What this tower proves / does not prove

## Proven by Gate C

```text
partition + overlap refresh + local actions + owned assembly compose
    into ONE global matrix-free operation
A_partitioned = A_global on four probes, including both interface
    basis vectors
the composite holds STRUCTURE and no mutable numerical state, so no
    halo can go stale - proved by five interleaved applications
    returning to the first answer
equivalence holds at solver % matvec, the seat GMRES consumes
production GMRES solves through the composite:
    q_partitioned = q_global = q*
a decomposition is bound to its graph: a same-sized foreign host is
    refused in domain(), before attach completes
```

## Proven by Gate B

```text
a real Class-1 operation, applied to overlap-complete part fields and
    assembled through owned outputs, reproduces the global operator
    EXACTLY
borrowed outputs are demonstrably NOT authoritative (17 vs 13, 5 vs 7)
borrowed INPUT numerically determines owned OUTPUT — proved by
    perturbation, not by the existence of a borrowed set
a state of the right SIZE on a foreign carrier is refused by identity
the implicit global solve returns q* through production GMRES
the graph a minimizer carries is LOAD-BEARING: same operation, same
    probe, different host topology, different mathematical action
```

## Proven by Gate A

```text
a linear partition of a chain yields the expected ownership
overlap contains exactly the borrowed neighbours a stencil needs
presence is not ownership: a member can be visible and unowned
global↔local maps are the only honest way to read part storage
a crossing edge is assembled exactly once
full vertex fields, full edge fields and proper subsets all survive
    partition and reassembly
```

## Not proven — anywhere in this tower

```text
MPI communication              asynchronous halo exchange
distributed vectors            distributed dot products
distributed Krylov orthogonalization
load balancing                 scalability
communication hiding           parallel performance
multi-rank failure handling
```

The phrase **"distributed solver" is not used** for what is a serial
semantic decomposition. The honest names are *partitioned matrix-free
solve* and *serial semantic model of a partitioned operator*.

---

# Code map

```text
test/partitioned-implicit-pde-tower/
├── README.md                     this document — the Rosetta stone
├── NUCLEUS-OBSERVATIONS.md       the evidence ledger (PIP-*)
├── run.sh                        gate-grouped runner
├── check_imports.sh              fail-closed per-gate allowlists
├── common/
│   └── partitioned_pde_assert.f90
│   └── shifted_laplacian_fixture.f90
│   └── partitioned_shifted_laplacian_fixture.f90
├── gate-a-partition/             test.f90
├── gate-b-operator/              test.f90 · refusal.f90
│                                 · check_refusals.sh
└── gate-c-statement/             test.f90 · refusal.f90
                                  · check_refusals.sh
```

`check_marker.sh` holds the result contract and self-tests before the
ladder runs.

The import gate keys its allowlists **per file** inside `common/`, so the
assert module's freedom from framework imports is proved mechanically
rather than asserted in a comment — the fixture's permissions do not
extend to it.

Development is grouped into three **gates**, not ten rungs: no empty
per-level directories exist, and the Rosetta table above maps each
gate's truths back onto the nucleus levels it consumes.

---

# Status

```text
partitioned implicit pde tower
├── Gate A · partition / ownership / transport ... PASS
├── Gate B · topology-consuming action .......... PASS
└── Gate C · partitioned implicit statement ..... PASS

solution field on V(G) ........................ q* = [1,2,4,7,11,16]
```

**The tower is complete and sealed**, with **zero production changes**
beyond a single corrected comment earned at Gate A.
