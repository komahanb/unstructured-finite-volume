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

# GATE B — Topology-consuming action *(unbuilt)*

Will ask whether applying a real Class-1 operation to overlap-complete
part fields, then assembling owned outputs, equals applying it globally —
and will put that operation behind production GMRES on the global graph.

# GATE C — Partitioned implicit statement *(unbuilt)*

Will ask whether the partitioned action can itself serve as the
matrix-free action of an implicit global solve.

---

# J. Graph roles at radius 2

This tower deliberately exercises the **legacy ordinary-graph / HPC
branch** of the nucleus. It introduces no `relational_graph`, because
its mathematics does not need one — and forcing one in merely because
other towers used one would be exactly the aesthetic reasoning the
reverse review forbids.

| Object | Type | Role |
|---|---|---|
| **global G** | `stored_graph` | whole topology; global field domain; target of assembly; (Gate B) GMRES contextual host |
| **parts G1, G2** | `stored_graph` with a part relation | **topology-consuming numerical operands**; local overlap carriers; global↔local maps; ownership environment |
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
| 9 statement | solve \(Aq=b\) by partitioned matrix-free action | Gate C |

---

# L. What this tower proves / does not prove

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
└── gate-a-partition/             test.f90
```

Development is grouped into three **gates**, not ten rungs: no empty
per-level directories exist, and the Rosetta table above maps each
gate's truths back onto the nucleus levels it consumes.

---

# Status

```text
partitioned implicit pde tower
├── Gate A · partition / ownership / transport ... PASS
├── Gate B · topology-consuming action .......... UNBUILT
└── Gate C · partitioned implicit statement ..... UNBUILT
```
