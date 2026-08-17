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

## Levels are the architecture; gates are only review checkpoints

The tower is **ten levels**, 0 through 9. Each level is a directory, a
test, an import ceiling and a distinct mathematical responsibility.

Review happens at three checkpoints — after Level 4, after Level 7,
after Level 9 — and that is *all* a gate is. Gates own no mathematics,
appear in no directory name, and never stand in for a level's own
result. In the runner they are horizontal separators between levels
that have each reported for themselves.

---

# A. The mathematical specimen

On the six-vertex chain \(G\), with \(L\) the production vertex
Laplacian:

\[
\boxed{A\,q = b},
\qquad
A = 2I - L .
\]

\[
q^{*}=\begin{bmatrix}1&2&4&7&11&16\end{bmatrix}^{T},
\quad
L\,q^{*}=\begin{bmatrix}1&1&1&1&1&-5\end{bmatrix}^{T},
\quad
b=\begin{bmatrix}1&3&7&13&21&37\end{bmatrix}^{T}.
\]

These are **oracles**. \(A\) is never built as a matrix: the road is
\(q \mapsto 2q - L(q)\), with the Laplacian traversing the graph.

**Honest naming.** A *shifted graph-diffusion (shifted-Laplacian)*
problem, not a complete finite-volume discretization. The continuous
analogue \(2q-q''=b\) is intuition only.

**Why this specimen.** The topology is genuinely traversed; the cut
falls between vertices 3 and 4; owned vertex 3's stencil needs borrowed
\(q(4)\) and owned vertex 4's needs borrowed \(q(3)\); the solution is
non-uniform; \(2I-L\) is non-singular; no boundary-condition machinery
needs inventing.

---

# B. The ten-level Rosetta table

| Level | Mathematical object | Framework object | Domains | Relations | Numerical meaning | Test | Truth established | Production consequence |
|---|---|---|---|---|---|---|---|---|
| **0** | \(V,E,K\) | `counted_set` | three carriers | — | none | `level-0-carrier/` | carriers precede structure; the integer 1 is a vertex, an edge *and* a part | none |
| **1** | \(\mathrm{Tail},\mathrm{Head}\subseteq E\times V\); \(\mathrm{Own}\subseteq K\times V\) | `csr_relation` | \(E\times V\), \(K\times V\) | primitive facts | none | `level-1-relation/` | the chain and its *intended* ownership, stated before any graph | none |
| **2** | \(A=\mathrm{Head}\circ\mathrm{Tail}^{T}\); \(\mathrm{TailOwner}=\mathrm{Own}^{T}\circ\mathrm{Tail}\); \(\mathrm{HeadOwner}=\mathrm{Own}^{T}\circ\mathrm{Head}\) | `transpose_of`, `compose_binary` | \(V\to V\), \(E\to K\) | **derived** | none | `level-2-relation-algebra/` | relation algebra exposes tail-based **and** head-based ownership policies; both are total functions, and the crossing edge distinguishes them | none |
| **3** | \(G\) | `stored_graph` | own carriers | realizes \(\mathrm{Tail},\mathrm{Head}\) | none | `level-3-graph/` | the ordinary graph realizes the relation structure — extensionally, not by identity | none |
| **4** | \(G\to\{G_1,G_2\}\) | `partitioner` | part carriers + global maps | `edge_owner_part` vs *both* Level-2 policies | none | `level-4-graph-calculus/` | owned / borrowed / overlap; **production partitioning realizes the tail-based policy** | none |
| **5** | fields under restriction and extension | `field`, `partition_data`, `assemble_data` | overlap and owned | — | values move | `level-5-field-calculus/` | **read = overlap, write-back = owned**; one entity, one contribution | none |
| **6** | \(L\) and \(A=2I-L\) | `laplacian`, `shifted_laplacian` | part carriers | incidence traversed | the discrete law | `level-6-discretization/` | local actions assemble to the global action; borrowed input is numerically load-bearing | none |
| **7** | solve \(Aq=b\) | `gmres` | \(V(G)\) | — | the implicit solve | `level-7-minimization/` | **the graph host is a conduit** — chain vs star changes the answer | none |
| **8** | \(A_{\text{part}}=A_{\text{global}}\) | `partitioned_shifted_laplacian` | \(V(G)\to V(G)\) | — | same map, alternate constitution | `level-8-constitution/` | structure once, overlap every apply, no cached state | none |
| **9** | the statement | GMRES + partitioned action | \(V(G)\) | — | \(q^{*}\) | `level-9-statement/` | \(q_{\text{part}}=q_{\text{global}}=q^{*}\) | none |

The road a reader can follow:

```text
L0  carriers
L1  primitive incidence + vertex ownership facts
L2  relational consequences AND candidate edge-ownership policies
L3  ordinary graph realization
L4  partition transform + empirical policy selection: production == TailOwner
L5  values transported and reconstructed under that selected policy
L6  the discrete operator over the local topology
L7  minimization, and what the graph host does
L8  the partitioned constitution
L9  the statement
```

---

# C. Topology, the cut, and the vocabulary

```text
GLOBAL G

  1 --e1--> 2 --e2--> 3 --e3--> 4 --e4--> 5 --e5--> 6
                       |
                      cut

PART 1                          PART 2

  1 -- 2 -- 3 -- (4)           (3) -- 4 -- 5 -- 6
                 borrowed       borrowed
```

**Parentheses mean BORROWED VISIBILITY, not ownership.**

| global | owner | G1 local | G1 status | G2 local | G2 status |
|---:|---:|---:|---|---:|---|
| 1 | part 1 | 1 | owned | — | absent |
| 2 | part 1 | 2 | owned | — | absent |
| 3 | part 1 | 3 | owned | 4 | **borrowed** |
| 4 | part 2 | 4 | **borrowed** | 1 | owned |
| 5 | part 2 | — | absent | 2 | owned |
| 6 | part 2 | — | absent | 3 | owned |

```text
G1 local order:  [1, 2, 3, 4]      G2 local order:  [4, 5, 6, 3]
```

G2's local order is deliberately **not** global order: local member 1 is
global vertex 4. Every value is read through `global_vertex_index` and
`local_index`.

| global edge | tail→head | owner | in G1 | in G2 |
|---|---|---:|---|---|
| e1 | 1→2 | part 1 | owned | absent |
| e2 | 2→3 | part 1 | owned | absent |
| **e3** | **3→4** | **part 1** | **owned** | **borrowed** |
| e4 | 4→5 | part 2 | absent | owned |
| e5 | 5→6 | part 2 | absent | owned |

The `owner` column is **production's**, and it follows the tail-based
policy — e3 leaves vertex 3, which part 1 owns. A head-based partitioner
would put e3 in part 2 and reconstruct just as exactly; Level 2 derives
both candidates and Level 4 establishes which one this partitioner
implements.

```text
OWNED      this part answers for it; its contribution is assembled
BORROWED   this part can SEE it; its contribution is discarded
OVERLAP    owned ∪ borrowed — everything locally present
```

\[
\boxed{
\begin{gathered}
\textbf{VISIBILITY}\ \text{governs what a local calculation may READ.}\\
\textbf{OWNERSHIP}\ \text{governs what it may authoritatively WRITE back.}
\end{gathered}}
\]

---

# Levels 0–2 — before there is a graph

**Level 0** declares \(V\) (6 vertices), \(E\) (5 edges) and \(K\) (2
part labels) and nothing else. Its sharpest truth is a hazard: all three
carriers enumerate from one, so **the integer 1 is a member of all
three**. Only identity separates them, and every level above depends on
that being true here.

**Level 1** states three primitive facts as binary relations:

```text
Tail ⊆ E×V     e_i → i          Head ⊆ E×V     e_i → i+1
Own  ⊆ K×V     part1 → 1,2,3    part2 → 4,5,6
```

Ownership here is an **intention over carriers** — no partitioner exists
to realize it. Every signature is pinned by identity, because with
colliding ids the *same* integer pair is a fact of two unrelated
relations:

```text
[1,1]   Tail : e1 leaves vertex 1        Own : part1 owns vertex 1
[1,2]   Head : e1 enters vertex 2        Own : part1 owns vertex 2
```

Same numerals, unrelated meanings; only the signature — \(E\times V\)
versus \(K\times V\) — tells them apart. **A raw integer tuple is not a
typed relational fact.** No pair belongs to all three relations and none
could: `[3,4]` is Head's alone, and `Own` refuses it because \(K\) has
no member 3 — a carrier truth, not an arithmetic one.

**Level 2** derives what follows, and this level is *inhabited*, not
N/A:

\[
A=\mathrm{Head}\circ\mathrm{Tail}^{T}:V\to V,
\qquad
\mathrm{TailOwner}=\mathrm{Own}^{T}\circ\mathrm{Tail}:E\to K,
\qquad
\mathrm{HeadOwner}=\mathrm{Own}^{T}\circ\mathrm{Head}:E\to K.
\]

The two ownership maps are the ones that matter, and deriving **both**
is the point. An edge has two vertices; each has an owner; composing
through either gives a perfectly good map \(E\to K\). Both are total,
both are single-valued, so *both* satisfy the reconstruction law — one
edge, one owner. On this chain they agree everywhere except at the
crossing edge:

```text
             e1     e2     e3     e4     e5
TailOwner  part1  part1  part1  part2  part2
HeadOwner  part1  part1  part2  part2  part2
                          ↑
                          the cut — the only place the choice shows
```

The distinction the level insists on:

\[
\boxed{
\begin{aligned}
\textbf{DERIVATION} &:\ \text{once a relational expression is chosen,
its extension follows algebraically.}\\
\textbf{POLICY SELECTION} &:\ \text{choosing Tail rather than Head as
the ownership anchor is an additional}\\ &\ \ \ \text{semantic decision.}
\end{aligned}}
\]

So \(\mathrm{TailOwner}:=\mathrm{Own}^{T}\circ\mathrm{Tail}\) is a
**definition — a policy expressed relationally**. *Given* that
definition, these are theorems: \(\mathrm{TailOwner}(e_3)=\text{part1}\),
and every edge has exactly one `TailOwner`. Likewise for `HeadOwner`.
What is *not* a theorem is that production must anchor at the tail.

The earlier gate-shaped tower found tail-ownership *operationally* — it
imposed an assembly law on a probe field and read back what the
partitioner had done. That was the right way to find it. What Level 2
adds is vocabulary: both candidates stated as mathematics **before any
partitioner exists**, so Level 4's reading of production becomes a
*selection between two named alternatives* rather than a bare
observation.

```text
        LEVEL 2                          LEVEL 4

  relational TailOwner  ─────┐
                             ├──?──  production edge_owner_part
  relational HeadOwner  ─────┘

  two candidates derived           one of them selected, empirically
```

---

# Level 3 — the ordinary graph realizes the structure

`stored_graph(6, tails=[1..5], heads=[2..6])`, checked against the
Level-1 oracle **extensionally and signature-aware**, never by carrier
identity: \(G\) builds its own carriers, and demanding they be the same
objects as \(V\) and \(E\) would be demanding that two parties who agree
must be the same party. Every `edge_tail` and `edge_head` must satisfy
`Tail` and `Head`, and nothing unlicensed may appear.

# Level 4 — the partition interpretation

Where the production partitioner first appears, and where *owned*,
*borrowed* and *overlap* first mean something. It pins both parts'
`global_vertex_index` maps, their owner functions, their
owned/borrowed/overlap subsets, and the crossing edge's presence in
**both** parts — then reads `edge_owner_part` against **both** Level-2
candidate policies, edge for edge.

The crossing edge is the entire discriminator, since it is the only
place the two policies differ:

```text
production(e3) = part1        TailOwner(e3) = part1        ✓ agrees
                              HeadOwner(e3) = part2        ✗ differs
```

\[
\boxed{\textbf{PRODUCTION SELECTS THE TAIL-BASED POLICY.}}
\]

Not *"production satisfies the uniquely forced ownership theorem"* —
there is no such theorem. The finding is empirical and specific, and it
names an old prose defect exactly: the partitioner's comment once
described a *hybrid* — tail owner normally, head owner when the tail is
borrowed — while the executable behaviour was uniformly `TailOwner`.
With both policies in hand, that sentence is no longer merely wrong; it
is wrong in a way this level can spell.

Graph-to-graph interpretation **only**: no field is transported here.

> **REVIEW GATE A** — after Level 4.

---

# Level 5 — field transport

A full global field becomes a **full field on each part's whole vertex
carrier**, borrowed member included:

```text
G1 globals [1,2,3,4]  →  q1 = [1, 2, 4, 7]
G2 globals [4,5,6,3]  →  q2 = [7, 11, 16, 4]     ← not global order
```

Assembling each part home contributes only what it **owns**, so the two
tile the whole exactly. The same law is checked on edges, where it bites
hardest — an edge probe \([10,20,30,40,50]\) must come home unchanged,
the crossing edge contributing 30 **once** — and on a proper subset
\(S=\{6,3,4\}\) declared in non-global order.

# Level 6 — the discrete operator

The production Laplacian traverses whatever graph it is handed, so each
part answers over its **own** incidence — and its answers at borrowed
members are wrong on purpose:

| part | global | status | q | L q | A q | authoritative? |
|---|---:|---|---:|---:|---:|---|
| G1 | 1 | owned | 1 | 1 | 1 | **yes** |
| G1 | 2 | owned | 2 | 1 | 3 | **yes** |
| G1 | 3 | owned | 4 | 1 | 7 | **yes** |
| G1 | 4 | *borrowed* | 7 | −3 | **17** | **NO** — global says 13 |
| G2 | 4 | owned | 7 | 1 | 13 | **yes** |
| G2 | 5 | owned | 11 | 1 | 21 | **yes** |
| G2 | 6 | owned | 16 | −5 | 37 | **yes** |
| G2 | 3 | *borrowed* | 4 | 3 | **5** | **NO** — global says 7 |

Keep only owned outputs and the pieces reconstruct the global action:

\[
\boxed{A_G\,q \;=\; \sum_p P_p^{T}\,O_p\,A_{G_p}\,P_p\,q}
\]

\(P_p\), \(O_p\) are *mathematics*, not production types. Level 5 stated
the read/write law structurally; **Level 6 proves it numerically**:

```text
G1: q(borrowed global 4) += 10  →  owned A at global 3:  7 → −3
G2: q(borrowed global 3) += 10  →  owned A at global 4: 13 →  3
restore the halo → the correct answers return
```

A state of the right *size* on a foreign carrier is refused by identity.

# Level 7 — minimization, and what the host does

Production GMRES attached to the Level-6 action on \(G\) solves
\(Aq=b\) to \(q^{*}\), with affine constant zero. Then the level's real
question, answered **behaviourally**:

```text
G      chain:  1→2→3→4→5→6
G_alt  star :  1→2, 1→3, 1→4, 1→5, 1→6      same 6 vertices, same 5 edges

solver_G   % matvec(q*)  =  b            ✓
solver_alt % matvec(q*)  ≠  b            ← nothing changed but the host
```

The same operation type, which stores no graph, gives different
mathematics on a different host.

> The finding is **not** "GMRES traverses the graph" — it does not. It
> is: **GMRES carries the graph to the attached operation, and that
> operation consumes the topology.**
>
> ```text
> minimizer               graph as conduit / context carrier
> differential operation  graph as numerical topology operand
> ```

> **REVIEW GATE B** — after Level 7.

---

# Level 8 — the partitioned constitution

\[
\boxed{\text{SAME MATHEMATICS, DIFFERENT COMPUTATIONAL CONSTITUTION.}}
\]

Level 6 established the local discrete law. Level 8 **composes** it into
a complete realization of the global operation, and the three things
routinely conflated as "partitioning" are kept apart by the type's very
shape:

\[
\boxed{
\begin{array}{lll}
\textbf{STRUCTURAL PARTITION} & G\to\{G_1,G_2\} & \textbf{once} \\
\textbf{NUMERICAL OVERLAP REFRESH} & q\to\{q_1,q_2\} & \textbf{every apply} \\
\textbf{OWNED ASSEMBLY} & \{A_1q_1,A_2q_2\}\to Aq & \textbf{every apply}
\end{array}}
\]

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
```

The composite owns structure and **no mutable numerical state** — no
cached `q`, halo, residual or previous answer. Nothing can go stale
because there is nothing to go stale, and that is *proved*: one instance
applied five times (\(q^{*}\), mixed, \(e_3\), \(e_4\), \(q^{*}\))
returns to its first answer, every intermediate matches the global
action at the time it was asked, and the interleaved probes genuinely
differ.

Extensional equality holds on four probes — \(q^{*}\), a mixed-sign
vector, and both **interface basis vectors** \(e_3,e_4\), which force
information across the cut.

**Two kinds of graph-dependence**, and Level 8 is the second:

| | Level 6 `shifted_laplacian` | Level 8 `partitioned_shifted_laplacian` |
|---|---|---|
| relationship to a graph | **graph-parameterized** | **decomposition-context-bound** |
| acts on | whatever graph it is handed | only the \(G\) its parts were cut from |
| Level 7's star | a legitimate different problem | **refused**, in `domain()` |

# Level 9 — the statement

The user says *solve \(Aq=b\)*; the implementation chooses the
partitioned realization. Equivalence is required first at
`solver % matvec` — the seat GMRES consumes — on all four probes, both
affine constants zero. Then both roads are solved independently:

\[
\boxed{q_{\text{partitioned}} = q_{\text{global}} = q^{*} = [1,2,4,7,11,16]}
\]

```text
PARTITIONED_PDE_RESULT =  1.0000000000000002E+00  2.0000000000000009E+00
                          4.0000000000000009E+00  7.0000000000000027E+00
                          1.1000000000000002E+01  1.6000000000000004E+01
```

One marker, six real tokens in global order, **unrounded** — the honest
image of the arithmetic. `check_marker.sh` validates shape and syntax
only; whether those numbers *are* \(q^{*}\) is the Level-9 test's
business.

> **REVIEW GATE C** — after Level 9. Tower sealed.

## Why this is NOT a distributed solver

```text
one process                    one global trial vector
partition_data reads global state
parts executed sequentially    assembly in-process
global-array inner products    global-array norms
global Krylov state
```

Not claimed: distributed GMRES, MPI solver, parallel halo exchange,
distributed vectors. The correct names are **partitioned matrix-free
solve** and *serial semantic model of a partitioned operator*. What a
genuinely distributed road would need is **derived, not implemented**,
in [`NUCLEUS-OBSERVATIONS.md`](NUCLEUS-OBSERVATIONS.md) under PIP-8.

---

# D. Graph roles at radius 2

This tower exercises the **legacy ordinary-graph / HPC branch** of the
nucleus and introduces no (S, P) view, because its mathematics
does not need one.

| Object | Type | Role |
|---|---|---|
| **global G** | `stored_graph` | topology **yes**; global domain source **yes**; GMRES host **yes**; downstream numerical influence **yes** (Level-7 star control) |
| **parts G1, G2** | `directed_stored_graph` + `partition_relation` | stands in r_p **yes**; local domain source **yes**; local topology operand **yes**; ownership/global maps **yes** — the graph carries r, the four verbs are handed it |
| **borrowed member** | a member of a part carrier | visible to the local operator **yes**; authoritative output contributor **NO** |
| **vertex/edge sets** | `counted_set` / `subset_set` | value domains — what fields live on |

```text
graph as TOPOLOGY         what the Laplacian traverses          (L6)
graph as CONDUIT          what the minimizer carries so its
                          action can traverse something         (L7)
graph as PARTITION FRAME  global↔local maps and ownership       (L4)
member_set as DOMAIN      what a field lives on; never a graph
```

---

# E. What this tower proves / does not prove

## Proven

```text
carriers precede structure, and numerals never establish identity
primitive incidence and intended ownership are relations, not a graph
relation algebra exposes TWO total edge-ownership policies, not one:
    TailOwner = Own^T ∘ Tail   and   HeadOwner = Own^T ∘ Head
unique ownership does NOT imply tail ownership - both policies have it
production SELECTS the tail-based policy, and the crossing edge proves
    which one it selected
the ordinary graph realizes the relation structure extensionally
presence is not ownership; overlap is visibility
a total single-valued ownership policy is what makes unique
    reconstruction possible - not tail-ownership specifically
read = overlap, write-back = owned; one entity, one contribution
local topology actions assemble to the global action, exactly
borrowed input is numerically load-bearing across the cut
the minimizer's graph host is a real conduit to a Class-1 consumer
a partitioned matrix-free action is extensionally the global operator
GMRES solves through it: q_part = q_global = q*
```

## Not proven — anywhere

```text
MPI communication              asynchronous halo exchange
distributed vectors            distributed dot products
distributed Krylov orthogonalization
load balancing                 scalability
communication hiding           parallel performance
multi-rank failure handling
```

---

# F. Code map

```text
test/partitioned-implicit-pde-tower/
├── README.md                     this document — the Rosetta stone
├── NUCLEUS-OBSERVATIONS.md       the evidence ledger (PIP-*), by level
├── run.sh                        level-by-level runner; gates are separators
├── check_imports.sh              fail-closed allowlists, PER LEVEL,
│                                 + its own --selftest
├── check_marker.sh               the result contract + its self-test
├── common/
│   ├── partitioned_pde_assert.f90                 (below everything)
│   ├── chain_carriers_fixture.f90                 earned at Level 0
│   ├── chain_relations_fixture.f90                earned at Level 1
│   ├── chain_algebra_fixture.f90                  earned at Level 2
│   ├── shifted_laplacian_fixture.f90              earned at Level 6
│   └── partitioned_shifted_laplacian_fixture.f90  earned at Level 8
├── level-0-carrier/              test.f90
├── level-1-relation/             test.f90
├── level-2-relation-algebra/     test.f90
├── level-3-graph/                test.f90
├── level-4-graph-calculus/       test.f90
├── level-5-field-calculus/       test.f90
├── level-6-discretization/       test.f90 · refusal.f90 · check_refusals.sh
├── level-7-minimization/         test.f90
├── level-8-constitution/         test.f90 · refusal.f90 · check_refusals.sh
└── level-9-statement/            test.f90
```

`common/` is **not** a hole in the stratification: the import gate keys
its allowlists **per file**, so each shared fixture is bound to the level
that earns it — Level 5 cannot reach the Level-6 operator, and Level 7
cannot reach the Level-8 composite.

The fixture ladder is the tower's own stratification applied to itself:

```text
Level 0    chain_carriers_fixture      declares V, E, K
Level 1    chain_relations_fixture     states Tail, Head, Own over them
Level 2    chain_algebra_fixture       composes what follows
```

The relation fixture does not *import* the carrier fixture — its
constructors receive \(V,E,K\) as arguments, because a Level-1 file may
state facts over sets but may not name a set into existence. The ladder
is realized by the tests' imports and by the per-file allowlist keying,
and it is enforced in **three** places at once: the import gate's
allowlists, each level's Makefile (Level 0 compiles and links the
carrier fixture and *never* the relation fixture), and
`check_imports.sh --selftest`, which asserts directly that a Level-0
source saying `use chain_relations_fixture` is refused.

---

# G. Status

```text
PARTITIONED IMPLICIT PDE TOWER

    L0 carrier ........................ PASS
    L1 relation ....................... PASS
    L2 relation algebra ............... PASS
    L3 graph .......................... PASS
    L4 graph calculus ................. PASS

    ===== REVIEW GATE A =====

    L5 field calculus ................. PASS
    L6 discretization ................. PASS
    L7 minimization ................... PASS

    ===== REVIEW GATE B =====

    L8 constitution ................... PASS
    L9 statement ...................... PASS

    ===== REVIEW GATE C =====

    PARTITIONED_PDE_RESULT = 1.0000000000000002E+00 ... 1.6000000000000004E+01

    TOWER SEALED.
```

**Zero production changes** beyond a single corrected comment earned at
what is now Level 4.
