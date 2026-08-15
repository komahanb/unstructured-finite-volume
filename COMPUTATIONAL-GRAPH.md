# COMPUTATIONAL-GRAPH.md — The Epistemic Pair

## Mission

Name the computational graph, and keep it distinct from the structural one.

This document reserves terminology only. It decides no levels, no storage,
no algorithms, no inheritance. Those decisions come later, and separately
(AGENTS.md, releveling report).

---

# 1. Two graphs, one word

The repository now holds two different mathematical objects that both
answer to the word "graph", and a third meaning that interprets the first:

```text
structural graph        Gamma = (S, P)     relational_graph
                                           (graph_structure.f90)

ordinary interpretations                   ordinary_graph_view,
of structural schemas                      directed_adjacency_view
                                           (graph_profile.f90), and the
                                           legacy abstract graph of
                                           graph_grammar.f90

computational graph     G = (Q, R)         computational_graph
                                           (graph_state.f90)
```

Do not write the bare word `graph` in new code or documentation when one
of these three meanings can be stated precisely.

The structural graph:

\[
\boxed{\Gamma=(\mathcal S,\mathcal P)}
\]

is a structured collection of member sets \(\mathcal S\) and typed
structural relations \(\mathcal P\). It says which members stand
together. It is the container of AGENTS.md 14 and nothing here changes
it.

The computational graph:

\[
\boxed{G=(Q,R)}
\]

is an epistemic pair over such a structure:

- \(Q\) — the **data**: all computational data participating in the
  problem;
- \(R\) — the **residual/operator**: the governing
  constraint/residual/operator over that data.

The letter \(R\) is reserved for this residual. Structural relations are
written \(P\), \(T\), \(H\), \(A\), \(\ldots\) — never \(R\) — and the
structural graph is never written \(G=(S,R)\).

The old incidence-first proposal \(G=(A,B,I)\) is superseded and stays
superseded: neither \(\Gamma\) nor \(G=(Q,R)\) revives a privileged
incidence representation. Incidence and adjacency remain ordinary
structural relations inside \(\mathcal P\).

---

# 2. Bottom is not empty

\(\bot\) means **unknown/unrealized** — a seat with no occupant. It never
means empty.

\[
\boxed{\bot\neq\varnothing}
\]

An allocated data object with zero entries is **realized data**: the
knowledge is present and its content happens to have zero length (the
empty subset has always been a valid domain — AGENTS.md 63). An
unallocated, undeclared seat is \(\bot\): nothing is known. The same
distinction holds for the residual seat.

Likewise `void` refers to computational knowledge, not topology. A void
computational graph may ride on a perfectly well-defined, fully
materialized \(\Gamma\). It is not an empty relational graph; structural
emptiness and epistemic unknownness are different absences.

---

# 3. The four states

Each seat is realized or it is \(\bot\). Two seats, four states, each
with exactly one canonical name:

| State | Pair | Canonical name |
|---|---|---|
| neither realized | \((\bot,\bot)\) | **void graph** |
| data only | \((Q,\bot)\) | **data graph** |
| operator only | \((\bot,R)\) | **operator graph** |
| both | \((Q,R)\) | **realized graph** |

Exactly one holds, always:

\[
\begin{aligned}
\text{void}      &\iff \neg Q\land\neg R,\\
\text{data}      &\iff Q\land\neg R,\\
\text{operator}  &\iff \neg Q\land R,\\
\text{realized}  &\iff Q\land R.
\end{aligned}
\]

Use these four names in comments, tests, documentation, constructors and
diagnostics. Do not coin synonyms. In particular, none of:

```text
empty graph      partial graph     physics graph     solution graph
complete graph   residual graph    inverse graph     forward graph
```

names any of the four states. Those words mean other things, or smuggle
in an interpretation the state does not carry.

## 3.1 Data graph

\((Q,\bot)\): data is realized; the governing residual is not.
Experimental observations, trajectories, measured fields, state
histories, parameter-response data. \(Q\) is not restricted to a PDE
solution. The data graph is the natural starting point for operator
discovery.

## 3.2 Operator graph

\((\bot,R)\): the residual is realized; the data it governs is not.
\(R(Q,\theta)=0\) with \(R\) known and \(Q\) sought. Not restricted to
differential equations: \(R\) may be any residual/constraint/operator.
The operator graph is the natural starting point for forward solution.

## 3.3 Realized graph

\((Q,R)\): both seats occupied. **Realized does not mean solved.** Do
not assume \(R(Q)=0\). A realized graph may hold a converged solution, an
unconverged iterate, observed data under a candidate operator, a
candidate discovered model, or a deliberately inconsistent pair. When
asserting conditions on \(R(Q)\), use separate vocabulary — `satisfied`,
`consistent`, `converged` — and assert them separately. Never rename
`realized_graph` to `solution_graph`.

---

# 4. One type, four states

The four states are one finite \(2\times2\) classification:

```text
has Q?   yes / no
has R?   yes / no
```

Under the admission laws this earns state constants, not a type
hierarchy:

```fortran
type :: computational_graph
contains
   procedure :: has_data
   procedure :: has_operator
   procedure :: state
end type
```

with

```fortran
GRAPH_STATE_VOID
GRAPH_STATE_DATA
GRAPH_STATE_OPERATOR
GRAPH_STATE_REALIZED
```

The public vocabulary provides constructors named for the states —
`void_graph`, `data_graph`, `operator_graph`, `realized_graph` — and all
four construct the same underlying `computational_graph`. Distinct
derived types for the four states are admitted only if a later
architectural analysis proves genuinely different roles or contracts.
Inheritance is never used merely to encode known/unknown status.

---

# 5. Public terminology

Mathematics writes \(Q\), \(R\), \(\Gamma\). Public Fortran interfaces
speak descriptive names:

```fortran
g % data()
g % residual()
g % structure()
```

never `g % q()` or `g % r()`. Single-letter notation belongs in
mathematics and local algebra, not the public semantic contract. For
\(R\), prefer `residual` or `residual_operator`; avoid bare `operator`,
which collides with the legacy generic `graph_operation`.

## 5.1 Q is not a graph_field

\(Q\) is conceptually all computational data participating in the
problem: state fields, parameters, observations, coordinates, history,
boundary data, experimental data. A `graph_field` may be one constituent
of \(Q\); the ontology does not define \(Q\) as one `graph_field`, and
the naming pass deliberately leaves the storage of \(Q\) unfixed.

## 5.2 R is not graph_operation

The existing `graph_operation` is the broader legacy role: a verb mapping
graph data to graph data. \(R\) is semantically the governing
residual/operator. Whether `residual_operator is-a graph_operation`, or
some other relationship holds, is an inheritance question deliberately
left open. This document reserves the words; it does not decide the
hierarchy. Do not globally rename `graph_operation` to `operator`.

---

# 6. Canonical computational directions

Processes move between states; states do not encode processes. `forward`
and `reverse` never enter the graph type names.

## Solution

\[
(\bot,R)\rightarrow(Q,R)
\]

Canonical wording: `solve for data`, `infer data`, `realize data`.
`solve` is appropriate when \(R(Q)=0\) — or another explicit governing
condition — is being satisfied.

## Operator discovery

\[
(Q,\bot)\rightarrow(Q,R)
\]

Canonical wording: `discover operator`, `infer operator`,
`realize operator`.

## Joint inference

\[
(\bot,\bot)\rightarrow(Q,R)
\]

when external statements, measurements, priors or constraints provide
sufficient information.

---

# 7. Composition, reserved

\[
(Q,\bot)\oplus(\bot,R)\rightarrow(Q,R)
\]

is the composition of compatible data and operator knowledge. The word
`compose` is reserved for it. No overloaded composition operator is
implemented until the compatibility laws are established:

```text
same structural host, or a defined structural mapping
compatible domains
compatible value spaces
compatible residual signature
```

Vocabulary now, semantics when the laws exist.

---

# 8. What this document does not decide

- **Levels.** The tower placement of the computational graph and any
  renumbering of existing levels belongs to the releveling report
  (doc/computational-graph-releveling.md) and a separate commit.
- **Storage of \(Q\).** See 5.1.
- **The operator hierarchy.** See 5.2.
- **Solving, fitting, minimization, adjoints, discovery.** None of it
  enters graph_state.f90; the module holds the vocabulary — the type,
  the state constants, the classification queries — and nothing else.
