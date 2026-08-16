# Identity and mutation
Spike only; `src/` is unchanged, and `fractal_graph.f90` is unchanged by this pass — still 67
code lines. The previous report, on the branch invariant, is at commit `073451d`.
Run by `run.sh`: 86 PASS, 0 FAIL. Fixtures compile against the shipped kernel, not a copy.

## 1. Is graph identity object identity or value identity?
**Object identity.** `declare` mints a token per object; `same_as` compares tokens and nothing
else. This is not a preference. Value identity would be structural equality, and structural
equality on a graph with cycles is not a component comparison — it is a bisimulation over the
reachable set, an algorithm the kernel does not contain and §13 would not admit. `same_as`
answers in constant time; value equality cannot. Model C therefore presupposes an equality this
ontology cannot compute.

Identity answers *which graph*. `status()` answers *what is this graph now*. The kernel keeps
them apart, and nothing below forces them together.

## 2. Does lawful branch mutation preserve identity?
**Yes.** `id(g)` is unchanged across every one of the nine branch transitions (`mutation.f90`, T)
and across mutation of both branches (I). `g % same_as(g)` holds throughout.

## 3. Is publication immutability mathematically necessary?
**No.** The previous pass argued it from caching. Two things dissolve that argument:

- A compiled representation is a **snapshot** (Q9): `C0` stays valid when `G` mutates, and `C1`
  is a second, independent value. Nothing about CSR requires `G` to stop changing.
- A cache that needs invalidation can carry whatever key it needs. That belongs to the cache,
  which knows its own validity conditions, not to a kernel that does not.

Q5 then shows the cost is not local: immutability is a property of an entire reachable set, not
of one graph. Nothing in the ontology requires paying that.

## 4. Can publication immutability be enforced while keeping public `branch(2)`?
**No.** Four candidate mechanisms, each built or rejected on every run (`fortran-mutation/`):

| candidate | result |
|---|---|
| graph-level defined assignment | **not invoked** by component assignment — a sealed graph had `branch(1)` assigned anyway |
| branch-level defined assignment | invoked, but its dummy arguments are two branches; no argument, host association or inquiry yields the owning graph, so the seal cannot be consulted |
| private branch, value-returning accessor | assignment **rejected** — *The function result on the lhs of the assignment must have the pointer attribute*. This does enforce it |
| private branch, pointer-returning accessor | assignment **accepted** — a pointer function result is a definable variable. Enforces nothing |

Only the value-returning accessor enforces anything, and the same language rule that enforces it
removes the navigation: with `branch` private, `g % branch(i) % status()` chains past a function
reference and is rejected (`accessor_navigation`; gfortran words the diagnostic by form, but the
rule is one — a function result is not a part-ref). Constness and the fractal API are not
separable trade-offs here; they are one rule read twice. §9 applies without being invoked.

Fortran's lack of a const reference is the conclusion, not a gap to be papered over.

## 5. Can a graph reached through `known()` be mutated transitively?
**Yes**, on the shipped kernel (`transitive_mutation`): `p => a % branch(1) % known()` then
`p % branch(1) = known_branch(c)` mutates `b`, whose identity is intact. Blocking this would
require every graph in the reachable set to refuse independently — the full cost of Q4's
value-accessor design paid on every graph, not on the one that was sealed.

## 6. Is `declare` identity assignment, publication, or both?
**Identity assignment only.** Cycle construction wires branches *after* both graphs declare:

```fortran
call a % declare(); call b % declare()
a % branch(1) = known_branch(b)
b % branch(1) = known_branch(a)
```

If `declare` froze the graph, that sequence could not exist. *Publication* does not name anything
the current architecture does, and this report does not use the word for `declare`.

## 7. Should `known_branch` require an already-declared graph?
**Yes, keep the requirement.** An undeclared token matches nothing, including itself, so an
undeclared graph is **not `same_as` itself** (`mutation.f90`, R). A KNOWN branch to such a graph
would put a non-reflexive element into the reachable set: equality would fail on it, and a map
keyed on identity would never find it while silently accepting new rows for it. Wiring before
declaring buys a construction-order convenience and costs reflexivity.

## 8. Are branch-state transitions unrestricted, monotone, or interpretation-dependent?
**Unrestricted in the kernel; monotone only where an interpretation says so.** All nine
transitions apply, each leaving the invariant and the identity intact (T). `KNOWN -> UNKNOWN` and
`KNOWN -> NULL` both execute — the adaptive example (X) attaches a cell so a boundary face
becomes interior (`NULL -> KNOWN`), then detaches it (`KNOWN -> NULL`), one face identity
throughout. A monotone-knowledge law would forbid the second transition, and the second
transition runs. The ontology states three states and no transition law; there is no evidence for
inventing one.

Calling that face *the same face in a later state* is coherent: its identity is what makes the
two states comparable at all.

## 9. Does CSR need graph immutability, or is CSR simply a snapshot?
**A snapshot** (`mutation.f90`, S). `G` compiles to `C0 = [1,2,4,5,5]`; `G` gains a connection;
`G` compiles to `C1 = [1,2,4,6,7]`. `C0` is unchanged, internally consistent, and does not follow
`G`. Neither carries a version, because neither needs one: each is the value compilation returned
when it ran. Versioning, if a client ever needs it, belongs to that client.

## 10. Does the kernel need any change?
**No.** `fractal_graph.f90` is byte-identical to `073451d` — 67 code lines. No `epoch`,
`generation`, `version`, `dirty`, `seal` or `frozen` field was added, and no counterexample
appeared that could not be handled outside it.

## Comparison
`t` tested by a fixture or assertion here; `a` argued.

| | A mutable | B sealed | C persistent |
|---|---|---|---|
| sharing | holds `t` | holds while constructing, frozen after `a` | degrades: every path to a change is rebuilt `t` |
| cycles | holds `t` | holds while constructing `a` | no fixpoint `t` |
| identity | stable object identity `t` | stable, plus a phase `a` | changes per mutation; needs an equality the kernel cannot compute `a` |
| direct branch API | `g % branch(i) % status()` `t` | requires private branch; navigation rejected `t` | as B `t` |
| caching | cache carries its own validity `a` | identity alone suffices `a` | each version is a new key `a` |
| CSR | snapshot `t` | snapshot `t` | snapshot per version `a` |
| adaptive graph | same graph, later state `t` | forbidden after the freeze `a` | rebuild per change `a` |
| implementation cost | **zero kernel change** `t` | private branch + seal field + checked setter + value accessor, on **every** graph in the reachable set, and navigation lost `t` | B's cost plus a computable value equality that does not exist `a` |

## Conclusion
Adopt **A**. A graph is a mutable object with stable identity; compiled representations are
snapshots; immutability and versioning belong to the clients that need them; `fractal_graph.f90`
requires no change.

B is not refuted as a *semantics* — it is refuted as something this language can enforce without
destroying the public fractal, and it was never shown necessary. C is refuted at the root: it
needs a value equality the ontology cannot compute.
