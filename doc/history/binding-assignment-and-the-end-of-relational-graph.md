# Binding assignment, then the end of relational_graph
From `50098a9`. `src/fractal_graph.f90` unchanged throughout, verified at each commit.
```
621fcb7  GATE       assignment is refused; the deep copy goes
997d8c8  A          construction + validity
870920d  B          access + identity: six towers at level 3
42b99ac  C+D+FINAL  the remaining consumers, the allowlists, the deletion
```

## 1. The lender-overwrite result
**The published law was stronger than the implementation.** `b = d` invoked a deep copy whose
left-hand side is `INTENT(OUT)`, so `b` was finalized on entry and its objects freed before the copy
began — `copy_binding` → finalization → `release_binding` → `free`. A pointer borrowed from `b`:
```
native   : associated, and answers R2      <- silently wrong, no diagnostic
valgrind : invalid read of a freed block, and answers R1
```
Two allocators, two wrong answers, and the native one is invisible to every tool: the freed block is
reused by the very next allocation inside the same copy, so the pointer is live and denotes a
different relation.

## 2. Extension, replacement, and whether Fortran can prohibit one
`bind_*` extends a binding and preserves what it lent. **No replacement can** — the objects it lent
are the ones it must free. So replacement is not an operation, and is refused. Five prohibition
mechanisms, compiled and rerun on every pass in `test/graph-relational/fortran-assignment`:

| mechanism | measured |
|---|---|
| private generic `assignment(=)` | **admits** — intrinsic assignment substitutes, silently |
| generic with no matching specific | **admits** — same |
| coarray ultimate component | **admits**, undiagnosed at `-Wall -pedantic`; and its pointer form rejects the *declaration* the callers need |
| `LOCK_TYPE` component | rejects the declaration, not the assignment |

The two silent fallbacks are worse than nothing: for a type whose rows hold pointers, intrinsic
assignment is a shallow copy and a double free. So the prohibition is at run time, and the law says
so. `intent_out_finalizes` holds the rule the refusal itself rests on — an `INTENT(OUT)` dummy of a
finalizable type is finalized on entry, so the refusal takes **`INTENT(INOUT)`**, or it would destroy
the lender before refusing to.

## 3. The storage law, and the copy that went with it
> `relational_binding` is not assignable — refused at run time, since no Fortran mechanism prohibits
> it at compile time; `bind_*` preserves every outstanding object pointer until it is destroyed.

No caller ever assigned a binding — `graph_profile` takes it `intent(in)`, the towers declare locals.
`copy_binding` existed only because individually allocated pointer objects would otherwise be
shallow-copied and freed twice; refusing removes the premise. `lifetime.f90`'s cases D and E existed
only to exercise it: D is now two coexisting bindings, and the refusal moved to `refusal.f90`.

## 4. Where create_graph's six laws went
The container refused at construction because it was a constructor. A view is not one, so the laws
divide by what they are about.

| old refusal | now |
|---|---|
| a foreign domain | `relational_valid` → `.false.`, already |
| the same domain twice, the same relation twice | `relational_valid` → `.false.` — **new** |
| a graph never signs twice | the kernel's `declare`, already |
| a domain that never signed, a borrowing view offered | `bind_*` refuses — **new** |

S and P are **sets**; the branches represent them as sequences, and a sequence may repeat what a set
cannot. That gap is a property of a graph already built, so it is answered. The last row is storage
law: an object with no identity cannot be compared, and a borrowing view copied into owned storage
carries a reference to a base the binding does not keep alive.

## 5. No shared builder, and one suite deleted
§12 asked that the fixture die without becoming one, so **each suite constructs the representation it
means, inline** — twelve lines of declare, bind and wire. The level-3 suites are *about* the
representation; hiding it behind a builder would make the recursion invisible in exactly the places
that test it. Zero helper LOC.

`test/graph-structure` is gone: its subject was the container's implementation. Every law it held has
a home — counts, seats by identity and `holds_set` were in `graph-relational` already;
ternary-beside-binary and same-signature coexistence moved there; stable borrows are `lifetime.f90`
case E; identity is `fractal-graph`. Retargeted: `graph-ordinary`, `graph-algorithms`, and every
level-3, -4, -6, -8 and -9 suite across six towers. Where `graph_grammar`'s own type named `graph` met
the kernel's, the kernel keeps the name.

## 6. The measurement
| | lines |
|---:|---:|
| `src/graph_structure.f90` removed | **344** (164 code) |
| `test/support/relational_fixture.f90` removed | **196** |
| shared helper added | **0** |
| `src/graph_relational_view.f90` | 371 → 408 (187 → 195 code) |
| **net src** / net test | **−308** / +1284, all rewrites |
| **types eliminated** | **3** — `relational_graph`, `held_set`, `held_relation`, plus `create_graph` |

## 7. Audit
`grep -R` over `src test` is empty for `use graph_structure`, `relational_graph`, `held_set`,
`held_relation`, `create_graph` and `relational_fixture`. The tower READMEs and AGENTS.md §0 and §58
describe what is there; the NUCLEUS-OBSERVATIONS ledgers keep the old names as the dated records they
are. `git diff 50098a9 -- src/fractal_graph.f90` is empty. Every suite matches; `graph-benchmark`
stays red on the pre-existing missing `class_graph_support.mod`. valgrind on `lifetime`: 279 allocs,
279 frees, 0 errors. `graph-relational` 30 → 49 assertions.

## 8. Next, and not yet
`member_set` and `relation` are the two remaining roots, both more substantial than the container
was. Each needs its own map before any code moves — is `member_set` ontology, a view of graph, an
attribute domain, or a compiled representation? is `relation` ontology, a view, or an operator? Then
`interface_graph.f90`'s abstract `graph`, an interface rather than a container.
