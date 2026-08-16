# Invariant closure
Spike only; `src/` is unchanged. Run by `test/fractal-graph/run.sh` (gfortran-15, `-std=f2023`).
Every claim is a fixture in `fortran-recursion/`, compiled against the shipped kernel, not a copy.

## 1. Final declarations
```fortran
integer, parameter :: GRAPH_NULL = 0, GRAPH_UNKNOWN = 1, GRAPH_KNOWN = 2

type :: graph_branch
   private
   integer              :: status_ = GRAPH_NULL
   type(graph), pointer :: known_  => null()
 contains
   procedure :: status, known
end type graph_branch

type :: graph
   type(graph_branch)   :: branch(2)
   type(token), private :: identity
 contains
   procedure :: declare, id, same_as
end type graph
```

## 2. Which components are private
`graph_branch % status_`, `graph_branch % known_`, `graph % identity`. Both branch components
must be private: public `status_` admits `status = KNOWN` with no reference, public `known_`
admits `nullify` under a KNOWN status.

`graph % branch` stays **public**: every value a caller can place there is lawful, since the
structure constructor is unavailable outside the module (`private_structure_constructor`) and the
three constructors are total on the invariant. Private `branch` would cost an accessor and, since
a function result cannot be a part-ref base (Q4), worsen navigation.

## 3. Public navigation syntax
```fortran
g % branch(i)                              ! public component
g % branch(i) % status()                   ! branch state
g % branch(i) % known()                    ! reference; disassociated unless KNOWN
p => g % branch(i) % known()                ! pointer binding, for depth
associate (x => g % branch(i) % known())    ! ASSOCIATE name, for depth
associated(g % branch(1) % known(), g % branch(2) % known())
```
`encapsulated_navigation` compiles all of these against the kernel.

## 4. Can chained recursive navigation compile?
**No**, and the obstruction is Fortran's. A function reference is not a part-ref, so no data-ref
may be built on one:

| form | gfortran-15 `-std=f2023` |
|---|---|
| `null_branch() % status()` | *The leftmost part-ref in a data-ref cannot be a function reference* |
| `g % branch(1) % known() % branch(1) % status()` | *Unclassifiable statement* — same rule, from the parser |

Depth is reached by a pointer or an ASSOCIATE name, one level at a time. Not the flattening §16
rejects: `branch` stays in the public structure, `type(graph_branch) :: branch(2)` stands in the
declaration, `g % branch(i) % status()` / `% known()` remain the navigation; only a multi-level
chain *within one expression* is unavailable. `graph_views` is unaffected — it passes
`g % branch(s) % known()` as an argument, one function reference, not a chain.

## 5. How legal branches are constructed
```fortran
null_branch()       status = NULL      known disassociated
unknown_branch()    status = UNKNOWN   known disassociated
known_branch(g)     status = KNOWN     known associated with g; refuses undeclared g
```

The only introductions of a branch value; `graph_branch(...)` is rejected outside the module, so
status and reference cannot be set independently. A branch is replaced whole,
`g % branch(1) = known_branch(h)`; `known_branch` cannot be `pure` (F2018 C1594).

## 6. How cycles are closed
By whole-branch assignment once both graphs exist — the construction mutation the previous report
established as necessary:

```fortran
call a % declare(); call b % declare()
a % branch(1) = known_branch(b)          ! a -> b
b % branch(1) = known_branch(a)          ! b -> a
```

No cycle-closing operation is introduced; sharing is the same act applied twice to one target.
Both survive encapsulation (tests B and C, `associated` included, so no copy).

## 7. Can published graphs be corrupted externally?
**The invariant, no. The shape, yes — deliberately.**

- `status_` and `known_` cannot be assigned, rebound, or constructed from outside; three fixtures
  assert the rejection and its reason. Every reachable branch satisfies
  `status() == GRAPH_KNOWN .eqv. associated(known())`, asserted over all nine states.
- A whole branch **can** still be replaced after publication (`g % branch(1) = null_branch()`) —
  semantic mutation, not invariant violation, since the replacement is lawful. §7 declines to
  police it here; the recommendation stays construction-mutable, publication-immutable.
- Lifetime is unchanged: branches do not own targets, `TARGET` belongs on the actual argument,
  reclamation is by region.

## 8. Kernel code lines
54 → **67** code (comment 35 → 49, blank 43 → 54), by
`grep -c -v -E '^[[:space:]]*(!|$)' fractal_graph.f90`. +13 against a limit of 80, for `private`
and two accessors. No type was added.

## 9. Counts
| | before | after |
|---|---:|---:|
| compile-time fixtures | 3 | **9** — 3 admitted, 6 rejected with their reason |
| runtime assertions | 39 | **45** — all 39 retained, 6 added for the invariant |
| runtime refusals | 8 | 8 |
| total PASS | 50 | **62**, 0 FAIL |

The three invariant protections are compile-time refusals — they never link, so they cannot join
`./refusal`; they are asserted as `graph_before_branch` is, by required failure with a diagnostic.

## 10. Recommendation
Adopt. Two types, no third; `branch(2)` public and visible; `status_`/`known_` private — the
smallest change that closes `KNOWN iff associated(known)`, at compile time rather than by
convention. Its one cost, a multi-level traversal needing a pointer or an ASSOCIATE name, is a
Fortran rule about function results, held by two fixtures that fail if it ever changes.

Rejected here against §14 and §17: private `branch`, a realize operation, lifecycle machinery.
The next question is `seal`, closing semantic mutation and making caching on identity sound.
