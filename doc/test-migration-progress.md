# Test migration: 18 of 33 suites, incremental

Branch `tower-domain-cutover`. **Not finished** — 15 suites remain, so the §8
deletion gate is not open and `graph_carrier` still stands. Every commit below
is one suite, built and run green before committing; nothing was swept.

## Migration table

| suite | files | obsolete fixtures | notes | result |
|---|---|---|---|---|
| graph-multigrid | 2 | — | production fix: count components lacked defaults | PASS |
| graph-robustness | 1 | — | needs **no map at all** | PASS |
| graph-constitution | 1 | — | local maps for a local carve | PASS |
| graph-mesh | 1 | — | local maps | PASS |
| graph-partition | 1 | — | binds carriers before transport | PASS |
| graph-field-transport | 1 | — | 2 `select type` collapsed | PASS |
| graph-field | 2 | — | transitive `declared_subobject` | PASS |
| graph-marching | 4 | — | production fix: `attach` re-enterable | PASS |
| graph-ordinary | 3 | **listed_set** dropped | fixture had nothing left to prove | PASS |
| fractal-map | 1 | — | scale argument restated | PASS |
| graph-characterization | 1 | — | carved edge sets | PASS |
| graph-algorithms | 2 | — | `sources`/`sinks` carve; shadowing bug | PASS |

Already green before this pass, untouched: `fractal-graph`, `graph-sequence`,
`graph-set-view`, `graph-state`, `graph-inclusion`, `graph-identity-map`.

**Remaining (15):** graph-minimization, graph-carrier, graph-binary,
graph-algebra, graph-relation, graph-contract, graph-relational, the seven
towers, and graph-benchmark.

## What the suites taught, that a sweep would have hidden

**Two production defects, both found by *running*, not compiling.**

- `n_unknown_domain` / `n_residual_domain` were added to `minimizer` without
  default initializers, so `jacobi()` and `gmres()` — structure constructors
  naming no component — stopped compiling. Fixed at the declaration.
- **`attach` is re-enterable.** Newton calls it once per iteration, and it
  declared the minimizer's number domain each time; a graph signs once, so the
  second iteration died. The *old* code minted a fresh domain per attach
  (`counted_set`'s constructor declared one), so that semantics is preserved by
  resetting the component to an unsigned graph before declaring. Which of the
  two readings was meant is now stated rather than implied.

**A map that describes nothing compiles perfectly well.** `graph-algorithms`
had a check given its own local `sets`, shadowing the program's — the map that
actually described the domain. It built cleanly and died at run time inside the
CSR build. This is the §5 scope law, and it is not enforceable by the compiler;
it needs the run.

**`graph-partition` failed first time with `no representation describes that
set`, correctly.** The test declared its domains without ever saying what they
contain, and `partition_data` asks where a global member stands. That is the
explicit-dependency law working, not a defect. All measured laws survive
unchanged: `sum_k A_k(P_k(D)) = D` at nparts 2 and 3 on vertex *and* edge data,
every cell owned exactly once, identity at nparts=1 and not at nparts=2.

## §2 — dead abstractions removed, not ported

- **`listed_set_fixture` dropped from graph-ordinary's build.** It existed to
  prove the relation generic over `member_set` concretions; that generality now
  lives in `set_representation`, of which the library ships two. Its sparse
  unordered-world domain is one `listed_set_representation` bind. (The file
  still exists because `graph-relation` also builds it — it dies with that
  suite.)
- **Two `select type` blocks collapsed** in graph-field-transport. They asked
  whether a transported domain was a `subset_set` *type*, then whether it
  embedded where it should. With one domain type the type arm is meaningless
  and the `class default` arm unreachable; what survives is
  `declared_subobject(dp_, part % vertex_set(), inclusions)` — the assertion the
  old code made *after* the type test, now carrying the whole meaning.

## §4 — the three distinctions, kept apart

No assertion was rewritten into a different question. `same_as` still answers
identity (graph-field's ambient domain, and its copy); `sets % index_in` answers
storage position; `declared_subobject` answers provenance, including
transitively — graph-field's `hot` descends through `walls` to `cells`, and the
inclusion map walks it.

Suites that ask only *which* and *how many* carry **no map at all**:
`graph-robustness` imports `set_graph` and nothing else, because the graph
already answers both.

## Gate status (§8) — not open

    member_set   src 30   test 372
    counted_set  src 12   test 531
    subset_set   src 17   test 262
    use graph_carrier   src 0   test 117

Every `src/` hit outside `src/graph_carrier.f90` is a comment, and **no src
module imports `graph_carrier`** — it is dead production code held alive
entirely by the 117 test files that still `use` it.

## §9 / §10

`graph-benchmark` fails identically to its pre-existing baseline —
`class_graph_support.mod` missing, from `a3817e3`, long before this branch. Not
a new failure hidden under the exemption; it is the same failure, and the suite
was not otherwise touched.

`src/fractal_graph.f90` is byte-identical to `57d8c51`.

## Stop reason

Neither §11 A nor B: no suite exposed a new production semantic dependency, and
no adapter was invented. I ran out of room with 15 suites left. The next pass
should continue in the same mode — the towers are the bulk, and their fixtures
implement `graph_operation`, whose shape `graph-marching` now demonstrates.
