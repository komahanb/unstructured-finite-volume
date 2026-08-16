# Test migration: 18 of 33 suites green — stop condition C

Branch `tower-domain-cutover`, head `279d869`. Twelve suites migrated across two
passes, each built, run and committed on its own. **No tower migrated this pass**
— the reason is a finding, below. `graph_carrier` still stands; the §9 gate is
shut.

## Migration table

| suite | files | obsolete fixtures | notes | result |
|---|---|---|---|---|
| graph-multigrid | 2 | — | production fix + both repairs pinned | PASS |
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

Green before these passes, untouched: `fractal-graph`, `graph-sequence`,
`graph-set-view`, `graph-state`, `graph-inclusion`, `graph-identity-map`.

## §5 and §6 — the two production repairs, now pinned

Both were found by running, and neither had a test naming it. Both are pinned in
`graph-multigrid`, the smallest suite that already exercised them:

- **`attach` is re-enterable.** Newton calls it once per iteration and it
  declared the solver's number domain each time; a graph signs once, so the
  second iteration died. The old `counted_set` constructor minted a *fresh*
  domain per attach, so that reading is kept — reset to an unsigned graph, then
  declare. The regression attaches a second time and asserts the same solution:
  *"attach is re-enterable: the second attach solves the same"*.
- **A count component needs a default.** `jacobi()` and `gmres()` name no
  component, so `n_unknown_domain` / `n_residual_domain` must have values without
  being named. Asserted directly: *"an unattached solver counts no unknowns, and
  says so"*.

No unrelated production cleanup was mixed in.

## Why no tower moved — a finding, not a shortfall

Each tower carries **`check_imports.sh`, an import gate** that names, per level
and per fixture, exactly which modules that level may import — `graph_carrier`
among them. It is the tower's layering law expressed as a script, and it is not
incidental machinery: migrating a tower means restating that law over
`fractal_graph` / `graph_set_representation` / `graph_set_map`, deciding for
each level whether it may now see a *map* at all.

That is a design question per tower, not a mechanical edit, and it is the reason
a tower is not simply "the biggest of the remaining suites". The fixture work
itself is straightforward and `graph-marching` remains the pattern: an operation
fixture stores `type(set_graph) :: state` plus its count, and `domain` answers
both. I began `time-integration-tower`'s common fixtures on that pattern, then
reverted them rather than leave a tower half-migrated and uncommitted.

**Recommendation for the next pass:** rule on the import gate first — whether a
level that declares carriers may import `graph_set_map` — then the towers become
mechanical. `time-integration-tower` is the smallest (16 files) and asks no
domain labels, so it needs `sets` and `inclusions` only.

## §7 — deleted, not ported

`listed_set_fixture` left `graph-ordinary`'s build: it existed to prove the
relation generic over `member_set` concretions, and that generality now lives in
`set_representation`. Two `select type` blocks collapsed in
`graph-field-transport` to `declared_subobject(...)` — the assertion the old code
made *after* the type test, now carrying the whole meaning.

## §8 — partition laws, unchanged

`sum_k A_k(P_k(D)) = D` at nparts 2 and 3, on vertex **and** edge data; every
cell owned exactly once; identity at one part and **not** at two. Nothing is
called an inverse.

## §9 — gate shut

    member_set   src 30   test 372
    counted_set  src 12   test 532
    subset_set   src 17   test 262
    use graph_carrier   src 0 files   test 117 files

Every `src/` hit outside `src/graph_carrier.f90` is a comment, and no src module
imports it: dead production code held alive entirely by 117 test files.

## §10, §11

`graph-benchmark` fails **exactly** the baseline — `class_graph_support.mod`
missing, from `a3817e3`. Same failure, not a new one.

`src/fractal_graph.f90` is byte-identical to `57d8c51`.

**Stop condition C.** Migrated 18, remaining 15: `adjoint-tower`,
`calculator-tower`, `derivative-action-tower`, `graph-algebra`,
`graph-benchmark`, `graph-binary`, `graph-carrier`, `graph-contract`,
`graph-minimization`, `graph-relational`, `graph-relation`, `learning-tower`,
`partitioned-implicit-pde-tower`, `time-integration-tower`,
`visualization-tower`.
