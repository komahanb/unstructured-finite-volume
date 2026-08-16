# Test migration: attempted, not completed

**I did not finish this.** `tower-domain-cutover` is unchanged at **`b2fdf03`**,
identical to what was pushed. The half-migration was reverted rather than
committed, so the branch is not worse than it was. Requirement 5's empty grep
is **not** met, and requirements 7 and 10 are therefore not met either.

What follows is what the attempt measured, because the work-list is now exact.

## Where it got to

Two mechanical passes (domain construction → identity + representation + label;
map declaration per scope) plus a third (map questions → the map, `graph`
ambiguity → the established `set_graph` rename, escaping domains → maps as
`intent(inout)` dummies) carried **7 of 33 suites**. The other 26 were left with
~1000 residual diagnostics in classes that need per-file judgement, and at that
point I was debugging regex collateral rather than migrating. That is the wrong
mode, so I stopped and reverted.

Suites already green on the branch, untouched by the attempt — these are the
worked examples for the rest:

    fractal-graph  graph-sequence  graph-set-view
    graph-state    graph-inclusion graph-identity-map

## The residual work-list, measured

| count | class | why it needs judgement |
|---|---|---|
| 33 | `'this' has no IMPLICIT type` | fixture modules whose `use graph_carrier` line carried other imports |
| 27 | `TYPE(graph) → CHARACTER(0)` | `field(...)` call sites still missing the `nentries` argument |
| 21 | `listed_set` undefined | `test/graph-relation/listed_set_fixture.f90` defines a `member_set` subtype that `listed_set_representation` now supersedes — delete, don't port |
| 17 | `sets has no IMPLICIT type` | scopes the declaration pass could not place a map in |
| 16 | syntax error in CALL | multi-line calls the argument rewrite split badly |
| 14+10+6 | wrong actual for `sets` | `, sets` appended to the wrong call by paren-matching |
| 12 | `SELECT TYPE` malformed | `select type` on a domain that is no longer polymorphic — these blocks collapse to straight code |
| 9 | missing `sets` | relation constructions not reached |
| 6 | missing `n_unknown_domain` | `minimizer % attach` call sites |

The two structural ones — `listed_set_fixture` deletion and the `select type`
collapse — are genuine simplifications the cutover *enables*, not damage: a
domain that is always `type(graph)` needs no dynamic dispatch to identify.

## Requirement 6 — the counts, and what they actually mean

| symbol | src | test |
|---|---|---|
| `member_set` | 30 | 396 |
| `counted_set` | 11 | 596 |
| `subset_set` | 17 | 289 |
| `graph_carrier` | 7 | 153 |

The src numbers overstate the situation and should be read with this:

- **Every src reference outside `src/graph_carrier.f90` is a comment.** Verified
  by stripping comment lines: the result is empty. The survivors are
  documentary — `graph_forms` explaining what it stopped extending,
  `graph_label_map` quoting the line that motivated it, `class_graph` naming what
  `carve` replaces.
- **No src module imports `graph_carrier`.** `grep -rn "use graph_carrier" src/`
  is empty. It is already dead production code, compiled only because it is in
  `src/OBJECTS`, and held alive solely by the test suites.

So `subset_set` deletion is blocked by nothing in production. It is blocked by
289 test references, and it is one commit with them.

## Requirements not at issue

**8.** `graph-benchmark` was not run and is not claimed either way; it was
already broken at HEAD (`class_graph_support`, deleted in `a3817e3`) and nothing
here touches it. No new failure is hidden under it, because no suite outside the
six above was run at all.

**9.** `src/fractal_graph.f90` is byte-identical to `57d8c51`.

**1–3.** The production architecture on the branch is untouched — the revert was
`git checkout -- test/` only, and the six passing suites were re-run from a clean
`make -C src clean` + `build.sh` afterwards to confirm it.

## What the next pass should do differently

Not one sweep across 132 files. Suite by suite, smallest first, compiling after
each: the per-suite residual is 10–40 diagnostics, which is reviewable, whereas
1000 at once is not. The three transform passes are sound for the first ~70% of
each file and should be kept; the last 30% is per-file and should be hand-edited,
not pattern-matched.
