# Test migration: complete — 31 of 32 suites green, `graph_carrier` deleted

Branch `tower-domain-cutover`. Every suite is off the retired carrier ontology
and onto set graph identity + representation + the three maps. The §9 gate is
**open and executed**: `src/graph_carrier.f90`, its `src/OBJECTS` entry and the
`test/graph-carrier` suite are gone.

The denominator moved from 33 to 32 by that deletion, not by drift.

## Where it stands

    green                31 / 32
    failing              graph-benchmark, at its baseline
    src/fractal_graph.f90 byte-identical to 57d8c51

`graph-benchmark` fails on `class_graph_support.mod`, deleted in `a3817e3`,
exactly as it did before this work began. Its four carrier lines were migrated
*before* the module was deleted, so the failure kept its shape rather than
becoming a dangling `use` — but that migration has therefore **never been
through the compiler**: the file dies on `class_graph_support` at line 33,
above every line that changed.

## §9 — the gate, now open

    use graph_carrier            src 0   test 0
    member_set                   src 0   test 0
    counted_set                  src 0   test 0
    subset_set                   src 0   test 0

Counts are of live code and of comments alike: `grep -rn graph_carrier src/`
returns nothing at all. The abstract interface that used to be called
`graph_carrier_interface` in `graph_grammar` — an identifier, never an import —
is now `set_graph_interface`.

Comments that name the retired types *historically* remain, and are still true:
what `graph_forms` once extended, what the old `counted_set` constructor used to
mint, what `graph-inclusion` restated. Comments that described the code **as it
stands** were restated in `586fdbc`.

## Where each carrier law went

| carrier law | successor |
|---|---|
| counted contract: size, membership, position, 1..n | `graph-set-view` §2 |
| structural identity: equal extensions, two sets | `graph-set-view` §1 |
| subsets, embedding order, declared provenance | `graph-inclusion` §1–3 |
| the name a set wears | `graph-identity-map/labels` |
| refusal `twice` | `fractal-graph`: identity is assigned once |
| refusal `unhosted` | `graph-inclusion`: `unsigned` |
| a graph mints its carriers once | `graph-contract` check 21 |

Two laws have **no** successor. Both are rulings, not omissions:

- **refusal `outsider`** — a subset member from beyond the ambient is no longer
  refused. Inclusion is *declared*, and `graph-inclusion` §2 rules on it
  directly: containment does not imply an inclusion edge, so its absence does
  not forbid one.
- **the carrier's own name** — a graph binds no label. It mints identity and
  counts; naming belongs to the label map.

The one law that moved rather than died is `graph-contract` check 21: a concrete
graph mints its two carrier identities once, hands out the same identity at
every asking, its two sides are two identities, and a twin built from identical
numbers owns its own. That is a law about `class_graph`, not about the retired
module. It is restated so the count comes from the graph and the size only from
a map that has been told.

## The two production repairs, pinned

Both were found by running, and neither had a test naming it. Both are pinned in
`graph-multigrid`:

- **`attach` is re-enterable.** Newton calls it once per iteration and it
  declared the solver's number domain each time; a graph signs once, so the
  second iteration died. The old `counted_set` constructor minted a *fresh*
  domain per attach, so that reading is kept — reset to an unsigned graph, then
  declare.
- **A count component needs a default.** `jacobi()` and `gmres()` name no
  component, so `n_unknown_domain` / `n_residual_domain` must have values
  without being named.

## The import gates

Each tower's `check_imports.sh` is its layering law as a script. All seven now
grant the four set capabilities separately rather than the single carrier grant:

    fractal_graph             identity — granted broadly, because WHICH set is
                              a question everything asks
    graph_set_representation  only where a representation is CONSTRUCTED
    graph_set_map             only where one is BOUND or QUERIED
    graph_inclusion_map       only where provenance is ASSERTED, or carve called
    graph_label_map           only where a label is READ, or carve called

Three gates additionally `refuses ... graph_carrier` by name. Those assertions
stay: they are the towers' own record that the module is retired, and they go on
passing.

## §8 — partition laws, unchanged

`sum_k A_k(P_k(D)) = D` at nparts 2 and 3, on vertex **and** edge data; every
cell owned exactly once; identity at one part and **not** at two. Nothing is
called an inverse.
