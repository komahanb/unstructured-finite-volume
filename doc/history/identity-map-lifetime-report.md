# The identity-map lifetime gate, and the measured shape of the cutover

Branch `tower-graph-as-sets-relations`, from `57d8c51`. `400a2c3` is the GATE.
Sections 4–14 are **measured and reported, not executed**: §8 says *"report the
measured choice before mass call-site edits"*, and the measurement changed the
answer.

## GATE — done

**1. Did set_map or inclusion_map have an internal lifetime dependency? Yes,
both.** Each row borrowed the caller's graph (`type(graph), pointer`) purely to
recognize it later, so each map outlived its own keys. Measured with a probe
that binds inside a scope, lets the binding graphs be deallocated, then asks the
map about **copies carrying the same tokens** — the only probe that separates
recognizing an *identity* from recognizing an *address*. Native (gfortran-15):
every proposition holds — **answered correctly off freed storage**. Under
valgrind: **170 invalid reads from 48 contexts**, in `row_of` in both maps and
through it all six set questions, `included`, `declared_subobject`,
`ambient_of`, and `bind` itself (its duplicate check scans existing rows). The
native run answering *correctly* is the finding. Per §0 no crash was sought; the
pointer established the dependency, the measurement only priced it.

**2. What now stores map keys?** `type(token)`, by value, copied at `bind` /
`include_in`. **3. Does any identity map contain a graph pointer?** No — the word
`pointer` survives in both files only in comments. Both dropped `TARGET` from
their graph dummies, so `test/graph-identity-map/lifetime.f90` (which declares
none) would not compile against the old maps: the compile-time half of the proof.
The law, recorded in tests (§3): *an identity map owns its keys by value; it
borrows no graph object merely to recognize it.*

**`ambient_of` is gone, not re-homed.** Its only caller wrote
`host = m % ambient_of(s); host % same_as(a)` — the identity predicate in two
steps. `declared_into(part, ambient)` answers it in one, the transitive walk
steps in tokens, and no registry was invented (§2).

Suites **32/33**. `test/graph-benchmark` was already broken at HEAD
(`class_graph_support`, deleted in `a3817e3`) — pre-existing, untouched.
**13. `src/fractal_graph.f90` is byte-identical to `57d8c51`.**

## §5 audit — what field-domain callers require

Every held `member_set` in `src/`, bucketed by the question asked of it. Tests
only price the edit; production determines the API.

| role (src/) | I ident | S size | M has | E enum | L local_idx | O subobj | N name | tot |
|---|---|---|---|---|---|---|---|---|
| **field-domain** | 13 | 6 | 1 | 1 | 2 | 4 | 2 | 29 |
| relation-domain | 1 | 7 | 2 | 5 | 4 | 0 | 0 | 19 |
| slot-carrier | 1 | 1 | 1 | 0 | 1 | 0 | 0 | 4 |
| other/constructed | 13 | 6 | 4 | 4 | 1 | 2 | 4 | 34 |

**4. What did field-domain callers require? All seven buckets** — not
identity+size. The non-identity ones concentrate in `class_graph_assembler`
(E, O, N, S) and `class_graph_partitioner` (M, L, O, N): the pair §17 forbids
touching. Two consequences outrank the storage question:

- `dom % name()` (assembler:315, partitioner:567) names a derived subset.
  **`type(graph)` has no name** — the kernel type is `branch(2)` plus a private
  token. Bucket N has *no home* in the §6 law. That is a gap in the law, not an
  implementation detail.
- `is_subobject_of` (4 sites) becomes `declared_subobject(dom, X, inclusion_map)`,
  and every E/L/M site needs a **set_map** in scope. Partition and assembly hold
  only `class(graph)` and `class(graph_field)`; neither map exists there. The one
  in-scope object that could supply them is the graph — i.e. the field-domain
  cutover forces a set_map onto `stored_graph`, the graph-host seam
  `doc/REVERSE-ARCHITECTURE-REVIEW.md` rejected on production evidence.
  **That call is yours, not mine.**

## §8 — field domain storage, measured

**6. Does a field duplicate O(N_extent) domain storage? Yes, today.** 40 fields
on one 200 000-member domain, RSS-measured:

| domain | batch cost | values alone | excess |
|---|---|---|---|
| counted | 62.8 MB | 61.0 MB | ~0 — O(1) domain |
| subset (listed) | 91.5 MB | 61.0 MB | **28.7 MB** |

Predicted for K copies of an N-integer roll: 30.5 MB; the measurement lands on
it. Per §8 that is not evidence the final architecture should. **B** and **C**
reproduce the 28.7 MB exactly, C with row storage on top; **A** and **D** are the
O(1) candidates.

**5. Measured choice: A for what a field retains, D for the rest.** A field
retains `type(graph) :: domain` + an integer `num_entries` snapshot — enough for
`num_entries`, the five value-vector width checks, and identity, standalone and
map-free. M/E/L are answered by the caller through a `set_map`, O through the
`inclusion_map`. **The snapshot changes no behaviour**, and the proof is one
line: `class_graph_field.f90:161` does `allocate(this % on, source=on)`, so a
field has held a *copy* of its domain since construction and `num_entries` is
already frozen against later mutation of the caller's set. The one mutable
member_set is a `form` (via `restrict`, sole caller `class_form_pruner`), and
mutating it never reached a field built on it. A alone is insufficient, and D is
not free — D is what puts a map into partition/assembly scope, i.e. the seam
above.

## §11 — the CSR hot path, measured

**9. How does CSR get local_index/member without semantic traversal?** Today it
need not: `csr_relation` holds slot-wrapped carrier copies and calls
`carrier % local_index(member)` **inside `image_view`, `preimage_view` and `has`**
(`graph_binary_relation.f90:410, 428, 453`) — per query, not only at build
(`create_csr:316-317`). A `type(graph)` signature would therefore put a map scan
plus a dispatch in the hot path, which §11 forbids: CSR must keep its own
numbering, and that cost is why §9 and §11 are one commit. Pre-existing and worth
stating: on a *subset* carrier `local_index` is an honest linear scan, so
`image()` is already **O(N_extent), not O(degree)**, there today.

## Not executed, and why

**7.** Signature is still `type(slot)` over `class(member_set)`. **8.** slot
remains. **10.** `form extends subset_set` still. **11.** subset_set remains.
**12. Legacy surface:** `subset_set` 52 src / 289 test, `member_set` 163 / 396,
`counted_set` 43 / 596, `slot` 54 / 184, `graph_carrier` 38 / 153 — nothing
opportunistically widened. Beyond volume, the blocker is `stored_relation`
validating through `carrier % has`: with graph signatures it needs a set_map, so
**every relation construction site in the repository changes**. That is what
makes A/B/C/D one atomic commit rather than four.

**14. Next atomic boundary.** Same membership — field domain + relation signature
+ form inheritance, one commit — with two prerequisites that measurement has now
made explicit and that are design calls, not implementation details:

1. **Where does a set_map live**, so partition, assembly, relation construction
   and CSR can reach one? The only in-scope candidate is the graph, and that seam
   was rejected once on production evidence. Answer before any edit.
2. **Bucket N has no home.** Decide whether a domain name is the field's, the
   call site's, or a role the §6 law must admit.
