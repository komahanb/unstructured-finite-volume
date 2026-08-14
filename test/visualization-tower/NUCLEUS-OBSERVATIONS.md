# Nucleus observations — Visualization Tower

The live evidence ledger of this orbital client, per the discipline of
[FRACTAL-GRAPH-ARCHITECTURE.md](../../FRACTAL-GRAPH-ARCHITECTURE.md):

\[
\boxed{\text{observation}\neq\text{immediate refactor}}
\]

*Local necessity* means necessity at the object's own radius; *global
necessity* means necessity at larger contextual radius. Cross-tower
recurrence is its own field and is never used as proof of global
necessity.

Every entry names the **LEVEL** that owns its evidence and the review
**gate** at which that level was reviewed. Levels are the architecture;
gates are checkpoints. **The ten levels of this tower are ONE client,
not ten.** Six of them are built; the tower is **not sealed**, and
Review Gate B is **not reached** - Gate B reviews Levels 5, 6 and 7
together, and two of those do not exist.

---

## A standing caution for this ledger

**VIZ-1 … VIZ-10 are Gate-A observations, and not one of them is
evidence for any seam of**
[the reverse architecture review](../../doc/REVERSE-ARCHITECTURE-REVIEW.md).

Seam **A2** — *operations should carry their domain rather than ask a
graph for one* — cannot have been exercised here, because **Gate A
attached no operation to anything.** No `graph_operation` was
constructed, none was applied, and the import gate refuses every
module that owns one at every level. The seam-A2 count stands where
the Time Integration Tower left it, at **three**, and nothing below
moves it.

Seam **B** — bidirectional/rectangular linearization — is likewise
untouched. This tower is *full* of rectangular structure, and none of
it is a linearization: `D_k^T` reverses the orientation of possible
dependence between two carriers, with no coefficient to transpose and
no vector to apply. **Structural transpose is not numerical adjoint**,
and calling the resemblance evidence would be the single easiest error
this tower could make.

**Gate A's production change is NONE**, which is itself the headline
observation and is recorded as VIZ-10. **Level 5's production change
is also NONE**, recorded in VIZ-11.

VIZ-11 through VIZ-15 belong to Level 5 and are likewise seam-free:
attaching a field to an occurrence carrier is not an operation taking
a domain from a graph, and a coefficient is not a linearization.

---

## OBSERVATION VIZ-1

```text
tower:                Visualization
level:                0-1  (carrier, relation)
review gate:          A
contextual radius:    1 (typed relations over declared carriers)

symptom / fact:       A STRUCTURAL DEPENDENCY NATURALLY HAS SEPARATE
                      INPUT AND OUTPUT CARRIERS.

                        D1 <= X0 x X1,  X0 /= X1
                        D2 <= X1 x X2,  X1 /= X2   (both size 3)
                        D3 <= X2 x X3

                      Every mesh client the repository has built reads
                      relations whose two ends land in ONE vertex
                      carrier. An operator chain does not: A1's input
                      domain and output domain are two different
                      worlds, and the tower declares them so before
                      any relation exists.

                      |X1| = |X2| = 3 and X1 is not X2. The specimen
                      plants that coincidence deliberately, so no
                      level can pass by comparing cardinalities.

exact caller:         test/visualization-tower/level-1-relation/test.f90
                      check_the_signatures_are_rectangular

evidence:             T1 % target() and H1 % target() are not same_as,
                      on all three operators. PASS.

confidence:           high - three independent instances, checked by
                      structural identity rather than by size.

action:               OBSERVE. The nucleus already permits it: a
                      relation signature holds one carrier per slot
                      and never asked for them to agree.
```

---

## OBSERVATION VIZ-2

```text
tower:                Visualization
level:                0  (carrier)
review gate:          A
contextual radius:    0 (declared sets, nothing between them)

symptom / fact:       THE OCCURRENCE CARRIER GIVES DEPENDENCY
                      OCCURRENCES FIRST-CLASS IDENTITY BEFORE ANY
                      COEFFICIENT EXISTS.

                        E1 = { e11 ... e15 }
                        E2 = { e21 ... e24 }
                        E3 = { e31 e32 e33 }

                      e11 is a member of a declared domain, not a
                      pair, not a nonzero, not an index into a table.
                      It can therefore be related, composed away, and
                      - at Level 5 - given a value. VIZ-12 records
                      that E_k turned out to be the ONLY carrier that
                      could seat the coefficient.

                      This is a MODELLING CHOICE the nucleus made
                      possible rather than a capability it supplied.
                      The same choice is what the mesh path makes when
                      it declares an edge carrier instead of storing
                      pairs.

exact caller:         test/visualization-tower/common/
                      visualization_carriers_fixture.f90
                      structural_carriers

evidence:             |E1|=5, |E2|=4, |E3|=3, all distinct from each
                      other and from the four state carriers. PASS at
                      Level 0.

confidence:           medium-high - one specimen, but the choice is
                      the reason Level 5's zero-coefficient experiment
                      was askable at all.

action:               OBSERVE. Nothing follows for production. E_k did
                      survive to Level 5 and did carry w_k; see
                      VIZ-12.
```

---

## OBSERVATION VIZ-3

```text
tower:                Visualization
level:                2  (relation algebra)
review gate:          A
contextual radius:    2 (composition over typed relations)

symptom / fact:       D_k IS DERIVED FROM PRIMITIVE INCIDENCE, NOT
                      INDEPENDENTLY STORED.

                        D_k = H_k o T_k^T
                            = compose_binary(T_k^T, H_k)

                      Nothing in the tower writes down a -> p. The
                      dependency is a consequence of two facts about
                      the occurrence e11, and the relation algebra
                      draws it.

                      This is a discipline the nucleus permitted but
                      did not enforce. Writing D1 out by hand beside
                      T1 and H1 would have compiled, passed, and
                      tested nothing about composition.

exact caller:         test/visualization-tower/common/
                      visualization_algebra_fixture.f90
                      derive_dependency

evidence:             Level 1 asserts that NO primitive relation has
                      signature X_(k-1) x X_k. Level 3 asserts that
                      NO graph-owned relation does either. Level 2
                      then derives all five dependencies and checks
                      their exact extensions. PASS.

confidence:           high.

action:               OBSERVE.
```

---

## OBSERVATION VIZ-4

```text
tower:                Visualization
level:                2  (relation algebra)
review gate:          A
contextual radius:    2 (composition along a three-operator chain)

symptom / fact:       STRUCTURAL COMPOSITION PRODUCES THE POSSIBLE-
                      DEPENDENCY SKELETON OF THE COMPOSED OPERATOR.

                        D2:1 = D2 o D1        : X0 -> X2
                        D3:1 = D3 o D2 o D1   : X0 -> X3

                      D3:1 is what a person would call the sparsity of
                      A3 A2 A1 - obtained with no A, no entries, and
                      no multiplication. It is BOOLEAN/RELATIONAL
                      composition, and the distinction is not
                      cosmetic: numerical matrix multiplication can
                      cancel, and relational composition cannot.

                      Associativity was checked rather than assumed:
                      (D3 o D2) o D1 and D3 o (D2 o D1) agree
                      extensionally.

exact caller:         test/visualization-tower/level-2-relation-algebra/
                      test.f90  check_the_full_composition

evidence:             |D3:1| = 6, exact extension
                      { a->m, b->m, c->m, b->n, c->n, d->n }. Both
                      bracketings agree. Both intermediate carriers
                      X1 and X2 appear in neither slot of the result.
                      PASS.

confidence:           high.

action:               OBSERVE. Note for a future Level 6 that a
                      numerically cancelling A would still have this
                      structural skeleton, and that the difference is
                      the point.
```

---

## OBSERVATION VIZ-5

```text
tower:                Visualization
level:                2  (relation algebra)
review gate:          A
contextual radius:    2 (composition through a shared middle carrier)

symptom / fact:       MULTIPLE INTERMEDIATE PATHS COLLAPSE TO ONE
                      ENDPOINT TUPLE UNDER SET SEMANTICS.

                      Seven walks run X0 -> X1 -> X2. D2:1 holds six
                      tuples. Both

                        b -> p -> u        and        b -> q -> u

                      witness the dependency b -> u, and it is held
                      ONCE.

                      The nucleus enforces this in the constructor -
                      csr_relation collapses duplicate tuples to first
                      appearance - so the level checks that it did
                      rather than doing it itself.

                      The dual fact matters as much: a->v has ZERO
                      witnesses and is absent. Composition adds
                      nothing it cannot walk to.

exact caller:         test/visualization-tower/level-2-relation-algebra/
                      test.f90  check_the_witness_collapse

evidence:             walks(D1, D2) = 7, |D2:1| = 6, and the tuple
                      table lists b->u exactly once. PASS.

confidence:           high.

action:               OBSERVE. This is the acceptance law that makes
                      D2:1 the structure of a composed operator rather
                      than an accounting of paths through it.
```

---

## OBSERVATION VIZ-6

```text
tower:                Visualization
level:                2  (relation algebra)
review gate:          A
contextual radius:    2 (transpose of typed relations)

symptom / fact:       TRANSPOSE GIVES REVERSE DEPENDENCY ORIENTATION
                      WITHOUT ANY NUMERICAL ADJOINT MACHINERY.

                        D1^T : X1 -> X0
                        D2^T : X2 -> X1
                        D3^T : X3 -> X2

                      The nucleus supplies transpose_of as an O(1)
                      VIEW. The tower materializes copies of those
                      views so they can be owned and passed, and
                      checks that materializing added ownership and
                      nothing else.

                      THE CAUTION THIS ENTRY EXISTS FOR: none of this
                      is an adjoint. There is no coefficient to
                      transpose, no vector to apply, no inner product
                      under which anything is dual, and no sensitivity
                      being propagated. D^T answers exactly one
                      question - which members could have reached this
                      one - and the resemblance to reverse mode is a
                      resemblance.

exact caller:         test/visualization-tower/level-2-relation-algebra/
                      test.f90  check_the_reverse_dependencies

evidence:             Every tuple x->y of D1 appears in D1^T as y->x,
                      on all five; counts unchanged; the materialized
                      copy is same_extension with the live
                      transpose_of view; (D1^T)^T = D1 extensionally.
                      PASS.

confidence:           high.

action:               OBSERVE. Explicitly: SEAM B GETS NO VOTE FROM
                      THIS TOWER. See the standing caution above.
```

---

## OBSERVATION VIZ-7

```text
tower:                Visualization
level:                2-3  (relation algebra, relational graph)
review gate:          A
contextual radius:    2 (the composition/transpose law)

symptom / fact:       (D3 o D2 o D1)^T AND D1^T o D2^T o D3^T AGREE
                      EXTENSIONALLY.

                      The law holds in the nucleus with no help:
                      compose_binary and transpose_of were written for
                      a mesh, and neither knew a chain of operators
                      would ask this of them.

                      Compared as EXTENSIONS - domains slot for slot
                      by same_as, cardinality, and membership both
                      ways round - never by tuple enumeration order,
                      because a relation is a set and how its tuples
                      were handed in is no part of what it is.

                      Orientation is part of the answer: D3:1 is NOT
                      equal to D3:1^T, and same_extension says so
                      because its first test is domain identity.

exact caller:         test/visualization-tower/level-2-relation-algebra/
                      test.f90  check_the_transpose_composition_law
                      and level-3-graph/test.f90
                      check_the_chain_is_recoverable

evidence:             same_extension(D_rev, D3:1^T) is true; both run
                      X3 -> X0; six tuples each; and the same law
                      closes again when every relation is re-derived
                      from GRAPH-OWNED primitives alone. PASS.

confidence:           high - proved twice, from two sources of the
                      same primitives.

action:               OBSERVE.
```

---

## OBSERVATION VIZ-8

```text
tower:                Visualization
level:                4  (structural interpretation)
review gate:          A
contextual radius:    3 (representation over graph-owned structure)

symptom / fact:       A STRUCTURAL RELATION ADMITS MULTIPLE
                      REPRESENTATIONS WITHOUT BECOMING THOSE
                      REPRESENTATIONS.

                      One relation, three renderings:

                        chain line      X0 --D1--> X1 --D2--> ...
                        sparsity grid   rows codomain, cols domain
                        listing         a -> p / b -> p q / ...

                      and the grid and the listing are required to
                      agree cell for cell, because they ask the
                      identical question in the identical order.

                      The renderer knows nothing about the specimen.
                      Order comes from member(i) and local_index;
                      content from relation % has; names from the
                      objects that carry them. A hostile carrier
                      declaring { 30, 10, 20 } renders 30 10 20.

                      THE FLOW IS ONE-DIRECTIONAL. The oracle pictures
                      live in the Level-4 test and are unreachable
                      from the renderer; a representation that could
                      change what is true would not be a
                      representation.

exact caller:         test/visualization-tower/common/
                      structural_renderer_fixture.f90, exercised by
                      level-4-graph-calculus/test.f90

evidence:             Five representations generated and matched
                      against hand-worked oracles; every cell of six
                      grids re-checked against relation % has; the
                      reverse skeleton drawn twice from two
                      independently composed relations and required to
                      agree grid for grid while differing in its name
                      line. PASS.

confidence:           high.

action:               OBSERVE. Do NOT promote the renderer to
                      production. It is a test fixture, and its whole
                      evidential value is that it may use nothing but
                      the nucleus.
```

---

## OBSERVATION VIZ-9

```text
tower:                Visualization
level:                4  (structural interpretation)
review gate:          A
contextual radius:    3 (the ordinary-graph profile's schema)

symptom / fact:       AN ORDINARY-GRAPH INTERPRETATION IS NOT REQUIRED
                      AT THIS RADIUS, AND WOULD BE INAPPROPRIATE IF
                      FORCED.

                      Both profile readings demand a SINGLE vertex
                      carrier:

                        ordinary_graph_view      T <= E x V
                                                 H <= E x V
                          refuses with 'the head relation must share
                          the tail''s domains'

                        directed_adjacency_view  A <= V x V
                          refuses with 'a directed adjacency runs over
                          one domain'

                      The specimen cannot satisfy either: T1 lands in
                      X0 and H1 in X1; D1's source is X0 and its
                      target X1.

                      Satisfying either would mean manufacturing
                      V = X0 u X1 u X2 u X3, and it is worth being
                      exact about the cost.

                      DIRECTION IS NOT WHAT WOULD BE LOST. An ordinary
                      directed graph preserves ordered-pair direction
                      perfectly well - directed_adjacency_view
                      documents that the tuple order of a same-domain
                      relation IS the direction, and (a,p) would still
                      be distinguishable from (p,a) inside a union.

                      WHAT WOULD BE LOST is the intrinsic distinct
                      TYPED source and codomain:

                        D : X_i -> X_j    against    D^T : X_j -> X_i

                      Under a union both become relations over one V.
                      The two declared carriers that made them
                      different KINDS of object are gone, recoverable
                      only from an offset convention the mathematics
                      never stated - and the first thing
                      same_extension tests, domain identity, has
                      nothing left to compare.

                      THIS IS NOT A DEFECT IN THE PROFILE. It is a
                      specialization the profile documents, correct
                      for the mesh path it was written for.

exact caller:         test/visualization-tower/level-4-graph-calculus/
                      test.f90  check_the_ordinary_graph_question,
                      against src/graph_profile.f90 create_view and
                      create_adjacency_view

evidence:             The premises are asserted executably; the
                      contracts were inspected rather than provoked,
                      per the brief. Level 3 asserts the graph holds
                      no twelve-member union carrier. The import gate
                      refuses graph_profile at EVERY level, so the
                      pictures cannot have leaned on one. PASS.

confidence:           high - one specimen, but the conclusion follows
                      from the profile's own documented schema rather
                      than from this specimen's shape.

action:               OBSERVE. Verdict: ordinary graph interpretation
                      is INAPPROPRIATE for rectangular typed
                      dependency, and unnecessary for its
                      visualization. No production change follows.
```

---

## OBSERVATION VIZ-10

```text
tower:                Visualization
level:                0-4  (all built levels)
review gate:          A
contextual radius:    3 (the whole Gate-A client)

symptom / fact:       VISUALIZATION ASKED FOR NO NEW MATHEMATICS. IT
                      ASKED ONLY FOR AN EXTERNAL RENDERER.

                      Five representations of a three-operator chain -
                      forward, reverse, three individual sparsities,
                      the composed skeleton and its transpose - were
                      produced, and NO RED OCCURRED AT ANY LEVEL.

                      git diff 44ae3da -- src/  is EMPTY.

                      AND THE RESULT IS SHARPER THAN "SUFFICIENT".
                      The two jobs used different amounts of nucleus.

                      DERIVING the structure needed the binary
                      specialization and the algebra:

                        binary_relation     source, target,
                                            transpose_of
                        relation algebra    compose_binary
                        relational_graph    ownership and closure

                      INTERPRETING it needed strictly less - the ROOT
                      relation contract and the carriers, and nothing
                      whatever below them:

                        relation            arity, domain(k), has,
                                            name
                        member_set          size, member, local_index,
                                            name

                      structural_renderer_fixture names no binary
                      relation, calls no source(), no target(), no
                      image_view and no transpose_of. It is handed
                      transposed and composed relations and cannot
                      tell them from any other. So RENDERING A
                      RECTANGULAR TYPED DEPENDENCY REQUIRED LESS OF
                      THE NUCLEUS THAN BUILDING ONE DID.

                      The absence of numbers is part of the claim, not
                      an omission - and it is ENFORCED, not asserted:
                      check_imports.sh refuses any real or complex
                      declaration and any literal carrying a decimal
                      point or an exponent, in every Gate-A source,
                      with its own selftest and a planted-violation
                      check. The structural picture preceded the
                      coefficients entirely.

exact caller:         the whole of test/visualization-tower/

evidence:             L0-L4 PASS. src diff empty. No production
                      abstraction added: no graph_visualization, no
                      class_graph_visualization, no graph_renderer, no
                      visualize()/print()/dependencies() on
                      graph_operation.

confidence:           high for the built radius; UNKNOWN beyond it.
                      Level 6 has not yet asked production what
                      dependencies() thinks a domain is, and that is
                      the level most likely to find pressure.

action:               OBSERVE, and BUILD NOTHING. Success here means
                      the nucleus was already sufficient - which is a
                      reason to add no abstraction, not a reason to
                      add one because the tower went well.
```

---

## OBSERVATION VIZ-11

```text
tower:                Visualization
level:                5  (field calculus)
review gate:          not yet reviewed - Gate B needs L5, L6, L7
contextual radius:    3 (fields over an occurrence carrier)

symptom / fact:       A RELATION AND A FIELD MAY SHARE AN
                      INTERPRETIVE PICTURE WITHOUT EITHER DEFINING
                      THE OTHER.

                      Three real fields were attached:

                        w1 : E1 -> R = [ 2, -1,  0,  3,  4 ]
                        w2 : E2 -> R = [ 1,  5, -2,  2 ]
                        w3 : E3 -> R = [ 3, -1,  4 ]

                      and afterwards D1, D2, D3, D2:1 and D3:1 held
                      exactly the tuples Gate A derived, and the
                      Level-4 pictures were the same strings
                      CHARACTER FOR CHARACTER - captured before the
                      fields existed and compared after.

                      The independence runs both ways, and a second
                      field on the same E1 shows the other direction:
                      w1_alt = [9, 8, 7, 6, 5] draws a different
                      coefficient picture over an identical structural
                      one. So D1 does not determine w1, and w1 does
                      not determine D1.

                      NO PRODUCTION CHANGE WAS NEEDED. The existing
                      field machinery already seats a coefficient on
                      an occurrence carrier, already takes its entry
                      count from the domain, and already answers
                      domain questions by identity.

exact caller:         test/visualization-tower/level-5-field-calculus/
                      test.f90
                      check_the_structural_pictures_are_unchanged,
                      check_the_derivation_is_unchanged,
                      check_value_and_structure_are_independent

evidence:             L5 PASS. git diff 29c0ccd -- src/ is EMPTY.

confidence:           high.

action:               OBSERVE.
```

---

## OBSERVATION VIZ-12

```text
tower:                Visualization
level:                5  (field calculus)
review gate:          not yet reviewed
contextual radius:    3 (where a coefficient can sit)

symptom / fact:       DEPENDENCY COEFFICIENTS NATURALLY INHABIT E_k,
                      NOT THE SOURCE CARRIER AND NOT THE TARGET
                      CARRIER.

                      The reason is forced by the specimen rather
                      than chosen:

                        e12 : b -> p    w1(e12) = -1
                        e13 : b -> q    w1(e13) =  0

                      Two occurrences read the same member b and
                      carry different numbers. A field on X0 could
                      hold only one of them. A field on X1 could hold
                      only one of them. E1 is where the two ends
                      meet, and it is the only one of the three
                      carriers that can seat both.

                      This is what Level 0's occurrence carrier was
                      FOR, five levels before anything needed a
                      value - and the justification only became
                      visible here.

exact caller:         test/visualization-tower/common/
                      visualization_values_fixture.f90, and
                      level-5-field-calculus/test.f90
                      check_the_domains_by_identity

evidence:             domain(w1) same_as E1, and NOT same_as X0 or
                      X1. Likewise w2/E2 and w3/E3. PASS.

confidence:           high.

action:               OBSERVE. Note for a future Level 6 that a
                      production operator's coefficients must also
                      have somewhere of this kind to live, and that
                      "somewhere" is the question Level 6 will ask.
```

---

## OBSERVATION VIZ-13

```text
tower:                Visualization
level:                5  (field calculus)
review gate:          not yet reviewed
contextual radius:    3 (a value on an existing structural seat)

symptom / fact:       NUMERICAL ZERO IS NOT STRUCTURAL ABSENCE.

                        e13 : b -> q        T1(e13) = b, H1(e13) = q
                        w1(e13) = 0
                        b -> q IS IN D1

                      D1 still holds five tuples, not four. The
                      structural picture still writes # at (b,q),
                      because it never asked what the value was.

                      And the two facts are DIFFERENT facts in the
                      representation:

                        (b,q)  ->  0     present, carrying zero
                        (a,q)  ->  .     absent, nothing to carry

                      A picture that wrote 0 for both would have
                      discarded the distinction.

                      THE FORBIDDEN INFERENCE, stated so it can be
                      checked for later:

                        if (coefficient /= 0) then present

                      appears nowhere. Presence comes from
                      relation % has and from nothing else; the value
                      is consulted only to decide how the number
                      LOOKS.

exact caller:         test/visualization-tower/level-5-field-calculus/
                      test.f90  check_the_zero_witness,
                      check_zero_is_not_absence

evidence:             |w1(e13)| < 1e-14; d1 % has([b,q]) true;
                      |D1| = 5; glyph_at(D1,b,q) = '#'; the rendered
                      q row carries both a 0 and a dot, in different
                      columns. PASS.

confidence:           high - this is the level's load-bearing
                      assertion and it is checked five ways.

action:               OBSERVE.
```

---

## OBSERVATION VIZ-14

```text
tower:                Visualization
level:                5  (field calculus)
review gate:          not yet reviewed
contextual radius:    3 (representation fusing two sources)

symptom / fact:       A VALUED SPARSITY PICTURE COMBINES TWO
                      INDEPENDENT QUERIES, IN ONE ORDER AND NEVER THE
                      OTHER.

                        structure   glyph_at(D, x, y), which is
                                    relation % has - the Level-4
                                    renderer's own primitive, CALLED
                                    here rather than reimplemented

                        value       the unique occurrence e with
                                    T(e)=x and H(e)=y, then w(e)

                      The second may never answer the first. The
                      valued renderer asks structure, and only if the
                      answer is yes does it look for a coefficient.

                      THE DEPENDENCY DIRECTION IS THE OTHER HALF OF
                      THE OBSERVATION:

                        L4 structural renderer
                                ^
                        L5 valued renderer

                      Level 4 was not made field-aware. It imports
                      nothing of Level 5, still cannot hold a real
                      number - the numberless law refuses it one - and
                      remains able to render mathematics that has no
                      coefficients at all.

                      Nothing in the Level-4 renderer was modified.
                      The valued picture takes its page width and its
                      stub POSITION from a structural picture of the
                      same relation, so no measurement is duplicated
                      as a literal and the two grids align row for
                      row.

exact caller:         test/visualization-tower/common/
                      valued_renderer_fixture.f90, exercised by
                      level-5-field-calculus/test.f90
                      check_the_two_views_agree_on_presence

evidence:             On every cell of D1, D2 and D3 the structural
                      view, the coefficient view and the relation
                      itself agree on presence - read by tokenizing
                      both rendered grids, so the two layouts may
                      differ freely. Still agree when the
                      coefficients are replaced wholesale. PASS.

confidence:           high.

action:               OBSERVE. Do NOT promote either renderer to
                      production.
```

---

## OBSERVATION VIZ-15

```text
tower:                Visualization
level:                5  (field calculus)
review gate:          not yet reviewed
contextual radius:    3 (field domain identity)

symptom / fact:       EQUAL CARDINALITY DOES NOT MAKE ONE CARRIER A
                      VALID COEFFICIENT DOMAIN FOR ANOTHER.

                      The specimen supplies the hostile case for
                      free:

                        |X0| = |E2| = 4

                      A four-valued real field on X0 was built. It is
                      a perfectly valid field - the nucleus accepted
                      it, it has four entries, one component, real
                      kind - and it is REFUSED as E2's coefficients,
                      because its declared domain is not E2.

                      NO NEW VALIDATION MECHANISM WAS INVENTED. The
                      check asks the field for its domain and the
                      domain for its identity, which is what
                      graph_carrier and class_graph_field have always
                      offered. Level 5 exposes those semantics; it
                      does not add to them.

                      A size-based check would have accepted the
                      decoy, which is exactly why the specimen was
                      built with the collision in it at Level 0.

exact caller:         test/visualization-tower/common/
                      valued_renderer_fixture.f90 coefficients_fit,
                      and level-5-field-calculus/test.f90
                      check_the_domains_by_identity

evidence:             coefficients_fit(decoy_on_X0, E2) is FALSE
                      while decoy % num_entries() = 4 = |E2|. Also
                      false for (w1,E2), (w2,E1), (decoy,E1). PASS.

confidence:           high.

action:               OBSERVE.
```

---

## Seam accounting through Level 5

| Seam | Before this tower | After Gate A | Why |
|---|---|---|---|
| **A1** graph host as conduit | CLOSED | CLOSED | nothing here reopens it; no operation was hosted anywhere |
| **A2** operations take domain from graph | 3 towers | **3 towers — unchanged** | Gate A attached no operation to anything, and Level 5 attached a FIELD to a CARRIER — no `graph_operation` exists anywhere in this tower, so the seam could not have been exercised |
| **A3** relational_graph ownership | KEEP | KEEP | one more successful typed-ownership pattern (7 carriers, 6 relations, full signature closure); no production change follows |
| **B** bidirectional/rectangular linearization | 2 towers | **2 towers — ZERO new votes** | this tower is full of rectangular structure and none of it is a linearization; structural transpose is not numerical adjoint, and Level 5 built no `w^T`, applied no `A^T v`, and composed no coefficients |

**This tower is ONE client, not six.** Level 5 produced no RED. If
Levels 6–9 are later built and do produce REDs, they will still be
**one** tower's vote.

---

## Frontier

Level 5 asked *"does a numerical zero erase structural dependency?"*
and answered **no**. `structure != value` is now a test rather than a
slogan.

```
NEXT FRONTIER = Level 6 — does a production operator expose
                          this same structural skeleton?
```

The question waiting there: how does the structural relation `D`
correspond to `discretization_operator % dependencies()`? Nothing is
pre-decided, `dependencies()` has not been generalized, and Level 6
must RED against the persistent typed chain before production is
touched.
