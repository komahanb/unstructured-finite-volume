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
not ten.** All ten are built and the tower is SEALED. Review Gates A,
B and C are behind it.

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

**VIZ-25 through VIZ-30 belong to Level 7, and VIZ-28 is this
ledger's first RED - captured, then RESOLVED.** It was a RED against
an executable mathematical contract - the diagonal of an unchanged
matrix - and not against a comment. The repair separated the
minimizer's execution context from its dependent-variable coupling;
it did NOT reach into the action for structure, for the reason
VIZ-29 records.

**A CORRECTION TO VIZ-20 AND VIZ-24.** Those entries read the two
concretes' answers as different structural AXES and inferred that the
family verb had no invariant meaning. That inference was wrong, and
the error was in step_dependencies rather than in the family: it
answered the succession 1->2->3, which is a chronology, where the
contract owes a stencil. Corrected, BDF2 answers the fan-in 1->3,
2->3, 3->3, and both citizens answer ONE thing - the stencil on the
axis the concrete type represents. See VIZ-31.

**VIZ-16 through VIZ-24 belong to Level 6, and they are the first
entries in this ledger about PRODUCTION.** They are still seam-free,
and for a reason worth stating: seam A2 is about an operation taking
its DOMAIN from a graph, and Level 6 never applied an operation to
anything. It constructed two and asked them what their dependency
pattern is. **Level 6's production change is NONE**, and no RED
occurred - see VIZ-23.

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

                        directed_incidence_view      T <= E x V
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
                      test.f90  check_the_directed_incidence_question,
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

                        w1 : E1 -> reals = [ 2, -1,  0,  3,  4 ]
                        w2 : E2 -> reals = [ 1,  5, -2,  2 ]
                        w3 : E3 -> reals = [ 3, -1,  4 ]

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
                      graph_carrier and field_stored have always
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

## OBSERVATION VIZ-16

```text
tower:                Visualization
level:                6  (discretization)
review gate:          not yet reviewed - Gate B needs L5, L6, L7
contextual radius:    4 (a typed relation against a production graph)

symptom / fact:       IDENTICAL BOOLEAN MATRICES CAN REPRESENT
                      DIFFERENT TYPED RELATIONS.

                      A production stencil carrying D2's occupancy
                      renders, through dependencies():

                        signature: vertices -> vertices
                                1 2 3
                        1       # # .
                        2       . # .
                        3       . . #

                      and the tower's D2 renders:

                        signature: X1 -> X2
                                p q r
                        u       # # .
                        v       . # .
                        w       . . #

                      Same grid, glyph for glyph. THREE DECLARED
                      CARRIERS stand behind those two pictures: X1,
                      X2, and the pattern's own 'vertices'. No two of
                      them are the same domain, and all three hold
                      three members.

                      So, with V the representation that forgets
                      carrier identity:

                        V(R1) = V(R2)   does NOT imply   R1 = R2

                      The tower's response is not to fix the
                      renderer but to SHOW THE SIGNATURE: every
                      Level-6 picture prints its signature above its
                      grid, because the grid alone is exactly the part
                      that cannot tell the two apart.

exact caller:         test/visualization-tower/level-6-discretization/
                      test.f90  check_the_typed_identity_differs,
                      check_the_visual_equality_theorem

evidence:             same_coordinate_pattern(D2, pattern) TRUE;
                      pattern's vertex_set() not same_as X1 and not
                      same_as X2; the two grids equal glyph for
                      glyph after labels and spacing are stripped.
                      PASS.

confidence:           high.

action:               OBSERVE. Classification: STENCIL-B.
```

---

## OBSERVATION VIZ-17

```text
tower:                Visualization
level:                6  (discretization)
review gate:          not yet reviewed
contextual radius:    4 (rectangular signature against one carrier)

symptom / fact:       THE PRODUCTION DEPENDENCY ANSWER CANNOT
                      INTRINSICALLY REPRESENT X -> Y WITH X /= Y.

                      D1 : X0 -> X1 runs 4 -> 3. The production
                      constructor takes ONE vertex count:

                        stored_graph(nv, tails=columns, heads=rows)
                        nv = size(constant)

                      so |V| would have to be 4 and 3 at once. Given
                      4, it holds every one of D1's five arrows - and
                      then reports a FOURTH ROW that D1's codomain
                      does not have:

                                1 2 3 4
                        1       # # . .
                        2       . # # .
                        3       . . . #
                        4       . . . .     <-- phantom

                      THE ARROWS SURVIVED. THE SIGNATURE DID NOT.

                      No union carrier was manufactured, no domain was
                      padded, and nothing was indexed out of range.
                      The phantom row is simply what one carrier
                      serving both axes looks like, read off safely.

                      THIS IS A SPECIALIZATION, NOT A DEFECT.
                      stencil_operator is a same-domain square-matrix
                      citizen, and its pattern contract is correct at
                      that radius. A defect would need an executable
                      promise that behaviour violates; see VIZ-22.

exact caller:         test/visualization-tower/level-6-discretization/
                      test.f90  check_the_rectangular_witness

evidence:             coordinate_shapes_fit(D1, pattern) FALSE;
                      pattern % num_vertices() = 4 while D1's codomain
                      holds 3; row 4 empty at every column;
                      same_coordinate_pattern refuses the comparison
                      rather than padding. PASS.

confidence:           high - the cardinality argument is exact, and
                      rests on the constructor's own single count.

action:               OBSERVE. Classification: RECT-B. No production
                      change. Rectangular support was NOT added
                      merely because this tower asked for it.
```

---

## OBSERVATION VIZ-18

```text
tower:                Visualization
level:                6  (discretization)
review gate:          not yet reviewed
contextual radius:    4 (what one concrete's verb means)

symptom / fact:       stencil_operator % dependencies() EXPOSES
                      ALGEBRAIC STATE SPARSITY.

                      Its graph is built at construction, one edge per
                      sparse coefficient, running

                        column -> row

                      over a single vertex carrier sized by the
                      constant vector. Read back through the graph's
                      own edge ends, the D2 stencil answers exactly
                      four arrows: 1->1, 2->1, 2->2, 3->3.

                      That is the sparsity of a matrix acting on a
                      state, and nothing else. It says which state
                      member feeds which.

exact caller:         test/visualization-tower/level-6-discretization/
                      test.f90  check_the_stencil_coordinate_pattern

evidence:             3 vertices, 4 edges, one per coefficient;
                      production_has agrees with the relation at every
                      cell; the rendered grid matches Level 4's. PASS.

confidence:           high.

action:               OBSERVE.
```

---

## OBSERVATION VIZ-19

```text
tower:                Visualization
level:                6  (discretization)
review gate:          not yet reviewed
contextual radius:    4 (what the other concrete's verb means)

symptom / fact:       step_operator % dependencies() EXPOSES THE
                      SCHEME'S TEMPORAL MOTIF.

                      BDF2 answers three vertices and two edges,
                      1 -> 2 -> 3: reach + 1 instants, each edge one
                      look further back. Rendered, it is a lower
                      subdiagonal - the shape of a HISTORY, not of a
                      sparsity.

                      It is manufactured fresh from `reach` at the
                      moment of asking, and it does not consult the
                      wrapped action at all.

                      THE PROBE THAT SETTLES IT. Two BDF2 steps were
                      built around two actions with genuinely
                      different state sparsity - D2's occupancy and a
                      diagonal - and dependencies() returned the
                      IDENTICAL motif for both. So the answer is
                      independent of what the step wraps.

                      Note the trap this specimen was built to spring:
                      the step's answer and its wrapped action's
                      answer BOTH have three vertices here. Equal
                      cardinality, entirely different semantics.

exact caller:         test/visualization-tower/level-6-discretization/
                      test.f90  check_the_step_pattern,
                      check_state_is_not_time,
                      check_the_motif_is_independent_of_what_it_wraps

evidence:             motif and wrapped pattern have equal vertex
                      counts and unequal adjacency; the stencil holds
                      1->1 and the motif does not; the motif holds
                      2->3 and the stencil does not; the two motifs
                      from two different wrapped sparsities are
                      extensionally identical. PASS.

confidence:           high.

action:               OBSERVE.
```

---

## OBSERVATION VIZ-20

```text
tower:                Visualization
level:                6  (discretization)
review gate:          not yet reviewed
contextual radius:    4 (one family verb, two concretes)

symptom / fact:       ONE FAMILY VERB DENOTES TWO DIFFERENT
                      STRUCTURAL AXES.

                        stencil % dependencies()  ->  state sparsity
                        step    % dependencies()  ->  temporal reach

                      Both are faithful to what their own concrete
                      represents. Neither is wrong. But the family
                      procedure does not have ONE invariant meaning
                      across its two citizens, and that was answered
                      from executable structure rather than from the
                      shared method name.

                      THE DEEPER READING, RECORDED AND NOT ACTED ON.
                      The pressure this exposes may not be

                        "we need one universal dependencies()"

                      but rather

                        "an operation may admit SEVERAL legitimate
                         structural interpretations"

                      D_state, D_time, D_space, D_block. The stencil
                      witness exposes one; the step witness exposes
                      another. That is evidence about the shape of the
                      problem, not a design.

                      DO NOT DESIGN THE ABSTRACTION. Two witnesses at
                      one radius is not enough to rule what the family
                      verb should mean, and this tower has inspected
                      no other operation family at all.

exact caller:         the whole of level-6-discretization/test.f90

evidence:             VIZ-18 and VIZ-19 together, plus the
                      independence probe. PASS.

confidence:           medium-high for the observation; LOW for any
                      design that might follow from it.

action:               OBSERVE. Classification: FAMILY-B.
```

---

## OBSERVATION VIZ-21

```text
tower:                Visualization
level:                6  (discretization)
review gate:          not yet reviewed
contextual radius:    4 (pressure toward the root)

symptom / fact:       THIS TOWER DOES NOT JUSTIFY STRUCTURAL
                      INTROSPECTION ON graph_operation.

                      Level 6 inspected discretization_operator and
                      its two concretes. It inspected NO other
                      operation family - not a linearization, not a
                      marcher, not a minimizer, not a partitioner.

                      The root question is:

                        Does EVERY graph_operation possess one
                        meaningful dependency structure?

                      Level 6 does not answer it, and VIZ-20 is if
                      anything evidence for caution: if two concretes
                      of ONE family already mean two different axes,
                      a single root-level verb across ALL families is
                      further from justified, not closer.

                      Nothing was added to graph_operation. No
                      dependencies(), no structure(), no visualize().
                      dependencies() was not moved upward and its
                      return type was not changed.

exact caller:         n/a - this observation is about what was NOT
                      done

evidence:             git diff b134a1f -- src/ is EMPTY.

confidence:           high.

action:               OBSERVE. Explicitly NOT REFACTOR.
```

---

## OBSERVATION VIZ-22

```text
tower:                Visualization
level:                6  (discretization)
review gate:          not yet reviewed
contextual radius:    4 (the contract's consumers)

symptom / fact:       dependencies() HAS TWO IMPLEMENTATIONS AND, AT
                      THE STARTING HEAD, NO CALLERS AT ALL.

                      The Level-6 census, run over every .f90 in the
                      repository at b134a1f:

                        discretization_operator      graph_calculus:219
                          deferred :: dependencies    interface :406-410
                          returns class(graph), allocatable

                        stencil_operator   class_graph_stencil.f90:47
                          stencil_dependencies                    :168
                          returns a copy of its own stored pattern
                          repo callers: 0

                        step_operator      class_graph_step.f90:46
                          step_dependencies                       :127
                          builds a fresh reach+1 chain
                          repo callers: 0

                        repo-wide callers of % dependencies():  ZERO

                      THREE FACTS THAT MUST NOT BE RUN TOGETHER:

                        the family contract exists
                        two implementations exist
                        no consumer exercises the contract

                      THE PROSE. graph_calculus says "the minimizers
                      one level up interrogate the pattern - the
                      diagonal, the colouring, the triangularity, the
                      Galerkin road - so it is exposed by law, never
                      by inspection." No minimizer does. That sentence
                      is recorded as

                        PROSE INTENTION - NO EXECUTABLE CALLER FOUND

                      and must not be repeated as an established
                      repository fact.

                      THIS IS NOT DEAD CODE. It is a LATENT contract:
                      two faithful implementations waiting for a
                      consumer. Level 6 is the first executable
                      consumer pressure it has ever had.

                      IT IS ALSO WHY NO RED WAS POSSIBLE. A RED needs
                      an executable promise that behaviour violates.
                      Design intent in a comment, with no caller
                      behind it, cannot justify modifying production.

exact caller:         the census itself, and
                      level-6-discretization/test.f90 as the first
                      caller

evidence:             grep over every .f90 in the repository returns
                      no `% dependencies(` outside the two procedure
                      bindings and the one deferred declaration.

confidence:           high - exhaustive over the tracked tree.

action:               OBSERVE. Records the exact count: ZERO.
```

---

## OBSERVATION VIZ-23

```text
tower:                Visualization
level:                6  (discretization)
review gate:          not yet reviewed
contextual radius:    4 (the whole Level-6 measurement)

symptom / fact:       A LIMITATION IS NOT A RED, AND LEVEL 6 PRODUCED
                      NONE.

                      Two findings could have been misread as defects
                      and are not:

                        "the ordinary dependency graph cannot
                         represent X0 -> X1 intrinsically"

                      is an architectural observation, because that
                      capability was never promised by executable
                      behaviour; and

                        "step dependencies means temporal reach while
                         stencil dependencies means state sparsity"

                      is not a failure merely because the two differ.
                      Each concrete faithfully implements the narrower
                      mathematics it actually represents.

                      A RED would require BOTH of:

                        1. an existing executable contract, test or
                           caller promising a behaviour
                        2. the implementation violating it

                      The census found no consumer at all (VIZ-22), so
                      condition 1 was never met.

                      NOTHING WAS APPLIED. stencil % apply and
                      step % apply were called zero times, and the
                      import gate refuses a `% apply(` in any tower
                      source so that is mechanical rather than
                      promised. No numerical composition of A3 A2 A1,
                      no transpose stencil, no A^T v.

exact caller:         the whole of level-6-discretization/

evidence:             git diff b134a1f -- src/ is EMPTY. L6 PASS.

confidence:           high.

action:               OBSERVE.
```

---

## OBSERVATION VIZ-24

```text
tower:                Visualization
level:                6  (discretization)
review gate:          not yet reviewed
contextual radius:    4 (multiple structural projections)

symptom / fact:       A SINGLE COMPUTATION MAY POSSESS SEVERAL
                      LEGITIMATE STRUCTURAL VIEWS.

                      The evidence, so far, is exactly two:

                        D_state   the stencil's sparsity
                        D_time    the step's motif

                      and the specimen showed them CO-EXISTING in one
                      composed object: a BDF2 step wrapping a stencil
                      has both, and dependencies() returns only the
                      temporal one.

                      Other axes are conceivable - D_space, D_block -
                      and this tower has seen none of them.

                      IF THIS HOLDS UP, the architectural question is
                      not "which one should dependencies() return"
                      but "how does an operation offer more than one".

                      DO NOT DESIGN THAT NOW. One tower, one radius,
                      two witnesses. Recording the shape of the
                      question is the whole of what is earned.

exact caller:         level-6-discretization/test.f90
                      check_state_is_not_time

evidence:             one composed object, two structures, and
                      dependencies() answering with one of them. PASS.

confidence:           medium - the observation is solid; its
                      generality is not established.

action:               OBSERVE. Do NOT promote a visualization or
                      introspection abstraction from this.
```

---

## OBSERVATION VIZ-25

```text
tower:                Visualization
level:                7  (minimization)
review gate:          NOT REACHED - Gate B needs L5, L6, L7
contextual radius:    5 (a solver over an operator over a host)

symptom / fact:       MINIMIZATION IS STRUCTURALLY INHABITED. LEVEL 7
                      IS NOT N/A.

                      Solver traversal already consumes graph
                      structure, and the census names the whole road:

                        minimizer % sweep_order()
                          -> walk(WALK_COLOURING)
                          -> colouring % apply(this % on)

                        minimizer % diagonal()
                          -> sweep_order()
                          -> one matvec per colour, with an
                             indicator vector

                        jacobi % solve()        -> diagonal()
                        gauss_seidel % solve()  -> diagonal(),
                                                   sweep_order()

                      So the question was never WHETHER structure is
                      consumed. It is WHOSE.

exact caller:         src/graph_minimization.f90 sweep_order :: diagonal;
                      src/class_graph_jacobi.f90 solve;
                      src/class_graph_gauss_seidel.f90 solve

evidence:             read from the local source at 7d4c501, and
                      exercised at
                      test/visualization-tower/level-7-minimization/
                      test.f90. PASS.

confidence:           high.

action:               OBSERVE. Level 7 is inhabited and must not be
                      marked N/A.
```

---

## OBSERVATION VIZ-26

```text
tower:                Visualization
level:                7  (minimization)
review gate:          NOT REACHED
contextual radius:    5 (where the colouring comes from)

symptom / fact:       sweep_order() DERIVES ITS COLOURING FROM THE
                      ATTACHED HOST GRAPH.

                        colouring % apply(this % on)

                      and `on` is the graph handed in at attach - the
                      legacy operation/minimizer context, which may be
                      unrelated to the operator's own coupling.

                      Three structures must be kept apart:

                        P_A     the operator's numerical coupling
                        H       the host/context graph
                        C(H)    what sweep_order() computes

                      and the present path is

                        H -> C(H) -> diagonal

                      never

                        P_A -> C(P_A) -> diagonal.

                      MEASURED, not inferred: one stencil attached
                      over two hosts gave colourings [1,2,1] and
                      [1,1,1]. Same operator, different host,
                      different sweep structure.

                      AND THE HOST'S COLOURING IS NOT A VALID
                      COLOURING OF THE OPERATOR. On the empty host,
                      unknowns 1 and 2 are coupled by A and share a
                      colour - which is exactly what a colouring is
                      supposed to prevent.

exact caller:         src/graph_minimization.f90 sweep_order;
                      level-7-minimization/test.f90
                      check_the_colouring_is_host_dependent

evidence:             both colourings recorded; properly_coloured(
                      col_match, P_A) TRUE and properly_coloured(
                      col_empty, P_A) FALSE. PASS.

confidence:           high.

action:               OBSERVE.
```

---

## OBSERVATION VIZ-27

```text
tower:                Visualization
level:                7  (minimization)
review gate:          NOT REACHED
contextual radius:    5 (operator against host)

symptom / fact:       A STENCIL'S NUMERICAL MAP IS INDEPENDENT OF ITS
                      HOST. This is the control the RED rests on.

                      stencil_apply computes entirely from

                        this % pattern

                      the stencil's OWN stored graph. The host is used
                      only for the output field's domain. So:

                        A_match x = A_empty x = [6, 14, 20]
                                                for x = [1,2,3]

                      measured through the minimizer's own matvec, not
                      by calling apply - the tower never calls apply,
                      and the import gate refuses one.

                      P_A is likewise host-free: dependencies() reads
                      the stored pattern, and answers the same graph
                      twice regardless of any attachment.

                      THE CONSEQUENCE, AND IT IS THE WHOLE POINT: any
                      difference caused by host topology alone belongs
                      to the MINIMIZER'S STRUCTURAL INTERPRETATION and
                      to nothing else. There is nowhere else for it to
                      come from.

exact caller:         level-7-minimization/test.f90
                      check_the_numerical_map_is_host_independent,
                      check_the_operator_owns_its_structure

evidence:             matvec equal under both hosts and equal to the
                      hand oracle; dependencies() idempotent. PASS.

confidence:           high.

action:               OBSERVE.
```

---

## OBSERVATION VIZ-28   *** RED ***

```text
tower:                Visualization
level:                7  (minimization)
review gate:          NOT REACHED
contextual radius:    5 (the reported diagonal)
status:               RED

status after repair:  RESOLVED - see the resolution block below.

symptom / fact:       THE REPORTED DIAGONAL CHANGES WHEN ONLY THE
                      IRRELEVANT HOST TOPOLOGY CHANGES.

                        A = [ 4 1 0 ]      true diagonal = [4, 5, 6]
                            [ 1 5 1 ]
                            [ 0 1 6 ]

                        on H_match = P_A     diagonal = [4, 5, 6]
                        on H_empty           diagonal = [5, 7, 7]

                      and [5,7,7] is exactly A times the all-ones
                      vector. The empty host colours everything 1, so
                      the coloured probe fires ONE indicator - all
                      ones - and reads a matvec where it expected a
                      diagonal.

                      THE FOUR FACTS, ESTABLISHED TOGETHER:

                        A_match x  =  A_empty x            YES
                        P_A same under both hosts          YES
                        H_match /= H_empty                 YES
                        d_match /= d_empty                 YES

                      That isolates the provenance error completely.
                      Nothing about the operator changed. Only the
                      graph the SOLVER was asked to interpret changed.

                      WHY THIS IS A GENUINE RED AND NOT A LIMITATION.
                      Level 6's findings were limitations because no
                      executable contract promised otherwise. This one
                      does: diagonal() means the diagonal of the
                      attached numerical action, and jacobi % solve
                      divides by it. An unchanged matrix has an
                      unchanged diagonal, and that is mathematics
                      rather than prose.

                      NOT WEAKENED. The oracle stays [4,5,6]. The
                      hostile host was not replaced, the expectation
                      was not lowered, no colour was hard-coded, and
                      the test asserts the truth and fails.

exact caller:         level-7-minimization/test.f90
                      check_the_diagonal
                      -> src/graph_minimization.f90 diagonal
                      -> sweep_order -> colouring % apply(this % on)

evidence:             FAIL : diagonal(A) on H_empty = [4,5,6] - THE
                      DIAGONAL OF AN UNCHANGED MATRIX DOES NOT DEPEND
                      ON AN IRRELEVANT CONTEXT GRAPH
                      FAIL : and the two agree with each other
                      captured verbatim before any production edit.

confidence:           high.

action:               RED, then REPAIRED.

resolution:           The minimizer gained a SECOND SEAT.

                        on         the execution context, handed to
                                   action % apply()
                        coupling   the dependent-variable stencil,
                                   the only thing sweep_order may
                                   colour

                      coupling is EXPLICIT AT ATTACH AND HAS NO
                      FALLBACK. `coupling := on` would have been
                      exactly the mistake the seat exists to prevent.
                      A structure-free minimizer never asks for one; a
                      structured one handed none now says so:

                        minimization: a sweep needs the
                        dependent-variable coupling - attach it with
                        coupling=

                      Callers that genuinely run on a graph which IS
                      the coupling of its unknowns - the mesh path -
                      now say so at their own call sites, where it is
                      a claim rather than an assumption.

                      AFTER: both hosts give colours [1,2,1] and
                      diagonal [4,5,6]. The hostile host was not
                      weakened; the solver simply stopped looking at
                      it. L7 PASS.
```

---

## OBSERVATION VIZ-29

```text
tower:                Visualization
level:                7  (minimization), resting on 6
review gate:          NOT REACHED
contextual radius:    5 (the obvious repair, refused)

symptom / fact:       dependencies() CANNOT BE SUBSTITUTED BLINDLY.

                      The tempting one-line repair is

                        colour action % dependencies()
                        instead of the host

                      and Level 6 already forbids it. VIZ-20 recorded:

                        stencil % dependencies() = state sparsity
                        step    % dependencies() = temporal motif

                      The family verb has NO INVARIANT MEANING across
                      its two concretes. Colouring its answer would
                      repair the stencil and silently break every
                      step - a BDF2 step would have its unknowns
                      coloured by a chain of instants that has nothing
                      to do with which unknowns are coupled.

                      Nor is a special case acceptable: no select type
                      on stencil_operator, because that would encode
                      "this family member happens to mean the right
                      thing" as architecture.

                      SO THE LEVEL-6 FINDING IS LOAD-BEARING HERE. It
                      is the reason a one-line fix is unavailable, and
                      the reason the RED is architectural rather than
                      local.

exact caller:         n/a - this observation is about what was NOT
                      done

evidence:             VIZ-20, plus git diff 7d4c501 -- src/ EMPTY.

confidence:           high.

action:               OBSERVE. Do NOT implement.
```

---

## OBSERVATION VIZ-30

```text
tower:                Visualization
level:                7  (minimization)
review gate:          NOT REACHED
contextual radius:    5 (what a consumer would have to ask for)

symptom / fact:       A CONSUMER MAY NEED TO REQUEST A PARTICULAR KIND
                      OF STRUCTURE, NOT MERELY "DEPENDENCIES".

                      The minimizer does not want "the operation's
                      dependency graph". It wants the ONE structure
                      its algorithm requires: which unknowns are
                      coupled, so that no two coupled unknowns share
                      a colour.

                      VIZ-24 already recorded that one computation may
                      admit several structural views - D_state,
                      D_time, and others unseen. Level 7 adds the
                      consumer's half of that: a caller has a
                      PARTICULAR view in mind, and today has no way to
                      name which.

                      Recorded as:

                        MINIMIZER STRUCTURAL PROVENANCE

                      and NOT as any of:

                        "dependencies() should move to
                         graph_operation"
                        "the ordinary graph needs replacing"
                        "add an explicit structural argument at
                         attach"
                        "add a typed operator-structure query"
                        "add named structural projections"

                      Those are candidate answers. Level 7 does not
                      choose among them, and one tower at one radius
                      is not enough to.

exact caller:         the whole of level-7-minimization/

evidence:             VIZ-26, VIZ-28 and VIZ-29 together.

confidence:           medium-high for the problem statement; NONE is
                      claimed for any remedy.

action:               OBSERVE. The next decision is architectural and
                      belongs to a review, not to this level.
```

---

## OBSERVATION VIZ-31

```text
tower:                Visualization
level:                6-7  (discretization, minimization)
review gate:          B
contextual radius:    5 (the family verb, correctly read)
supersedes:           the FAMILY-B reading in VIZ-20 and VIZ-24

symptom / fact:       dependencies() IS AXIS-RELATIVE, AND THAT IS ONE
                      MEANING RATHER THAN TWO.

                        dependencies() = the stencil on the axis this
                                         concrete type represents

                        stencil_operator -> DEPENDENT-variable stencil
                        step_operator    -> INDEPENDENT-variable
                                            stencil

                      Level 6 first read this as two unrelated
                      meanings, and that reading was wrong. The fault
                      was in step_dependencies, which answered

                        1 -> 2 -> 3

                      the SUCCESSION of instants. A BDF2 residual at
                      the newest instant reads three instants, so its
                      stencil is a FAN-IN

                        1 -> 3,  2 -> 3,  3 -> 3

                      and backward euler answers 1->2, 2->2. The
                      self-arrow is the implicit part - what makes the
                      newest instant an unknown rather than data.

                      A STENCIL IS NOT A CHRONOLOGY. Succession is a
                      true relation and the time integration tower
                      describes it correctly; it is simply not what
                      this contract owes, and that tower's own
                      chronology was left untouched.

                      THE INDEPENDENT AXIS NEED NOT BE TIME. A
                      continuation coordinate, a parameter, a spatial
                      sweep direction all take the same seat. The
                      concrete type carries the context so the verb
                      does not have to.

                      FAMILY-B is therefore WITHDRAWN and replaced by
                      FAMILY-C: the apparent semantic difference
                      disappears under the axis-relative reading.

exact caller:         src/class_graph_step.f90 step_dependencies;
                      level-6-discretization/test.f90
                      check_the_step_pattern

evidence:             BDF2 answers three vertices and three edges,
                      fanning in on the newest; the same stencil under
                      two differently-coupled wrapped actions; the
                      succession arrow 1->2 is absent and the
                      self-arrow present. PASS.

confidence:           high.

action:               PRODUCTION CORRECTED. This is the tower's second
                      production change and the only one to alter
                      behaviour rather than a seat.
```

---

## OBSERVATION VIZ-32

```text
tower:                Visualization
level:                7  (minimization)
review gate:          B
contextual radius:    5 (who supplies structure to whom)

symptom / fact:       THE CONSUMER IS HANDED ITS STRUCTURE; IT DOES
                      NOT REACH FOR IT.

                      VIZ-30 asked where a consumer should ask for a
                      particular structural projection. The answer
                      this tower reached needs no new production
                      vocabulary at all:

                        stencil_operator
                              |  dependencies()
                              v
                        dependent-variable stencil
                              |
                              v
                        minimizer % attach(..., coupling = ...)

                      The CALLER knows which object owns the dependent
                      axis. Nothing in the minimizer inspects an
                      action's type, and nothing asks a graph_operation
                      for a structure it may not have.

                      WHAT WAS NOT BUILT, and was not needed:

                        dependencies() on graph_operation
                        a structural-axis enum
                        graph_observation
                        graph_visualization
                        a structural_projection type
                        a generic visualization abstraction

                      A step_operator's dependencies() is an
                      INDEPENDENT-axis stencil and is therefore never
                      passed as a minimizer's coupling. That is not a
                      special case; it is a caller reading a type it
                      already knows.

exact caller:         src/graph_minimization.f90 attach :: sweep_order;
                      level-7-minimization/test.f90 attach_to

evidence:             L7 PASS with coupling = P_A under two unrelated
                      execution contexts; the no-fallback refusal
                      fires when a structured solver is handed none.

confidence:           high at this radius.

action:               OBSERVE. The seam is closed for this consumer
                      without a root abstraction.
```

---

## OBSERVATION VIZ-33

```text
tower:                Visualization
level:                8-9  (constitution, statement)
review gate:          C
contextual radius:    5 (the whole sealed client)

symptom / fact:       THREE STRUCTURES COEXIST IN ONE COMPOSITION
                      WITHOUT BEING CONFLATED, AND NO NEW ABSTRACTION
                      WAS NEEDED TO KEEP THEM APART.

                        independent-variable stencil   step's own
                        dependent-variable stencil     stencil's own
                        execution context              neither

                      all three carrying three members in the Level-8
                      specimen - the coincidence that makes the
                      distinction worth proving rather than asserting.

                      The solver, handed the dependent stencil
                      explicitly, produced a colouring that is valid
                      for THAT structure, invalid for the independent
                      one, and not the flat colouring the empty
                      context would have given. Diagonal [4,5,6].

                      NO KRONECKER PRODUCT WAS FORMED and no
                      product-space type exists. The axes were shown
                      to be distinct and composable; building the
                      product would have been choosing an abstraction
                      on one specimen's evidence.

                      WHAT THE TOWER DID NOT BUILD, across all ten
                      levels:

                        graph_visualization
                        graph_observation
                        graph_representation
                        graph_interpretation
                        graph_renderer
                        visualize() / print() / structure() on any root
                        dependencies() on graph_operation
                        a structural-axis enum
                        a structural_projection type
                        a product-space type

                      The renderers are test fixtures and stay test
                      fixtures. Their whole evidential value is that
                      they may use nothing but the nucleus.

exact caller:         level-8-constitution/test.f90,
                      level-9-statement/test.f90

evidence:             L8 PASS, L9 PASS, tower sealed.

confidence:           high at this radius.

action:               OBSERVE, AND BUILD NOTHING. A visualization
                      abstraction was not earned.
```

---

## Seam accounting at the seal

| Seam | Before this tower | After Gate A | Why |
|---|---|---|---|
| **A1** graph host as conduit | CLOSED | CLOSED | nothing here reopens it; no operation was hosted anywhere |
| **A2** operations take domain from graph | 3 towers | **3 towers — unchanged** | Gate A attached no operation to anything; Level 5 attached a FIELD to a CARRIER; Level 6 interrogated operations without applying one. Levels 7–9 applied one through a minimizer — and the fault found there was the solver taking its COLOURING from the host, which is a different question from an operation taking its DOMAIN from one. The seam is untouched |
| **A3** relational_graph ownership | KEEP | KEEP | one more successful typed-ownership pattern (7 carriers, 6 relations, full signature closure); no production change follows |
| **B** bidirectional/rectangular linearization | 2 towers | **2 towers — ZERO new votes** | this tower is full of rectangular structure and none of it is a linearization; structural transpose is not numerical adjoint, and Level 5 built no `w^T`, applied no `A^T v`, and composed no coefficients |

**This tower is ONE client, not eight.** Level 7 produced this
tower's first RED, and it is resolved. It is still **one** tower's
evidence, and it is not a vote for seam A2: the minimizer's fault was
taking its COLOURING from the host, which is a different question from
an operation taking its DOMAIN from one.

---

## Frontier

Level 5 asked *"does a numerical zero erase structural dependency?"*
and answered **no**. Level 6 asked whether a production discretization
operator exposes the same structural skeleton, and answered
**STENCIL-B / RECT-B / FAMILY-B**: coordinate-equivalent but not
typed-equivalent; unable to hold a rectangular signature; and one verb
denoting two different axes.

Level 7 asked whether minimization consumes discretization structure.
It does — through colouring — and it colours the **host**, not the
operator. The diagonal of an unchanged matrix moved when only an
irrelevant context graph changed.

The decision it demanded turned out to need no new vocabulary: the
CALLER knows which object owns the dependent axis, and hands that
structure to the minimizer at attach. See VIZ-32.

```
TOWER SEALED.

    relation        ->  structure
    field           ->  values
    dependencies()  ->  the stencil on the concrete type's axis
    renderer        ->  representation
    minimizer       ->  explicit dependent-variable coupling

    context graph != independent stencil != dependent stencil
```

Production changed in three places, all narrow and all earned:
class_graph_step (a stencil is not a chronology), graph_calculus
(the contract's prose), graph_minimization (a seat for the coupling),
plus the caller adjustments the new argument required. No root was
generalized and no visualization abstraction was promoted.
