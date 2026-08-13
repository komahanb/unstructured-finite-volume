# Nucleus observations — Adjoint Sensitivity Tower

The live evidence ledger of this orbital client, per the discipline of
[FRACTAL-GRAPH-ARCHITECTURE.md](../../FRACTAL-GRAPH-ARCHITECTURE.md):

\[
\boxed{\text{observation}\neq\text{immediate refactor}}
\]

Entries record load-bearing architectural evidence, not a diary.
*Local necessity* means necessity at the object's own radius; *global
necessity* means necessity at larger contextual radius. Recurrence
across sibling towers is its own field and is never used as proof of
global necessity.

---

## OBSERVATION AD-0

```text
tower:                Adjoint Sensitivity
gate / level:         A / 0
contextual radius:    0

symptom / fact:       |Q| = |Y| = 2 and Q is NOT Y. The state domain
                      and the residual domain of a square implicit
                      problem have equal cardinality and no
                      relationship whatever; only same_as separates
                      them. Four roles - P, Q, Y, Z - are ordinary
                      subobjects of two parents, and no
                      parameter_field/state_field/residual_field/
                      adjoint_field type was needed or invented.

exact caller:         level-0-carrier/test.f90
                      (check_role_identities)

mathematical concept: role as domain; identity vs cardinality

local necessity:      yes - every level above must keep them apart,
                      and a solver that collapsed them would answer
                      a plausible wrong number (A is non-symmetric)
global necessity:     unknown at larger radius by this gate

cross-tower recurrence: learning distinguished Theta from Y with
                      |Theta| = |Y| = 1; this tower raises the stakes
                      because BOTH domains are square-sized and the
                      operator between them is not symmetric

graph role:           not applicable at this rung

comparison:           the first tower where confusing two domains
                      would produce a wrong ANSWER rather than a
                      compile error or an obvious mismatch

suspected nucleus implication: none - carriers and subobjects
                      sufficed; the temptation lives in client code,
                      and the nucleus already refuses it by identity

confidence:           high
action:               observe
```

---

## OBSERVATION AD-1

```text
tower:                Adjoint Sensitivity
gate / level:         A / 1
contextual radius:    0

symptom / fact:       this specimen's dependency relation is DENSE:
                      R_dep = T x V, all nine facts. There is no
                      in-range absence to report, so membership
                      queries carry no structural information at all
                      here.

                      Worse (and instructive): the two carriers
                      enumerate from one, so their raw ids OVERLAP -
                      TGT_R1 and VAR_P are both the integer 1 - and a
                      tuple written in the wrong orientation,
                      [VAR_U, TGT_R1] = [2,1], is a perfectly good
                      member reading as (r2, p). Membership cannot
                      even detect orientation.

                      Orientation is carried by the SIGNATURE,
                      checked by identity.

exact caller:         level-1-relation/test.f90 (check_absences) -
                      an assertion that ORIGINALLY claimed the
                      reversed pair was absent, and failed; the code
                      was right and the assertion was wrong

mathematical concept: extension vs signature as carriers of meaning

local necessity:      the signature check is locally necessary; the
                      membership check is locally worthless here
global necessity:     unknown

cross-tower recurrence: the opposite of Derivative Action, whose
                      whole message came from a SPARSITY pattern
                      ((x,z) absent yet reachable). Two derivative-
                      family clients, two different sources of
                      structural content: sparsity there, role
                      partition and orientation here

graph role:           not applicable

comparison:           a dense relation is not a degenerate one - it
                      says every derivative block is structurally
                      full, and the blocks are then told apart by
                      their domains alone

suspected nucleus implication: none. But a caution worth carrying:
                      any future tooling that reports structure by
                      probing membership across carriers with
                      overlapping id spaces can be silently fooled;
                      signatures cannot.

confidence:           high
action:               observe
```

---

## OBSERVATION AD-2

```text
tower:                Adjoint Sensitivity
gate / level:         A / 2
contextual radius:    0

symptom / fact:       all four derivative supports came from ONE
                      dependency relation and the roles' own
                      relational faces:

                        J_Q = I_Q^T o R_dep o I_Y   <= Y x Q
                        J_P = I_P^T o R_dep o I_Y   <= Y x P
                        F_Q = I_Q^T o R_dep o I_Z   <= Z x Q
                        F_P = I_P^T o R_dep o I_Z   <= Z x P

                      (relational composition, read right to left;
                      in code the same road is written left to right
                      as compose_binary(compose_binary(I_Y,R_dep),
                      I_Q^T), since compose_binary(R_AB,R_BC) =
                      R_BC o R_AB)

                      inclusion_of() turned each subobject's
                      membership into algebra, and the right-hand
                      selectors were used TRANSPOSED - borrowed
                      views, never rebuilt. Four blocks, one stored
                      truth, no hand-maintained structure.

exact caller:         level-2-relation-algebra/test.f90

mathematical concept: restriction as composition with an inclusion

local necessity:      yes - four independently authored blocks could
                      drift from the source and from each other
global necessity:     unknown

cross-tower recurrence: first client to use inclusion_of at all.
                      The calculator/learning/derivative towers
                      restricted with restrict_slot on ports; this
                      tower restricts by ROLE, and the subobject's
                      relational face is exactly the right tool -
                      it was already in production, unexercised by
                      any prior tower

graph role:           not applicable

comparison:           the blocks are distinguishable ONLY by domain
                      identity: J_Q and F_Q share their second slot,
                      J_Q and J_P share their first, and no size
                      comparison separates any pair

suspected nucleus implication: none - the earned primitives
                      sufficed. Worth recording that inclusion_of
                      has now found its first real caller.

confidence:           high
action:               observe
```

---

## OBSERVATION AD-3

```text
tower:                Adjoint Sensitivity
gate / level:         A / 3
contextual radius:    0-1

symptom / fact:       a relational_graph seated SUBOBJECTS as
                      ordinary citizens beside their parents: six
                      carriers (V, T and the four roles) and five
                      relations, with signature closure holding over
                      MIXED ownership - R_dep's slots resolve to the
                      parents, the four blocks' slots to the
                      subobjects, and every one is a carrier the
                      graph holds.

exact caller:         level-3-graph/test.f90
                      (check_roles_are_citizens,
                      check_signature_closure)

mathematical concept: ownership closure over a domain lattice

local necessity:      yes - the blocks cannot be owned unless their
                      role domains are
global necessity:     unknown

cross-tower recurrence: prior towers owned only parent carriers;
                      this is the first graph in the repository whose
                      carrier list contains both a parent and its
                      subobjects

graph role:           owner / structural closure. Nothing more was
                      asked of it at this rung: no adjoint seat, no
                      transpose seat, no derivative metadata

comparison:           the container needed no change to accept a
                      subobject - a subset_set is a member_set, and
                      held_set took it without ceremony

suspected nucleus implication: none. Evidence that the domain
                      lattice and the ownership contract compose
                      cleanly, which a coupled/partitioned client at
                      larger radius will lean on harder.

confidence:           high
action:               observe
```

---

## OBSERVATION AD-4

```text
tower:                Adjoint Sensitivity
gate / level:         A / 4
contextual radius:    0-1

symptom / fact:       THE structural break with every sibling tower.
                      Calculator, Learning and Derivative Action all
                      found an execution order at level 4. This
                      tower's state is implicitly coupled:

                        C_Q = J_Q o J_Q^T = {(u,u),(u,v),(v,u),(v,v)}

                      (the path Q -> Y -> Q; as a Boolean matrix
                      pattern the same object reads J_Q^T J_Q)

                      derived from the same J_Q. The directed view
                      builds and is perfectly valid; reachability
                      answers both ways; there are NO sources and NO
                      sinks; and topological_order REFUSES with
                      'a topological order needs an acyclic graph'.

                      A valid directed graph is not a DAG. No acyclic
                      fiction was invented to make the walk succeed.

exact caller:         level-4-graph-calculus/test.f90 and
                      refusal.f90 (case cyclic-order)

mathematical concept: implicit coupling vs executable dependency

local necessity:      yes - an implicit system HAS no execution
                      order, and the level's whole content is saying
                      so honestly
global necessity:     unknown; larger implicit systems (blocked,
                      partitioned, multiphysics) will meet the same
                      wall

cross-tower recurrence: first tower whose level-4 truth is NEGATIVE.
                      Derivative Action's DA-4 kept 'dependency
                      relation / directed interpretation / derivative
                      traversal' distinct; this tower adds a fourth
                      distinction: interpretation exists, and
                      EXECUTION does not follow from it

graph role:           interpretation of graph-owned structure; the
                      algorithms are here to refuse, not to schedule

comparison:           the calculator walked [plus, times]; this
                      tower cannot walk at all, and that is the
                      correct answer. What solves an implicit system
                      is minimization (level 7), and the ladder's
                      shape says so

suspected nucleus implication: none needed - the nucleus ALREADY
                      refuses correctly. Recorded as strong evidence
                      that graph calculus and equation solving are
                      genuinely different contracts, and that a
                      level may be inhabited by a refusal

confidence:           high
action:               observe
```

---

## OBSERVATION AD-5

```text
tower:                Adjoint Sensitivity
gate / level:         A / 5
contextual radius:    0

symptom / fact:       the first numbers are p = [2] on P and a
                      deliberately WRONG q0 = [0,0] on Q, both
                      ordinary production fields distinguished only
                      by their domains. Nothing else was built - in
                      particular NO lambda = 0.

exact caller:         level-5-field-calculus/test.f90

mathematical concept: initial data vs fabricated answer

local necessity:      yes for the two fields; the absent ones are
                      locally unnecessary and would be actively
                      harmful
global necessity:     unknown

cross-tower recurrence: third confirmation of roles-are-domains
                      (learning L5, derivative L5, here). New here:
                      the refusal to pre-create an ADJOINT field.
                      A zero adjoint is not an empty adjoint - it is
                      a claim about an uncomputed solution, and a
                      later bug that never overwrote it would pass
                      unnoticed

graph role:           none - a field needs a domain, not a graph
                      (confirmed a fourth time)

comparison:           learning refused to fabricate values on its
                      computed domain U; this tower refuses to
                      fabricate a solution field for an equation not
                      yet posed

suspected nucleus implication: none

confidence:           high
action:               observe
```

---

## OBSERVATION AD-6

```text
tower:                Adjoint Sensitivity
gate / level:         A / 6
contextual radius:    0-1

symptom / fact:       the adjoint support is the transpose VIEW of
                      the primal support, and the transpose swaps
                      DOMAIN IDENTITIES, not array indices:

                        J_Q   <= Y x Q   supports  Rq   : Q -> Y
                        J_Q^T <= Q x Y   supports  Rq^T : Y -> Q

                      Both slots hold two members; the domains are
                      still different; J_Q is materialized and its
                      transpose is not. No second reverse relation
                      exists anywhere.

exact caller:         level-6-discretization/test.f90
                      (check_orientation_is_identity,
                      check_one_stored_truth)

mathematical concept: transpose as reorientation between domains

local necessity:      yes - and uniquely consequential here: with
                      |Q| = |Y| and A non-symmetric, a reversed
                      adjoint would pass every dimension check and
                      return a wrong, plausible number
global necessity:     unknown

cross-tower recurrence: FOURTH tower using transpose-as-view
                      (calculator L6, learning L6, derivative L6,
                      here). Strongest recurring behaviour in any
                      ledger so far. New here: the two ends are
                      different domains of the SAME size, so the
                      view's value is semantic, not merely storage
                      economy

graph role:           view over a materialized citizen

comparison:           a support relation is written (row, column),
                      so its first slot is the operator's CODOMAIN.
                      Reading a support as the operator's direction
                      is the easiest way to get an adjoint backwards

suspected nucleus implication: none yet. Gate B will test whether
                      the NUMERICAL reverse action can be expressed
                      as cleanly as this structural one - the seam
                      Derivative Action flagged (DA-8B) and this
                      tower will press at solver radius.

confidence:           high
action:               observe; press at Gate B
```

---

## OBSERVATION AD-7

```text
tower:                Adjoint Sensitivity
gate / level:         B / 7
contextual radius:    1 (an equation inside a solver)

symptom / fact:       ONE minimizer family governed both orientations.
                      The same gmres TYPE was attached twice:

                        unknown Q, residual Y  ->  q      = [2,4]
                        unknown Y, residual Q  ->  lambda = [-0.4,0.6]

                      No adjoint_solver, transpose_gmres or
                      reverse_gmres was needed or written. The solver
                      never learns the word adjoint: it knows an
                      unknown domain, a residual domain, and an
                      operation that answers on the latter.

exact caller:         level-7-minimization/test.f90
                      (check_primal_solve, check_adjoint_solve,
                      check_one_solver_family)

mathematical concept: the transposed problem is a different EQUATION,
                      not a different kind of solving

local necessity:      not applicable - this is a NEGATIVE result about
                      what was NOT needed
global necessity:     unknown; a symmetric/SPD family or a Krylov
                      method requiring the transpose action internally
                      might yet distinguish them

cross-tower recurrence: learning L7 and derivative-action never posed
                      a transposed equation at all; this is the first
                      evidence in any tower that unknown/residual
                      orientation can be exchanged with no solver
                      change

graph role:           see AD-9

comparison:           the minimizer's explicit unknown_domain and
                      residual_domain - earned during the learning
                      tower's Phase 5B - is exactly what makes this
                      possible. A solver that inferred its domains
                      from the host could not have done it

suspected nucleus implication: NONE. The hypothesis "an adjoint solver
                      is a new solver family" is contradicted for this
                      radius; record and move on.

confidence:           high
action:               observe
```

---

## OBSERVATION AD-8

```text
tower:                Adjoint Sensitivity
gate / level:         B / 7
contextual radius:    1

symptom / fact:       Gate A's theorem became numerical and PRODUCTION
                      enforced it. A right-hand side on the wrong
                      domain has exactly the right SIZE when
                      |Q| = |Y| = 2, and graph_minimization refuses it
                      by identity:

                        'a right-hand side lives on the residual domain'

                      Two such refusals (primal rhs on Q, adjoint rhs
                      on Y) are pinned by message, plus two more from
                      the test-local equations refusing states on the
                      wrong unknown domain.

exact caller:         level-7-minimization/refusal.f90 and
                      check_refusals.sh; the production check lives in
                      src/graph_minimization.f90 (solver_apply)

mathematical concept: dimension agreement vs domain identity

local necessity:      yes, and sharply: nothing but same_as can reject
                      these fields
global necessity:     unknown, but the risk grows with radius - every
                      square system offers this trap

cross-tower recurrence: the minimizer's domain check existed before
                      this tower and was never exercised adversarially;
                      this is its first hostile test

graph role:           not applicable

comparison:           had the check compared sizes, all four wrong
                      fields would have been accepted and the tower
                      would have returned plausible wrong numbers

suspected nucleus implication: none needed - production already holds.
                      Recorded as evidence that the identity-based
                      contract is load-bearing, not decorative.

confidence:           high
action:               observe
```

---

## OBSERVATION AD-9

```text
tower:                Adjoint Sensitivity
gate / level:         B / 7-8
contextual radius:    1

symptom / fact:       the legacy graph host is COMPATIBILITY at this
                      radius. The graph_operation face requires a
                      class(graph); a five-vertex stored_graph was
                      supplied whose vertex set is neither Q nor Y and
                      not even their size. It is never queried for
                      topology, never influences q or lambda, and
                      every operation ignores it via
                      `associate (u1 => input_graph); end associate`.
                      GMRES itself needs only matvec, inner_product
                      and norm.

exact caller:         level-7-minimization/test.f90
                      (check_host_is_nobodys_domain); every apply() in
                      both fixtures

mathematical concept: graph as operand vs context vs compatibility

local necessity:      NO as an operand - the mathematics never reads it
global necessity:     unknown - deliberately. Partitioned, coupled or
                      distributed adjoints may need context; nothing
                      here recommends removing anything

graph role classification (the evidence, not a verdict):
                        operand?               no
                        context?               not used as one here
                        ownership environment? not at this radius
                        compatibility host?    YES - this is what it is

cross-tower recurrence: THIRD tower to report it - learning L7/L9,
                      derivative-action (DA-8F), and now here. The
                      strongest recurring seam in the derivative
                      family. New this time: the host is provably
                      distinct from BOTH the unknown and residual
                      domains, which sharpens the earlier reports

comparison:           derivative action found the graph locally
                      unnecessary as an operand but load-bearing as an
                      ownership environment at its statement (DA-9A);
                      Gate C of this tower is where the same question
                      must be asked again

suspected nucleus implication: the operation-host contract remains the
                      repository's oldest live seam. Three towers now
                      agree it is not an operand at radius 0-1. Still
                      NOT a refactor: the Case-III caution applies, and
                      the field_operation proposal drafted earlier
                      should be re-examined only with a higher-radius
                      client in hand.

confidence:           high (the fact); none offered (the design)
action:               test at larger radius
```

---

## OBSERVATION AD-10

```text
tower:                Adjoint Sensitivity
gate / level:         B / 8
contextual radius:    1

symptom / fact:       ONE coefficient table generated everything. Each
                      entry of A, b, c, d is written once, keyed by the
                      MEMBERS it relates, and the following all read
                      it: the primal residual, the response, the
                      forward state action, the reverse state action,
                      the parameter action, and both readings of the
                      response block.

                      There is no A^T in the file. The reverse action
                      walks the same J_Q with the same coeff_state,
                      accumulating with += into the state slots; and
                      f_q^T is obtained as the response block's reverse
                      action with a unit seed, never written down.

exact caller:         level-8-constitution/adjoint_constitution_fixture.f90
                      (coeff_state, rq_forward, rq_reverse, fq_reverse)

mathematical concept: one linearization, read in two directions

local necessity:      yes - a second transposed table could drift from
                      the first, and with a non-symmetric A the drift
                      would be silent
global necessity:     unknown

cross-tower recurrence: SECOND derivative-family client to establish
                      it (derivative action's DA-8B was the first).
                      New here: the reverse reading now feeds a SOLVE
                      rather than a single evaluation, which is the
                      stronger demand

graph role:           none - the actions take relations, domains and
                      arrays

comparison:           derivative action applied one local linearization
                      per operation along an execution order; this
                      tower applies one coefficient table across a
                      structural support with no order at all. Same
                      law, different traversal

suspected nucleus implication: see AD-12 - two clients now share a
                      PROTOCOL shape, which is worth watching but not
                      yet worth promoting

confidence:           high
action:               observe
```

---

## OBSERVATION AD-11

```text
tower:                Adjoint Sensitivity
gate / level:         B / 8
contextual radius:    1

symptom / fact:       every numerical action is driven by the Gate-A
                      STRUCTURAL SUPPORTS, and every vector is indexed
                      through its domain's own local_index. Nothing
                      loops over "two rows and two columns".

                      The whole Level-8 battery therefore runs twice:
                      canonically, and with BOTH two-member roles
                      independently reversed (Q' = [v,u],
                      Y' = [r2,r1]). Every answer is checked by member
                      and none moves.

                      This was verified adversarially: replacing the
                      member-keyed lookup with a positional one
                      (reading A as a literal 2x2 array indexed by
                      position in Y and Q) passes the canonical run
                      and breaks SEVEN truths in the permuted one.

exact caller:         level-8-constitution/test.f90 (run_battery,
                      called twice with asserted storage positions)

mathematical concept: structure decides participation; enumeration
                      decides only storage

local necessity:      yes
global necessity:     unknown, but any client whose roles are built by
                      partitioning, sorting or renumbering will meet
                      exactly this

cross-tower recurrence: the strongest form yet of the roles-are-domains
                      discipline: not merely "role is a domain" but
                      "no answer may depend on how a role is
                      enumerated"

graph role:           the supports are relations; the graph owns them
                      at Gate A

comparison:           earlier towers pinned declaration order (K={y,x},
                      X={y,x}) to defeat raw-index assumptions in ONE
                      domain. This tower permutes TWO domains
                      independently, which is what catches a
                      cross-domain positional assumption

suspected nucleus implication: none - the nucleus already keyed
                      everything by member. The lesson is for CLIENTS.

confidence:           high
action:               observe
```

---

## OBSERVATION AD-12

```text
tower:                Adjoint Sensitivity
gate / level:         B / 8
contextual radius:    1

symptom / fact:       TWO questions, kept apart deliberately.

                      (a) The existing production
                      difference_linearization builds its perturbed
                      state and its Jv result on
                      input_graph % vertex_set() - it is same-domain
                      and graph-vertex specialized. It cannot express
                      Rq : Q -> Y with Q not same_as Y. This tower did
                      NOT contort Q and Y into vertices to reuse it,
                      and did NOT pre-emptively generalize it.

                      (b) Across the two derivative-family towers, the
                      duplicated material is now visible. What is
                      duplicated is APPLICATION LAW - coefficients,
                      residual and response definitions - which is
                      healthy. What RECURS is a protocol shape:

                        an input domain
                        an output domain
                        a forward action
                        a reverse action reading the same law
                        the transpose orientation between the two

exact caller:         src/class_graph_linearization.f90 (lines ~156-170)
                      for (a); the two towers' fixtures for (b)

mathematical concept: a general linearization contract between
                      DIFFERENT domains

local necessity:      the test-local adapters sufficed here
global necessity:     unknown

cross-tower recurrence: TWO independent clients have now written the
                      same protocol shape by hand. Per the fractal
                      document's scale that is a "recurring seam", not
                      yet "strong architectural evidence" - three
                      independent towers, or one high-radius failure,
                      would raise it

graph role:           difference_linearization's dependence on
                      vertex_set is precisely a graph-as-domain
                      assumption, and it is what makes it unusable
                      here

comparison:           the seam derivative action flagged (DA-8B) as
                      "can the numerical reverse action be expressed as
                      cleanly as the structural one" is answered YES -
                      but by a test-local adapter both times

suspected nucleus implication: a future general linearization/action
                      contract carrying its own input and output
                      domains is the candidate. It is NOT earned yet.
                      Explicitly NOT added during Gate B:
                      apply_transpose, adjoint_operation,
                      bidirectional_operator - naming a thing is not
                      the same as needing it.

confidence:           medium
action:               observe; revisit if a third client or a
                      high-radius tower repeats the shape
```

---

## OBSERVATION AD-13

```text
tower:                Adjoint Sensitivity
gate / level:         C / 9
contextual radius:    1 (the complete statement)

symptom / fact:       "GRAPH" CARRIES TWO DIFFERENT SOFTWARE ROLES,
                      and this statement holds both at once.

                      MODEL (relational_graph): owns V,T,P,Q,Y,Z and
                      R_dep, J_Q, J_P, F_Q, F_P. The statement
                      locates its citizens by identity, DESTROYS
                      every construction selector, and every number
                      afterwards arrives through model-owned
                      relations. Load-bearing.

                      HOST (stored_graph): seven vertices in a chain,
                      provably neither Q nor Y nor their size, not
                      owned by the model, never queried for topology
                      or coefficients. It exists because
                      graph_operation demands a class(graph).
                      Compatibility scenery.

                      Therefore both of these are true at once:

                        the MODEL graph is necessary as an
                            ownership environment
                        the SOLVER HOST is unnecessary as a
                            numerical operand

exact caller:         level-9-statement/test.f90
                      (check_model_ownership,
                      check_host_is_not_the_model)

mathematical concept: ownership environment vs operand

local necessity:      model: YES at statement radius - the selectors
                      are gone and the mathematics continues
                      host: NO - nothing reads it
global necessity:     unknown for both

cross-tower recurrence: FOURTH tower to meet the operand question and
                      the SECOND to answer it in both directions at a
                      statement:
                        learning L9      ownership load-bearing
                        derivative DA-9A operand no / owner yes
                        adjoint AD-9     host = compatibility (B)
                        adjoint AD-13    model = owner (C)
                      Two independent derivative-family clients now
                      report the same split at the same radius. On
                      the fractal document's scale this is
                      approaching STRONG architectural evidence

graph role:           both roles, simultaneously, in one program -
                      which is precisely why "is the graph
                      necessary?" has no answer until the question
                      names WHICH graph

comparison:           earlier towers could conflate the two because
                      their model structure and their solver host
                      never coexisted in one statement. Here they do,
                      and they are different types holding different
                      things

suspected nucleus implication: the operation-host contract is the
                      repository's oldest live seam and is now a
                      REVERSE-MODE REFACTOR CANDIDATE - to be
                      considered AFTER this tower is sealed, with the
                      field_operation sketch re-examined in the light
                      of four towers rather than one. NOT during
                      Gate C, and not as a mutation justified by a
                      2x2 specimen.

confidence:           high (the fact); medium (the readiness)
action:               reverse-mode refactor candidate, after the tower
```

---

## OBSERVATION AD-14

```text
tower:                Adjoint Sensitivity
gate / level:         C / 9 (whole-tower)
contextual radius:    1

symptom / fact:       the complete adjoint road required NO
                      adjoint-specific production ontology:

                        no adjoint field type      no adjoint solver
                        no adjoint graph           no adjoint Jacobian
                        no adjoint statement class no sensitivity type

                      lambda is an ordinary field distinguished by
                      living on Y; the adjoint equation is an
                      ordinary operation with its unknown and
                      residual domains exchanged; the adjoint solve
                      is the ordinary minimizer; the adjoint support
                      is a transpose view.

exact caller:         the whole of level-9-statement/test.f90; the
                      absence of any new src/ file across three gates

mathematical concept: "adjoint" as a USE of general structure

local necessity:      not applicable - a negative result
global necessity:     not applicable

cross-tower recurrence: the fourth consecutive tower to close with
                      zero production changes (calculator, learning,
                      derivative action, adjoint)

graph role:           see AD-13

comparison:           the adjective names a mathematical use of
                      ordinary domains, fields, actions, minimization
                      and transpose interpretation - not a family of
                      software nouns. Naming a thing is not the same
                      as needing it

suspected nucleus implication: none. Recorded as the tower's headline
                      architectural result.

confidence:           high
action:               observe
```

---

## OBSERVATION AD-15

```text
tower:                Adjoint Sensitivity
gate / level:         C / 9
contextual radius:    1

symptom / fact:       the one-law-two-directions protocol has now
                      been used at three increasing radii:

                        derivative action   evaluation only
                        adjoint Gate B      inside a SOLVE
                        adjoint Gate C      composed into a complete
                                            statement, with f_q^T
                                            itself generated by the
                                            reverse reading

                      In every case the duplicated MATERIAL was only
                      application law. What recurred was the protocol
                      shape: an input domain, an output domain, a
                      forward action, a reverse action reading the
                      same law, and the transpose orientation between
                      them.

exact caller:         adjoint_constitution_fixture.f90 used by both
                      level-8 and level-9; derivative action's
                      fixture for the earlier data point

mathematical concept: a general linearization contract between
                      DIFFERENT domains

local necessity:      test-local adapters sufficed every time
global necessity:     unknown

cross-tower recurrence: TWO independent towers, THREE radii. By the
                      ledger's own threshold ("two unrelated towers:
                      recurring seam; three or more independent
                      towers: strong architectural evidence") this is
                      a recurring seam that has now been stressed at
                      a larger radius rather than a third independent
                      client. It does NOT yet clear the bar

graph role:           the existing production
                      difference_linearization cannot serve because
                      it builds on input_graph % vertex_set()
                      (AD-12) - a graph-as-domain assumption

comparison:           the seam has strengthened without a third
                      client appearing; that is exactly the situation
                      the discipline says to record rather than act on

suspected nucleus implication: a general linearization/action
                      contract carrying its own input and output
                      domains remains THE candidate. Marked a
                      reverse-mode refactor candidate for AFTER the
                      tower. Explicitly not added during Gate C:
                      apply_transpose, adjoint_operation,
                      bidirectional_operator.

confidence:           medium
action:               reverse-mode refactor candidate, after the tower
```
