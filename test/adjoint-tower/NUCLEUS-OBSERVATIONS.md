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

                        J_Q = (I_Y o R_dep) o I_Q^T   <= Y x Q
                        J_P = (I_Y o R_dep) o I_P^T   <= Y x P
                        F_Q = (I_Z o R_dep) o I_Q^T   <= Z x Q
                        F_P = (I_Z o R_dep) o I_P^T   <= Z x P

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

                        C_Q = J_Q^T o J_Q = {(u,u),(u,v),(v,u),(v,v)}

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
