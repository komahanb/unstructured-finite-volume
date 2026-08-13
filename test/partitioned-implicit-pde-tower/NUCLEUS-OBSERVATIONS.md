# Nucleus observations — Partitioned Implicit PDE Tower

The live evidence ledger of this orbital client, per the discipline of
[FRACTAL-GRAPH-ARCHITECTURE.md](../../FRACTAL-GRAPH-ARCHITECTURE.md):

\[
\boxed{\text{observation}\neq\text{immediate refactor}}
\]

*Local necessity* means necessity at the object's own radius; *global
necessity* means necessity at larger contextual radius. Cross-tower
recurrence is its own field and is never used as proof of global
necessity. **The three gates of this tower are one client, not three.**

This tower was named by
[the reverse architecture review](../../doc/REVERSE-ARCHITECTURE-REVIEW.md)
as the discriminating experiment: the first deliberate **radius-2**
client, the first outside the derivative family, and the first that will
put a topology-consuming operation behind a minimizer.

---

## OBSERVATION PIP-1

```text
tower:                Partitioned Implicit PDE
gate:                 A
contextual radius:    2 (an object inside a partition of a larger whole)

symptom / fact:       PRESENCE IS NOT OWNERSHIP. Global vertex 3 is
                      present in both parts - owned by part 1, borrowed
                      by part 2 - and global vertex 4 likewise the
                      other way round. Global edge e3 = 3->4 is
                      present in both parts and owned by exactly one.

                      A part carrier therefore holds members it does
                      not answer for, and no count of members can tell
                      you which is which. Only vertex_owner_part /
                      edge_owner_part can.

exact caller:         gate-a-partition/test.f90
                      (check_part_structure,
                      check_crossing_edge_is_present_twice)

mathematical concept: overlap decomposition; ownership as a function
                      on local members

local necessity:      yes - assembly would double-count without it
global necessity:     yes at this radius and above; a partition with
                      no ownership function cannot reconstruct a whole

cross-tower recurrence: FIRST tower to meet it. Every prior tower's
                      domains were global and undivided

graph role:           partition frame - the part graph carries the
                      ownership function and the global map

comparison:           the adjoint tower distinguished DOMAINS of equal
                      size by identity (AD-0); this tower distinguishes
                      MEMBERS of one domain by ownership. Different
                      question, same discipline: a count never settles
                      identity

suspected nucleus implication: none - the production part relation
                      already carries exactly this and needed nothing

confidence:           high
action:               observe
```

---

## OBSERVATION PIP-2

```text
tower:                Partitioned Implicit PDE
gate:                 A
contextual radius:    2

symptom / fact:       OVERLAP IS VISIBILITY. A full global vertex field
                      partitions to a FULL field on each part's whole
                      vertex carrier - q1.domain same_as
                      G1.vertex_set(), four members including the
                      borrowed one - not to a field on the owned
                      subset.

                      That is the shape a stencil needs: owned vertex 3
                      cannot be evaluated without seeing q(4).

exact caller:         gate-a-partition/test.f90
                      (check_vertex_transport, check_one_transport)

mathematical concept: V_part = owned union borrowed as the numerical
                      INPUT carrier

local necessity:      yes - PROVED NUMERICALLY at Gate B, not merely
                      predicted: perturbing ONLY the borrowed copy
                      moves an OWNED result.
                        G1: q(borrowed global 4) += 10
                            -> owned A at global 3:  7 -> -3
                        G2: q(borrowed global 3) += 10
                            -> owned A at global 4: 13 ->  3
                      Restore the halo and the correct answers
                      return. Structural visibility and numerical
                      dependence are different claims, and this is
                      the second one
global necessity:     unknown at larger radius (deeper stencils need
                      wider overlap; this tower tests one ring only)

cross-tower recurrence: first tower to need it

graph role:           the part graph defines the visibility, and the
                      field's domain is the part's carrier

comparison:           learning L5 established that a field's domain is
                      a member_set and roles are domains; here the SAME
                      contract carries a subtler distinction - the
                      domain is the overlap, while the contribution
                      domain will be the owned subset

suspected nucleus implication: none. Recorded because the READ domain
                      and the WRITE domain of a partitioned operation
                      differ, which is the tower's central law and is
                      easy to get wrong in either direction.

confidence:           high
action:               observe
```

---

## OBSERVATION PIP-3

```text
tower:                Partitioned Implicit PDE
gate:                 A
contextual radius:    2

symptom / fact:       ONE GLOBAL ENTITY, ONE ASSEMBLED CONTRIBUTION.
                      The law was imposed BEFORE any ownership array
                      was inspected: an edge probe [10,20,30,40,50] was
                      partitioned to both parts, each part's field
                      assembled home, and the two summed. The result is
                      the probe exactly - the crossing edge contributes
                      30 once.

                      Only then was edge_owner_part consulted, and it
                      documents part 1 (the TAIL's owner) as canonical,
                      uniformly for every edge.

                      A PROSE DEFECT was found beside that code. The
                      comment claimed an edge is owned by its tail's
                      part "unless the tail is borrowed, in which case
                      the head's owner answers for it" - but both arms
                      of the if/else assign the tail's owner. The
                      behaviour is unconditionally tail-owned, which
                      satisfies the law; the comment described a rule
                      that was never implemented.

exact caller:         gate-a-partition/test.f90
                      (check_edge_assembly_law,
                      document_canonical_edge_owner);
                      src/class_graph_partitioner.f90 (the eowner
                      assignment)

mathematical concept: a partition of a set of entities into ownership
                      classes; assembly as the inverse of restriction

local necessity:      yes
global necessity:     yes - reconstruction is meaningless without it

cross-tower recurrence: first tower to test it on a CROSSING entity

graph role:           not applicable

comparison:           the law is what makes ownership canonical. Had
                      the tower begun by reading the comment, it would
                      have "fixed" working code to match wrong prose

suspected nucleus implication: PROSE ONLY, and applied: the comment now
                      describes tail-ownership and names the vestigial
                      branch as vestigial. The dead if/else was NOT
                      removed - RED exposed no defect, so the code was
                      left alone. It is recorded here as a hazard: both
                      arms are identical, so a future editor could
                      "repair" one of them and silently break unique
                      reconstruction. Removing it is a candidate for a
                      later cleanup pass, not for this tower.

confidence:           high
action:               observe (prose corrected; dead branch flagged)
```

---

## OBSERVATION PIP-7

```text
tower:                Partitioned Implicit PDE
gate:                 A
contextual radius:    2

symptom / fact:       LOCAL NUMBERING IS A CONTEXTUAL COORDINATE, not a
                      shared one. G2's local order is [4,5,6,3]: its
                      local member 1 is global vertex 4, and its local
                      member 4 is the borrowed global vertex 3.

                      No value in this tower can be read positionally
                      across the global/local boundary. Every read goes
                      through global_vertex_index and local_index, and
                      a proper subset S = {6,3,4} declared in
                      non-global order survives partition and
                      reassembly with its values intact.

exact caller:         gate-a-partition/test.f90 (every value check;
                      check_proper_subset_transport in particular)

mathematical concept: a domain's enumeration is storage, not identity

local necessity:      yes - the maps exist precisely because the
                      numbering differs
global necessity:     unknown, but the pressure only grows: more parts
                      means more independent numberings

cross-tower recurrence: THIRD distinct form of the same discipline -
                      learning/derivative pinned declaration order
                      within one domain; adjoint permuted TWO domains
                      independently; this tower has two carriers with
                      genuinely different numbering systems related by
                      a map. The lesson has now appeared in three
                      independent clients in three different shapes

graph role:           the part graph owns the map

comparison:           the adjoint tower's permutation test was
                      artificial hostility, deliberately constructed.
                      Here the mismatch is NATURAL - it is what
                      partitioning does - which is stronger evidence
                      that the discipline is not merely defensive

suspected nucleus implication: none - the maps are production and they
                      work. Recorded as the first NATURAL occurrence of
                      a hazard the earlier towers had to manufacture.

confidence:           high
action:               observe
```

---

---

## OBSERVATION PIP-4

```text
tower:                Partitioned Implicit PDE
gate:                 B
contextual radius:    2

symptom / fact:       THE GRAPH HOST IS A REAL CONDUIT, and it is
                      proved BEHAVIOURALLY rather than inferred from
                      production call sites.

                      One operation type (shifted_laplacian) that
                      stores NO graph; two gmres instances; two hosts
                      with the SAME six vertices and SAME five edges
                      but different topology - the chain and a star;
                      one numerical probe q*:

                        solver_G   % matvec(q*)  =  b
                        solver_alt % matvec(q*)  /= b

                      Nothing changed but the graph the solver carries.
                      If the host were scenery the two would agree.

                      The finding is NOT "GMRES traverses the graph" -
                      it does not; it reads no topology at all. It is
                      that GMRES CARRIES the graph to an attached
                      operation which does. Two roles:

                        minimizer     graph as conduit / context carrier
                        differential  graph as numerical topology
                        operation     operand

exact caller:         gate-b-operator/test.f90
                      (check_host_is_load_bearing); the chain reaches
                      production laplacian via
                      common/shifted_laplacian_fixture.f90

mathematical concept: an operator parameterized by the structure it
                      is applied over

local necessity:      YES for the attached action; the minimizer
                      itself still reads nothing
global necessity:     yes at this radius - and the reverse review's
                      Case-III caution is now settled empirically for
                      Class-1 clients

cross-tower recurrence: the FOUR sealed towers reported the host as
                      scenery (calculator, learning, adjoint) or took
                      no host at all (derivative action). Their
                      observation is not overturned - it was true of
                      THEIR actions, all Class-2/3. What is overturned
                      is the generalization. Independent tower count
                      for "host unused" does NOT change; its
                      INTERPRETATION does

graph role:           conduit (minimizer) and topology operand
                      (differential operator) in one call chain

comparison:           the reverse review inferred this from
                      test/graph-minimization/test.f90:345 and
                      class_graph_multigrid.f90:123. This tower
                      measures it: same operation, same probe,
                      different host, different answer

suspected nucleus implication: the reverse review's REJECT of
                      "remove class(graph) from graph_operation" is
                      now supported by evidence rather than by
                      inference. Removing it would break this call
                      chain outright. Seam A1 is CLOSED as a refactor
                      candidate.

confidence:           high
action:               observe - and close A1
```

---

## OBSERVATION PIP-5

```text
tower:                Partitioned Implicit PDE
gate:                 B
contextual radius:    2

symptom / fact:       THE PART GRAPH IS A GENUINE NUMERICAL OPERAND.
                      The production Laplacian is applied to G1 and G2
                      directly and traverses each part's own local
                      incidence, producing local answers that differ
                      from the global ones exactly where the part's
                      topology is incomplete:

                        G1 borrowed global 4:  local L = -3, global 1
                        G2 borrowed global 3:  local L =  3, global 1

                      So one object - the part graph - carries four
                      roles at once: partition frame, ownership
                      environment, local field-domain source, and
                      numerical topology operand.

exact caller:         gate-b-operator/test.f90 (check_local_actions,
                      one_local_action)

mathematical concept: restriction of an operator to a subdomain with
                      overlap

local necessity:      yes
global necessity:     unknown beyond this radius

cross-tower recurrence: first tower in which a graph is BOTH a frame
                      and an operand. In the adjoint tower the model
                      graph was an owner and never an operand
                      (AD-13); here the same object is both

graph role:           all four, simultaneously

comparison:           the borrowed rows are the evidence: a part's
                      answer is authoritative only where its topology
                      is complete, which is precisely what ownership
                      records

suspected nucleus implication: none - production needed nothing. The
                      observation matters for what it forbids: a
                      future "partitioned operator" abstraction must
                      not assume a part's output is globally valid.

confidence:           high
action:               observe
```

---

## OBSERVATION PIP-6

```text
tower:                Partitioned Implicit PDE
gate:                 B
contextual radius:    2

symptom / fact:       UNRESOLVED - and deliberately so.

                      The reverse review's seam A2 is that Class-2
                      operations (fit, reduction, broadcast,
                      difference_linearization) obtain their DOMAIN
                      from a graph's vertex_set() instead of carrying
                      it. Gate B was watched for natural pressure of
                      that kind.

                      None arose. Every operation this gate needed was
                      Class-1: it wanted the topology, not merely a
                      domain, and asking the graph for its vertex set
                      was therefore not a workaround but the correct
                      thing to do.

                      No Class-2 caller was imported to manufacture a
                      vote. difference_linearization, graph reduction
                      and broadcast are absent from this tower.

exact caller:         none - this entry records an ABSENCE

mathematical concept: n/a

local necessity:      n/a
global necessity:     n/a

cross-tower recurrence: seam A2 remains at TWO towers (derivative
                      action, adjoint), both in the derivative family.
                      This tower adds nothing to it

graph role:           n/a

comparison:           a tower should answer its own problem, not
                      collect evidence for someone else's. Recording
                      the absence honestly is the result

suspected nucleus implication: A2 stays in BUILD-ANOTHER-TOWER-FIRST.
                      A partitioned client did not supply the missing
                      independent vote; some other client must.

confidence:           high (that no pressure arose)
action:               observe; A2 still unresolved
```

---

## OBSERVATION PIP-8

```text
tower:                Partitioned Implicit PDE
gate:                 C
contextual radius:    2

symptom / fact:       A PARTITIONED MATRIX-FREE SOLVE WORKS, AND IT IS
                      NOT DISTRIBUTED. Production GMRES solved
                      A q = b through a composite that decomposes,
                      acts locally on each part's own topology, and
                      reassembles owned outputs - reaching the same q*
                      as the global road, and agreeing with the global
                      action at solver % matvec on four probes
                      including both interface basis vectors.

                      What remains CENTRALIZED, exactly:

                        one process, one global trial vector
                        direct global access inside partition_data
                        parts executed sequentially, in one address
                          space
                        assembly performed in-process
                        inner products over a global array
                        norms over a global array
                        Krylov state held globally

                      So the honest names are PARTITIONED MATRIX-FREE
                      SOLVE and serial semantic model of a partitioned
                      operator. The words distributed GMRES, MPI
                      solver, parallel halo exchange and distributed
                      vector are NOT used anywhere in this tower.

exact caller:         common/partitioned_shifted_laplacian_fixture.f90
                      (part_apply); gate-c-statement/test.f90

mathematical concept: operator decomposition vs execution
                      distribution - two independent axes

local necessity:      n/a - this entry records a BOUNDARY, not a need
global necessity:     n/a

cross-tower recurrence: first tower to reach radius 2 at all

graph role:           G1/G2 as topology operands (PIP-5); G as the
                      recorded decomposition context

comparison:           the SEMANTIC decomposition is complete: nothing
                      in the operator's mathematics assumes a shared
                      address space. What is missing is entirely
                      EXECUTION machinery. That is a useful thing to
                      have learned separately

DERIVED frontier - what a genuinely distributed road would need, in
the order the current code path would demand it:

                        1. an owned distributed field/vector - today
                           partition_data reads a global array that a
                           rank would not possess
                        2. overlap refresh as COMMUNICATION rather than
                           as a read: neighbour exchange of borrowed
                           values, replacing step (2) of part_apply
                        3. a distributed assembly that sums owned
                           contributions without materializing a global
                           array
                        4. a global inner product - the minimizer's
                           inner_product currently reduces over one
                           array
                        5. a global norm, likewise
                        6. Krylov synchronization: orthogonalization
                           over distributed vectors, and the collective
                           points that implies

                      Items 4-6 sit in graph_minimization, NOT in this
                      tower's fixtures - which is itself informative:
                      the OPERATOR decomposed cleanly with no
                      production change, while the SOLVER is where
                      distribution would actually bite.

                      These are FRONTIER observations. No MPI API is
                      designed here, no distributed_graph is
                      introduced, and none should be until a tower
                      genuinely needs one.

confidence:           high (the boundary); medium (the ordering of the
                      frontier items)
action:               observe; candidate subject for a future tower
```

---

## Reverse-evidence status after this tower

```text
A1  graph_operation host          CLOSED at Gate B - measured, not
                                  inferred. Gate C reinforces it and
                                  does not reopen it.

A2  operations carrying their     UNRESOLVED. Still TWO clients, both
    own domain                    derivative-family. This tower was the
                                  candidate for an independent third
                                  vote and honestly did not supply one
                                  (PIP-6).

B   bidirectional linearization   NO new vote. This tower is not a
                                  derivative client and needed no
                                  reverse action.

partition semantics               strong evidence ACROSS INCREASING
                                  RADIUS within ONE tower. Gates A, B
                                  and C are one client, not three, and
                                  are counted as one.
```

Nothing in this tower justifies a nucleus refactor. Its most useful
contribution to the reverse ledger is a rejection made safe: the seam
with the most tower-votes (A1) is now closed on measurement rather than
on inference.
