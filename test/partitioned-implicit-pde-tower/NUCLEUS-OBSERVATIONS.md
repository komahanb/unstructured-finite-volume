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

local necessity:      yes - Gate B's stencil will starve without it
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

## Deferred to later gates

```text
PIP-4   graph host as a real conduit to a Class-1 consumer   Gate B
PIP-5   the part graph as a genuine numerical operand        Gate B
PIP-6   whether Class-2 domain-from-graph pressure arises    Gate B
PIP-8   serial partitioned action vs distributed execution   Gate C
```

These are named now so that Gate B and Gate C are answering questions
that were asked in advance, rather than reporting whatever they happen
to notice.
