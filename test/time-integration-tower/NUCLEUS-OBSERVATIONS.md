# Nucleus observations — Time Integration Tower

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
not ten** — and Levels 5–9 do not exist yet, so nothing below claims
evidence from them.

---

## A standing caution for this ledger

Every entry at Gate A is `action: observe`, and **none of them is
evidence for seam A2** of
[the reverse architecture review](../../doc/REVERSE-ARCHITECTURE-REVIEW.md)
— *operations should carry their domain rather than ask a graph for
one*.

A2 concerns what an **operation** does. Gate A has attached no
operation to anything. It has not asked a graph for a domain, because it
has nothing that would need one. The seam is untouched here, and
recording these six structural observations as though they bore on it
would be counting a tower that has not run the experiment.

---

## OBSERVATION TI-1

```text
tower:                Time Integration
level:                0  (carrier)
review gate:          A
contextual radius:    0 (three declared sets, nothing between them)

symptom / fact:       Q AND T ARE INDEPENDENT CARRIERS, and the tower
                      declares both before anything can conflate them.

                        Q = { x, y }                state coordinates
                        T = { t0 ... t4 }           time instants
                        E = { e1 ... e4 }           time steps

                      Q % same_as(T) is false, and so are the other
                      two pairings. All three enumerate from one, so
                      the integer 1 is a member of each - meaning x,
                      t0 and e1 - and no count and no numeral
                      separates them.

                      Q is declared at Level 0 deliberately, before
                      any value exists to live on it, because the
                      collapse this tower refuses

                        state coordinate ~ time instant ~ vertex

                      has to be refused before there is anything to
                      collapse them WITH.

exact caller:         level-0-carrier/test.f90
                      (check_identities, check_boundaries)

mathematical concept: a carrier is an identity, not a cardinality;
                      two axes of one problem are two objects

local necessity:      yes - nothing else distinguishes the sets
global necessity:     unknown at this radius. The claim that KEEPING
                      them apart is load-bearing is a Level-5+
                      question, not settled here

cross-tower recurrence: SECOND independent form. The partitioned
                      tower met it on V/E/K, where all three belonged
                      to one mesh. Here the two axes belong to
                      DIFFERENT mathematical questions - state and
                      time - which is a stronger separation and a
                      weaker excuse for confusing them

graph role:           none yet; no graph exists at this level

comparison:           the partitioned tower's L0 proved the same
                      discipline on carriers that were all spatial.
                      This is the first tower whose carriers span two
                      independent axes

suspected nucleus implication: none. counted_set already carries
                      identity and needed nothing.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-2

```text
tower:                Time Integration
level:                1  (relation)
review gate:          A
contextual radius:    0

symptom / fact:       TEMPORAL DIRECTION IS RELATION STRUCTURE, NOT
                      LOOP INDEX ORDER.

                        Tail <= E x T    e_i -> t_(i-1)
                        Head <= E x T    e_i -> t_i

                      A step knows which end it leaves because Tail
                      says so and Head says otherwise - not because
                      one integer exceeds another. Pinned: every step
                      has exactly one tail and exactly one head, and
                      NO step leaves the instant it enters. Had the
                      two agreed anywhere, that step would be a loop
                      and time would stand still on it.

                      Q participates in no relation at this level.
                      That absence is the level's content, checked
                      positively: no slot of Tail or Head answers Q.

                      A numbering hazard worth the ledger's space:
                      over this specimen Tail's extension is
                      { [1,1] [2,2] [3,3] [4,4] }, tuple for tuple
                      what a six-vertex chain's tail map looked like
                      in the partitioned tower over entirely
                      different carriers. The integers carry none of
                      the meaning; the signature carries all of it.

exact caller:         level-1-relation/test.f90
                      (check_direction_is_structure,
                      check_state_axis_is_untouched)

mathematical concept: orientation as a pair of maps from an incidence
                      carrier, rather than as an ordering on values

local necessity:      yes
global necessity:     unknown - production's marcher currently
                      expresses time as a loop, and whether that is
                      equivalent is a Gate-B question

cross-tower recurrence: first tower to state TIME this way

graph role:           none yet

comparison:           the partitioned tower derived spatial adjacency
                      from the same two primitives. That the identical
                      shape carries temporal direction is either a
                      real generality or a coincidence of specimens;
                      one tower cannot tell

suspected nucleus implication: none. The observation is about what
                      this client SAYS, not about what the nucleus
                      lacks.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-3

```text
tower:                Time Integration
level:                2  (relation algebra)
review gate:          A
contextual radius:    0

symptom / fact:       ONE-STEP AND TWO-STEP REACH ARE GENERATED
                      RELATIONS, not stored ones.

                        A1 = Head o Tail^T : T -> T
                             { t0->t1, t1->t2, t2->t3, t3->t4 }
                        A2 = A1 o A1       : T -> T
                             { t0->t2, t1->t3, t2->t4 }

                      Neither is written down as data. A2 has THREE
                      facts, not four - two steps do not fit from t3 -
                      and that arithmetic falls out of composition
                      rather than being maintained by anyone.

                      A1 and A2 differ by IDENTITY and by EXTENSION
                      both, though they share a signature. A shared
                      signature has never been an address.

                      AND THE BOUNDARY THIS ENTRY EXISTS TO HOLD:

                        A2 IS NOT BDF2.

                      A2 says only that an instant two steps later is
                      structurally reachable. The executable form of
                      the distinction is a refusal - A2 does NOT
                      relate t0 to t3 - because two-step reach lands
                      only on even offsets, while a two-step SCHEME
                      would have to consume the history at t1 and t3.
                      Whatever BDF2 later consumes, it must
                      CONSTITUTE out of this structure; it will not
                      find itself lying here.

                        temporal reach != temporal discretization

exact caller:         level-2-relation-algebra/test.f90
                      (check_two_step_reach, check_the_two_differ,
                      check_reach_is_not_a_scheme)

mathematical concept: composition of a relation with itself; reach
                      versus dependency

local necessity:      yes
global necessity:     unknown - whether a scheme can be constituted
                      from these two is exactly the Level-6 question

cross-tower recurrence: the calculator tower derived a dependency by
                      restrict/project/compose; this is the first
                      client to compose a relation WITH ITSELF, and
                      the first to need a name for what the result is
                      NOT

graph role:           none yet

comparison:           no union and no transitive closure were added.
                      Writing A1 U A2 because the notation exists
                      would be inventing algebra for a reader rather
                      than for a caller, and the nucleus has held that
                      line since the calculator tower

suspected nucleus implication: none. compose_binary was sufficient,
                      including composing a relation with itself,
                      which no prior caller had asked for.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-4

```text
tower:                Time Integration
level:                3  (relational graph)
review gate:          A
contextual radius:    1 (a container holding both axes)

symptom / fact:       THE RELATIONAL GRAPH MAY OWN A STATE CARRIER
                      THAT NO RELATION NAMES, and this is lawful
                      rather than tolerated.

                        G_time = ( {Q,T,E}, {Tail,Head,A1,A2} )

                      Q is owned. No slot of any owned relation
                      resolves to Q - checked positively, by walking
                      every relation's every slot through the graph's
                      own accessors.

                      The law is ASYMMETRIC, and that asymmetry is
                      the finding. create_graph validates RELATIONS
                      against SETS: every slot must answer an owned
                      carrier. It never validates sets against
                      relations, and it should not - a relational
                      graph is a collection of member sets and typed
                      relations over them, not a connected object.

                      So the state axis waits inside the same
                      structure as the time axis, answering nothing,
                      until Level 5 puts a field on it. The
                      alternative - inventing an incidence to attach
                      Q to the chain - would manufacture the exact
                      conflation this tower refuses, in the name of
                      tidiness.

                        Q is unrelated here.  Unrelated is not absent.

exact caller:         level-3-graph/test.f90
                      (check_state_axis_is_owned_and_unrelated,
                      check_signature_closure)

mathematical concept: a relational structure need not be connected;
                      signature closure is a one-way obligation

local necessity:      yes - the alternative is a fabricated relation
global necessity:     unknown, but the pressure is legible: any client
                      with two independent axes will meet it

cross-tower recurrence: FIRST tower to own a carrier no relation
                      names. Calculator, Learning, Derivative Action
                      and Adjoint all named every carrier they owned

graph role:           structural owner of two independent axes

comparison:           the reverse review's A3 found relational_graph
                      load-bearing as an owner across three towers and
                      recommended KEEP. This is a fourth, and a new
                      shape of the same finding: it holds correctly
                      even when part of what it owns is, so far,
                      inert

suspected nucleus implication: none, and specifically NOT a request
                      for a "disconnected graph" refusal. The
                      container is right as it stands; this entry
                      records that its asymmetry is deliberate and now
                      has a client that depends on it.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-5

```text
tower:                Time Integration
level:                4  (graph calculus)
review gate:          A
contextual radius:    1

symptom / fact:       CAUSALITY IS A GRAPH-CALCULUS INTERPRETATION OF
                      A1, not a property A1 carries.

                        sources(A1)           = { t0 }
                        sinks(A1)             = { t4 }
                        reachable(t0, t4)     = true
                        reachable(t4, t0)     = false
                        topological_order(A1) = [t0 t1 t2 t3 t4]

                      Nothing at Level 2 called t0 "first". A1 is a
                      relation; it does not know it is time. What
                      makes this list a CAUSAL order is the reading
                      imposed by a profile and four algorithms that
                      would say the same things about a dependency
                      between calculator operations.

                      The order comes from the carrier's declaration
                      order through local_index, never from
                      arithmetic on member values - and the test
                      checks it against t % member(i) rather than
                      against a written [1,2,3,4,5], because here the
                      two coincide and a literal would hide the
                      dependency that makes the coincidence
                      meaningful.

                      Also pinned: A1 and A2 are DIFFERENT STRUCTURAL
                      VIEWS OVER THE SAME T CARRIER - both
                      view % domain() answer T by identity. One time
                      axis, several readings of it, which is the
                      Rosetta truth a later scheme can consume.

exact caller:         level-4-graph-calculus/test.f90
                      (check_causal_ends, check_causal_reachability,
                      check_forward_causal_order,
                      check_two_views_one_carrier)

mathematical concept: interpretation as a separable layer above
                      structure; topological order as causal order

local necessity:      yes
global necessity:     unknown - whether production's marcher agrees
                      with this order is a Gate-B question, and the
                      marcher is deliberately not imported here

cross-tower recurrence: the calculator tower read a dependency this
                      way at its own Level 4, and the adjoint tower
                      after it. This is the third client, and the
                      first whose subject is TIME - the profile and
                      the four algorithms needed no change to say
                      something true about it

graph role:           the graph owns A1 and A2; the views BORROW
                      them. The selectors are separate objects
                      carrying only the relations' identities, and
                      they are destroyed before any question is asked
                      of the views - so the borrow is demonstrated,
                      not described

comparison:           the calculator tower killed its selector the
                      same way at its own Level 4. That the idiom
                      transfers unchanged to a client with a wholly
                      different subject is evidence the ownership
                      contract is general, which supports A3's KEEP

suspected nucleus implication: none. graph_profile and
                      graph_algorithms answered a subject they were
                      not designed for - time - without modification.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-6

```text
tower:                Time Integration
level:                4  (graph calculus)
review gate:          A
contextual radius:    1

symptom / fact:       REVERSE CAUSAL ORDER EXISTS STRUCTURALLY
                      BEFORE AN ADJOINT EXISTS.

                        forward   [t0 t1 t2 t3 t4]
                        reverse   [t4 t3 t2 t1 t0]

                      The reverse order is the forward order read
                      backwards, and that is the whole of it. No new
                      algorithm was written. No adjoint was
                      implemented. No derivative or adjoint module is
                      imported at this level or any level below it,
                      and the import gate refuses them all.

                      The observation is worth a number because of
                      what it implies about sequencing: whoever later
                      builds an adjoint through time will find the
                      ordering ALREADY HERE, established by causal
                      structure, and will not have to invent it
                      alongside the mathematics. Reverse traversal is
                      a property of the dependency; the adjoint is a
                      separate thing that happens to need it.

exact caller:         level-4-graph-calculus/test.f90
                      (check_forward_causal_order, final assertion)

mathematical concept: the reverse of a topological order on a DAG;
                      traversal direction as distinct from the
                      mathematics traversing it

local necessity:      yes for the statement; nothing consumes it yet
global necessity:     unknown. The adjoint tower built reverse-mode
                      sensitivity WITHOUT a time dimension, so it has
                      no opinion on this

cross-tower recurrence: first tower to separate reverse ORDER from
                      reverse-mode MATHEMATICS. The adjoint tower
                      needed the second and never had the first

graph role:           the order is read off the A1 view; nothing
                      further is asked of the graph

comparison:           this entry deliberately claims LESS than it
                      could. It would be easy to write "the nucleus
                      already supports time-adjoint traversal". It
                      does not: no adjoint has been attached to
                      anything here, and one tower stating an ordering
                      is not evidence that a reverse-mode calculation
                      through it would work

suspected nucleus implication: none yet. Recorded so that a future
                      time-adjoint client can point at where its
                      traversal order came from.

confidence:           medium — the FACT is certain; its usefulness is
                      a prediction about levels that do not exist
action:               observe
```

---

# Frontier

```text
Levels 5-9 are UNBUILT.

WHAT GATE B WILL ASK

    When q becomes a field on Q, can the existing temporal
    discretization and minimization stack preserve Q as the state
    domain, independently of whatever graph supplies structural
    context?

WHAT WILL BE TESTED THERE, recorded neutrally

    the current step / march code - class_graph_step,
    class_graph_marcher - is deliberately NOT imported at Gate A,
    and will be exercised at Gate B.

WHAT THIS LEDGER DOES NOT SAY

    that the production marcher is wrong
    that SEAM A2 must be closed
    that graph_operation must change
    that seam A2 has gained a third independent tower

A NOTE ON THE LABEL, because this ledger uses it both ways:

    SEAM A2   the reverse review's "operations should carry their
              own domain rather than ask a graph for one"
    A2        this tower's two-step reach relation, A1 o A1

They share two characters and nothing else. Every entry above means
the RELATION; the standing caution and this block mean the SEAM.

None of those is established by Levels 0-4, and the tower will be
worth less as evidence if it claims them early.
```
