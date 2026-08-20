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
not ten.** The tower is sealed: all ten levels are built.

---

## A standing caution for this ledger

**TI-1 … TI-6 are Gate-A structural observations, and none of them is
evidence for seam A2** of
[the reverse architecture review](../../doc/REVERSE-ARCHITECTURE-REVIEW.md)
— *operations should carry their domain rather than ask a graph for
one*. Gate A attached no operation to anything, so it could not have
exercised the seam.

**TI-8 is the seam-A2 observation.** It was earned at Level 6, by an
experiment that failed first and was recorded before production was
touched. It is ONE tower's evidence, not three levels' worth: Levels
0–9 of this tower are ONE client, and the count below moves by one.

**This tower is ONE client, and it produced THREE REDs.**
`class_graph_step` at Level 6 (TI-8), `class_graph_marcher` at Level 8
(TI-14), and `class_graph_linearization` at Level 8 (TI-16). They are
one tower's evidence, not three. **Time contributes ONE tower vote to
seam A2**, and the count moves from 2 to 3 — never to 5.

What the three REDs increase is evidence *quality*, not count: they
show the same seam in a discretization, in a march, and in a
finite-difference linearization, reached through one production call
chain. That is recorded in TI-16 and in the seal summary at the end.

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
                      equivalent remains UNANSWERED: Gate B solved
                      single steps by hand and never imported the
                      marcher

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
                      relation_algorithms answered a subject they were
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

## OBSERVATION TI-7

```text
tower:                Time Integration
level:                5  (field calculus)
review gate:          B
contextual radius:    0

symptom / fact:       THE STATE FIELD HAS ITS OWN DOMAIN, and no
                      graph exists anywhere in the Level-5 program.

                        q0   : Q -> reals [2, 0]
                        time : T -> reals [0, 1/2, 1, 3/2, 2]
                        h    : E -> reals [1/2, 1/2, 1/2, 1/2]

                      Each pinned by domain IDENTITY, and the three
                      domains pairwise distinct. q0 lives on Q
                      because a field is a function over one
                      member_set - which is what production field
                      calculus already says, in those words.

                      This is not a new capability. It is an
                      EXERCISE of an existing one by a client whose
                      state domain is emphatically nobody's vertex
                      set, and the exercise is what makes Level 6's
                      question askable at all: if q0 had needed a
                      graph to exist, the seam would have been
                      invisible one rung later.

                      Values are not structure. The instant t2 is a
                      MEMBER of T; the real 1.0 is a VALUE at t2.
                      Their consistency is proved rather than
                      assumed, over the relations Level 1 earned:

                        time(head(e)) - time(tail(e)) = h(e)

                      for every step - field calculus over
                      established structure, and not a time scheme.

exact caller:         level-5-field-calculus/test.f90
                      (check_domains_by_identity,
                      check_coordinates_agree_with_structure)

mathematical concept: a field is a function over a member_set; the
                      domain is the mathematics, the container is not

local necessity:      yes
global necessity:     unknown at this radius

cross-tower recurrence: every prior tower put fields on domains too.
                      What is new is that this client's state domain
                      could not be mistaken for a vertex set even by
                      accident - |Q| = 2 against a 5-vertex host one
                      level later

graph role:           NONE. No graph is constructed at this level,
                      and that absence is the observation

comparison:           the learning tower established that a field's
                      domain is a member_set and roles are domains.
                      This client needed that to already be true, and
                      found it was

suspected nucleus implication: none. field_stored needed
                      nothing; it was already right.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-8

```text
tower:                Time Integration
level:                6  (discretization)
review gate:          B
contextual radius:    1

symptom / fact:       *** THE SEAM-A2 OBSERVATION. RED FIRST. ***

                      TEMPORAL DISCRETIZATION SUBSTITUTED THE HOST'S
                      VERTICES FOR THE STATE DOMAIN.

                      The experiment: an action S : Q -> Q that
                      carries its own domain and stores no graph;
                      a compatibility host H_t with FIVE vertices;
                      |Q| = 2. The mismatch is deliberate - at equal
                      cardinality the substitution would have
                      produced plausible numbers and hidden.

                      The direct action was already correct and
                      needed no production change:

                        decay % domain(H_t)  ->  Q          PASS
                        decay % apply(H_t, q0) on Q         PASS
                        S(q0) = [2, -2]                     PASS

                      The step operator built FROM it was not. Run
                      against src as it stood at 6b60879 - BEFORE
                      the correction this commit carries - VERBATIM.
                      Every file:line below is a PRE-FIX line number
                      and will not match the current source:

                        FAIL : the backward-euler STEP answers Q when
                               asked its domain: TEMPORAL
                               DISCRETIZATION PRESERVES THE DOMAIN OF
                               THE ACTION IT DISCRETIZES
                        FAIL : and it does NOT answer the host's five
                               vertices - the state domain is not
                               inferred from the conduit

                        ERROR STOP field_stored: a value vector
                        must fill its domain exactly

                        Error termination. Backtrace:
                        #3 __class_graph_field_MOD_field_set_real_vector
                           at src/field_stored.f90:305
                        #4 __class_graph_step_MOD_step_apply
                           at src/class_graph_step.f90:200

                      Two distinct defects in one operator:

                        step_domain  answered
                                     input_graph % all_vertices(...)
                                     - a 5-member carrier for a
                                     2-member unknown

                        step_apply   took its width as
                                     size(y) / num_vertices() = 2/5 = 0
                                     and landed the answer on
                                     input_graph % vertex_set()

                      The residual of a step is a statement about the
                      same unknown the action is a statement about.
                      Neither the domain nor the width was ever the
                      host's business.

exact caller:         level-6-discretization/test.f90
                      (check_step_domain_is_the_action_s,
                      check_backward_euler_residual);
                      src/class_graph_step.f90 (step_domain:163,
                      step_apply:179/198/199 as they stood)

mathematical concept: an operator built from another operator
                      inherits that operator's domain, not its
                      container's

local necessity:      yes - the level cannot pass otherwise
global necessity:     YES for temporal discretization. A step is
                      ALWAYS built from an action, so the delegation
                      is right for every caller, not only this one

cross-tower recurrence: *** THIS IS THE THIRD INDEPENDENT TOWER FOR
                      SEAM A2, AND THE FIRST OUTSIDE THE DERIVATIVE
                      FAMILY. ***

                        before   Derivative Action, Adjoint
                                 = 2 towers, one family
                        after    + Time Integration
                                 = 3 towers, two families

                      which reaches the repository's strong-evidence
                      bar. Counted as ONE client: this tower's ten
                      levels are one client, not ten

graph role:           compatibility conduit. H_t is present because
                      the graph_operation contract requires a graph,
                      and its topology is never read by this action

comparison:           the reverse review predicted A2 would need a
                      client "outside the derivative family" whose
                      domains are not built by a graph's vertex set.
                      It proposed a partitioned or coupled tower.
                      Time integration turned out to be a cleaner
                      instrument: the state domain is not merely
                      DIFFERENT from the host's vertices, it is a
                      different KIND of thing, and no partitioning
                      relates them

suspected nucleus implication: APPLIED, and deliberately narrow -
                      see TI-9. Nothing broader was performed.

confidence:           high
action:               ACT (narrow), then observe
```

---

## OBSERVATION TI-9

```text
tower:                Time Integration
level:                6  (discretization)
review gate:          B
contextual radius:    1

symptom / fact:       DISCRETIZATION PRESERVES THE ACTION'S DOMAIN -
                      the correction TI-8's RED earned, and the whole
                      of it.

                      In src/class_graph_step.f90 ONLY:

                        step_domain   delegates to
                                      this % action % domain(...)

                        step_apply    takes ncomp from
                                      input_data(1) % num_components()
                                      validates the input's domain
                                      validates the action's answer
                                      lands the residual on the
                                      action's domain

                        nv            deleted, with both its uses

                      The rule, stated once:

                        TEMPORAL DISCRETIZATION PRESERVES THE DOMAIN
                        OF THE ACTION IT DISCRETIZES.

                      REGRESSION, verified rather than argued. Every
                      action on the ordinary-graph road answers
                      input_graph % all_vertices(domain) - checked
                      directly in mandelbrot_law, vdp_law,
                      vdp_tangent_law and vdp_adjoint_law - and the
                      test pins that all_vertices(H) and
                      H % vertex_set() are the same carrier. So
                      delegation returns exactly what asking the graph
                      returned, and test/graph-marching passes
                      unchanged, including its two-numbers-wide cell
                      case, which is the one that exercises ncomp.

                      class_graph_marcher is the ONLY consumer of
                      step_operator in the repository, and it was not
                      touched.

                      WHAT WAS NOT DONE, deliberately:

                        graph_operation root      unchanged
                        graph_grammar             unchanged
                        graph_fitting             unchanged
                        reduction / broadcast     unchanged
                        difference_linearization  unchanged
                        class_graph_marcher       unchanged

                      Strong evidence makes the broader A2 migration
                      ELIGIBLE for reverse review; it does not make
                      every possible implementation MANDATORY. None
                      of those four call sites was exercised by this
                      tower, and a refactor of an untested call site
                      is a speculation wearing evidence's clothes.

exact caller:         src/class_graph_step.f90 (step_domain,
                      step_apply); level-6-discretization/test.f90

mathematical concept: delegation of a domain question to the object
                      that owns the mathematics

local necessity:      yes
global necessity:     yes for this operator; UNKNOWN for the other
                      four A2 call sites, which remain untested

cross-tower recurrence: the correction is new; the seam is TI-8's

graph role:           unchanged - the host is still passed to the
                      action, which may consume its topology or not

comparison:           graph_minimization already had this contract
                      (TI-11). One module asked the host, its
                      immediate collaborator asked the action, and
                      only a client with |Q| != |V(H)| could tell
                      them apart

suspected nucleus implication: applied and complete for this
                      operator. The broader question is recorded for
                      reverse review, not answered.

confidence:           high
action:               ACT - applied, narrow, regression-verified
```

---

## OBSERVATION TI-10

```text
tower:                Time Integration
level:                6  (discretization)
review gate:          B
contextual radius:    1

symptom / fact:       STRUCTURE AND SCHEME, JOINED WITHOUT BEING
                      CONFUSED - the other half of TI-3.

                      Level 2 derived A1 and A2 and REFUSED to call
                      either one BDF2. Level 6 supplies what was
                      missing. At instant t2:

                        A1-predecessor of t2  =  t1
                        A2-predecessor of t2  =  t0

                      and bdf(2, ...) carries reach = 2 with

                        a0 = 3/2  at t2      the present
                        a1 = -2   at t1      one-step history
                        a2 = 1/2  at t0      two-step history

                      So:

                        A1, A2   supply STRUCTURAL REACH - which
                                 instants a two-step scheme may look
                                 at

                        bdf-2    supplies NUMERICAL COEFFICIENTS on
                                 exactly those roles

                      Neither contains the other. A2's extension is
                      fixed by composition and knows no coefficient;
                      the scheme's table is fixed by accuracy order
                      and knows no carrier. Level 2 was right to
                      refuse the name, and this is where the two
                      halves meet.

exact caller:         level-6-discretization/test.f90
                      (check_reach_supplies_the_history_roles)

mathematical concept: a discretization is a weighting of a
                      structurally available dependency

local necessity:      yes
global necessity:     unknown - one scheme family on one specimen

cross-tower recurrence: first tower to hold reach and scheme apart
                      and then join them explicitly

graph role:           none at this junction; the reach relations are
                      graph-owned but the scheme reads none of them

comparison:           TI-3 recorded the refusal; this records the
                      join. Together they are one argument in two
                      parts, and the parts are one level apart on
                      purpose

suspected nucleus implication: none. Recorded because "A2 is BDF2"
                      is the single most natural wrong sentence
                      anyone reading this tower could write.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-11

```text
tower:                Time Integration
level:                7  (minimization)
review gate:          B
contextual radius:    1

symptom / fact:       MINIMIZATION WAS ALREADY DOMAIN-EXPLICIT, and
                      needed no change whatever.

                      Production GMRES solved two implicit temporal
                      steps on Q while the host carried five
                      unrelated vertices:

                        backward euler  c = -q0 = [-2, 0]
                                        rhs = [2, 0]
                                        q1 = [4/3, 4/9]

                        bdf-2           c = -2q1 + q0/2
                                        rhs = [5/3, 8/9]
                                        q2 = [5/6, 47/72]

                      Both right-hand sides are fields ON Q; both
                      solutions come back ON Q; the affine constants
                      were measured, not assumed.

                      graph_minimization takes its unknown domain as
                      an EXPLICIT argument and asks the ACTION for
                      the residual domain. Its own comment says so:
                      "no hidden fallback to the host's vertices; a
                      caller that means vertices says so at its own
                      call site."

                      THE CONTRAST IS THE EVIDENCE. Two collaborating
                      modules, one seam:

                        class_graph_step      asked the HOST     RED
                        graph_minimization    asked the ACTION   fine

                      Seam A2 is therefore not a uniform defect in
                      the nucleus. It is a contract that some modules
                      already keep and others did not, and the ones
                      that keep it show what the correction should
                      look like rather than needing one.

exact caller:         level-7-minimization/test.f90
                      (check_backward_euler_solve, check_bdf2_solve,
                      check_unknown_domain_is_the_caller_s_word);
                      src/graph_minimization.f90 (attach:159-166)

mathematical concept: an unknown domain is the caller's declaration,
                      not the container's property

local necessity:      yes
global necessity:     yes - and already satisfied

cross-tower recurrence: the four sealed towers all attached
                      minimizers whose unknown domain HAPPENED to be
                      the host's vertex set. This is the first client
                      for which it demonstrably is not, so the first
                      to prove the explicit argument does real work

graph role:           conduit only. The solver passes H_t to the
                      step, which passes it to the action, which
                      ignores it

comparison:           seam A1 stays CLOSED. The partitioned tower
                      settled on production evidence that the host is
                      a real conduit for topology-consuming actions.
                      This action consumes no topology because a
                      triangular 2x2 decay has none - a property of
                      THIS action, not a counterexample to that
                      finding

suspected nucleus implication: NONE, and specifically no
                      "improvement" to graph_minimization. Nothing
                      here found anything to improve.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-12

```text
tower:                Time Integration
level:                8  (constitution)
review gate:          C
contextual radius:    1

symptom / fact:       THE PRODUCTION MARCHER'S CONTROL CHAIN IS
                      EXTENSIONALLY A1.

                      clock % instants(4) generates a five-vertex,
                      four-edge chain. Step i joins instant i to
                      instant i+1 - exactly the pairs A1 holds, read
                      through the TIME carrier's own members rather
                      than through the integers the chain happens to
                      use.

                      The bridge is EXTENSIONAL and deliberately not
                      identity:

                        chain % vertex_set()  same_as  T     NOT
                                                            required
                        incidence(chain)      =        A1    required

                      Two realizations of one structure. Demanding
                      they be the same object would be demanding
                      that two parties who agree must be the same
                      party - the discipline Level 3 of the
                      partitioned tower already established for an
                      ordinary graph against its relations.

exact caller:         level-8-constitution/test.f90
                      (check_control_chain_realizes_a1)

mathematical concept: two realizations of one relational structure;
                      extensional agreement without identity

local necessity:      yes
global necessity:     unknown - see TI-18. This holds because the
                      specimen's time graph IS a simple chain

cross-tower recurrence: the partitioned tower checked an ordinary
                      graph against Tail/Head the same way. This is
                      the second client, and the first where the
                      realization is generated by production rather
                      than declared by the test

graph role:           the control chain is a stored_graph the
                      marcher makes for itself

comparison:           Gate A recorded that production "regenerates a
                      linear chain from nsteps" as a neutral fact.
                      This measures it

suspected nucleus implication: none. Recorded because it is the
                      bridge that lets Levels 0-4 mean something
                      about production rather than beside it.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-13

```text
tower:                Time Integration
level:                8  (constitution)
review gate:          C
contextual radius:    1

symptom / fact:       THE OPERATION HOST IS NOT THE TIME GRAPH, and
                      at Level 8 there are THREE five-element objects
                      in one program to keep apart:

                        T                       the instant carrier
                        chain from instants(4)  the marcher's own
                                                control chain
                        V(H_context)            the operation host's
                                                vertices

                      plus the two-member Q, which is none of them.
                      No two are same_as, and the level pins all of
                      it.

                      THAT ALL THREE HAVE FIVE ELEMENTS IS A
                      COINCIDENCE OF THIS SPECIMEN. The assertions
                      refuse to lean on it, and the naming was
                      corrected at this level: what Gate B called
                      H_t is H_context, the compatibility conduit
                      the graph_operation contract requires - never
                      the clock.

exact caller:         level-8-constitution/test.f90
                      (check_the_three_carriers_stay_apart)

mathematical concept: distinct identities behind equal cardinalities

local necessity:      yes - the whole tower is an argument against
                      collapsing these
global necessity:     unknown, but the pressure is legible: any
                      client with a time axis and a state axis meets
                      it

cross-tower recurrence: the fourth distinct shape of "a count never
                      settles identity" - calculator/learning
                      (declaration order), adjoint (two permuted
                      domains), partitioned (two numbering systems),
                      and now THREE five-element objects that are
                      three different things

graph role:           two graphs in one program with different jobs:
                      one is the clock, one is the conduit

comparison:           Gate B chose a five-vertex host to make |Q|
                      differ from |V(H)|. At Level 8 that same
                      choice creates a NEW hazard - the host now
                      matches the control chain's size - and the
                      level pins them apart rather than renumbering
                      to avoid the question

suspected nucleus implication: none. Recorded as the naming
                      correction it forced.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-14

```text
tower:                Time Integration
level:                8  (constitution)
review gate:          C
contextual radius:    1

symptom / fact:       *** SEAM-A2 RED, SECOND SITE. RED FIRST. ***

                      THE MARCHER RECONSTRUCTED THE EVOLVING STATE
                      FROM THE HOST'S VERTICES.

                      Same experiment as TI-8, one level up: the
                      action carries its own domain, |Q| = 2, the
                      host has five vertices. All nine structural
                      assertions of Level 8 passed first - the
                      control-chain bridge and the carrier
                      distinctions needed no production change - and
                      then the first march died.

                      Run against src at f7fb641, VERBATIM:

                        ERROR STOP field_stored: a value vector
                        must fill its domain exactly

                        Error termination. Backtrace:
                        #3 __class_graph_field_MOD_field_set_real_vector
                           at src/field_stored.f90:305
                        #4 read_statement
                           at src/class_graph_marcher.f90:234
                        #5 __class_graph_marcher_MOD_march
                           at src/class_graph_marcher.f90:124

                      (PRE-FIX line numbers; they will not match the
                      current source.)

                      read_statement did what step_apply had done
                      one level down:

                        nv    = on % num_vertices()      = 5
                        ncomp = size(q) / max(nv, 1)     = 2/5 = 0
                        state = field(..., on % vertex_set(), ...)

                      A march is a repeated application of ONE
                      action, so the thing being marched inhabits
                      that action's domain. The host is the conduit
                      the action is reached through, and was never
                      the seat of the mathematics.

exact caller:         level-8-constitution/test.f90
                      (check_forward_march);
                      src/class_graph_marcher.f90
                      (read_statement:230-233, and the implicit
                      attach at 161-162, as they stood)

mathematical concept: a repeated operation evolves a state in the
                      operation's own domain

local necessity:      yes
global necessity:     YES for marching. A march is ALWAYS a repeated
                      action, so the delegation is right for every
                      caller

cross-tower recurrence: SAME TOWER, second site. This is NOT a
                      second tower vote - see the standing caution

graph role:           conduit, and still passed to every action

comparison:           TI-8 found the seam in a discretization; this
                      finds it in the thing that repeats the
                      discretization. Two sites, one contract, one
                      client

suspected nucleus implication: APPLIED, narrow - see TI-15.

confidence:           high
action:               ACT (narrow), then observe
```

---

## OBSERVATION TI-15

```text
tower:                Time Integration
level:                8  (constitution)
review gate:          C
contextual radius:    1

symptom / fact:       THE MARCHER PRESERVES THE ACTION'S DOMAIN -
                      the correction TI-14's RED earned.

                      In src/class_graph_marcher.f90 ONLY:

                        state_seat()  new, private: asks
                                      action % domain(on, ...),
                                      refuses an empty domain,
                                      refuses a state that does not
                                      divide evenly, and answers the
                                      seat and the width together

                        read_statement  builds the state on that
                                        seat; validates the action
                                        answered there; validates
                                        the answer's width

                        march (implicit) attaches the governed solve
                                        with unknown_domain =
                                        state_domain and ncomp from
                                        the same source

                      The rule:

                        A MARCH EVOLVES A STATE IN THE DOMAIN OF THE
                        ACTION IT MARCHES.

                      The equality check is by IDENTITY, not length:
                      q <- q - h s is an equation between two states,
                      and equal length would let a foreign carrier
                      through.

                      REGRESSION: every action on the ordinary-graph
                      road answers input_graph % all_vertices, and
                      Level 6 pins that all_vertices and vertex_set
                      are the same carrier. graph-marching was run
                      IMMEDIATELY after this change and passed
                      unchanged, including its two-numbers-wide cell.

                      WHAT WAS NOT DONE:

                        the public march signature   unchanged
                        the graph_operation root     unchanged
                        the host                     STILL PASSED to
                                                     every action

                      Host removed? No. Host as state domain? No.
                      Host as context? Yes. Seam A1 is not reopened.

exact caller:         src/class_graph_marcher.f90 (state_seat,
                      read_statement, march)

mathematical concept: delegation of a domain question to the object
                      that owns the mathematics

local necessity:      yes
global necessity:     yes for this citizen

cross-tower recurrence: the same correction shape as TI-9, one level
                      up

graph role:           unchanged

comparison:           TI-9 corrected the step; this corrects the
                      march. Both delegate rather than infer, and
                      neither removes anything

suspected nucleus implication: applied and complete for this
                      citizen.

confidence:           high
action:               ACT - applied, narrow, regression-verified
```

---

## OBSERVATION TI-16

```text
tower:                Time Integration
level:                8  (constitution)
review gate:          C
contextual radius:    2 (a production composition four deep)

symptom / fact:       *** THE CLASS-2 WITNESS, AND IT IS NATURAL. ***

                      DIFFERENCE_LINEARIZATION BUILT ITS STATES ON
                      THE HOST'S VERTICES.

                      This is the seam the reverse architecture
                      review named as Class-2, and it was reached
                      WITHOUT this tower ever naming the module.
                      The path is the one an implicit march
                      requires:

                        marcher % march
                          -> newton % solve
                            -> minimizer % attach
                              -> raw_apply
                                -> derivative_apply

                      Run after the TI-15 correction, VERBATIM:

                        ERROR STOP field_stored: a value vector
                        must fill its domain exactly

                        Error termination. Backtrace:
                        #3 __class_graph_field_MOD_field_set_real_vector
                           at src/field_stored.f90:305
                        #4 __class_graph_linearization_MOD_derivative_apply
                           at src/class_graph_linearization.f90:163
                        #5 raw_apply
                           at src/graph_minimization.f90:202
                        #6 __graph_minimization_MOD_attach
                           at src/graph_minimization.f90:173
                        #7 __class_graph_newton_MOD_solve
                           at src/class_graph_newton.f90:100
                        #8 __class_graph_marcher_MOD_march
                           at src/class_graph_marcher.f90:172

                      (PRE-FIX line numbers.)

                      derivative_domain had ALWAYS delegated
                      correctly. derivative_apply had not:

                        nv    = input_graph % num_vertices()   = 5
                        ncomp = max(size(at)/max(nv,1), 1)
                              = max(2/5, 1) = 1
                        state = field(..., input_graph % vertex_set(),
                                      ncomp=1)

                      NOTE THE max(...,1) FLOOR. It converts a
                      divide-to-zero into a silently WRONG width of
                      one, so the failure surfaces at the field
                      constructor as a shape error rather than at
                      the division as a domain error. A module that
                      says the right thing about its domain and the
                      wrong thing about where it builds its states
                      is exactly the shape this seam takes.

                      WHY THIS IS NOT MANUFACTURED: the import gate
                      refuses class_graph_linearization at EVERY
                      level of this tower, and --selftest asserts
                      that refusal. No level could have named the
                      module. It was reached by the production call
                      chain and by nothing else.

exact caller:         level-8-constitution/test.f90
                      (check_backward_march, through the governor);
                      src/class_graph_linearization.f90
                      (derivative_apply:137/156/162/170 as it stood)

mathematical concept: a finite difference of an operation is a
                      statement about that operation's domain

local necessity:      yes
global necessity:     YES for this citizen

cross-tower recurrence: SAME TOWER, third site. NOT a third vote.
                      What it adds is KIND, not count: seam A2 now
                      has a Class-2 production witness reached
                      through a four-deep composition, which no
                      derivative-family tower had produced

graph role:           conduit throughout the composition

comparison:           the reverse review listed
                      difference_linearization among A2's affected
                      files and marked global necessity "unknown",
                      with "test-local operations carrying explicit
                      domains" as the workaround. This is the first
                      time production ITSELF walked into it

suspected nucleus implication: APPLIED, narrow - see TI-17.

confidence:           high
action:               ACT (narrow), then observe
```

---

## OBSERVATION TI-17

```text
tower:                Time Integration
level:                8  (constitution)
review gate:          C
contextual radius:    2

symptom / fact:       A SAME-DOMAIN FINITE DIFFERENCE CAN PRESERVE
                      THE OPERATION'S DOMAIN WITHOUT BECOMING
                      RECTANGULAR - the correction TI-16 earned, and
                      the limit of it.

                      In src/class_graph_linearization.f90 ONLY:

                        derivative_apply asks
                          this % of % domain(input_graph, on)
                        and builds the frozen state, the perturbed
                        state and the J v answer on `on`

                        ncomp comes from size(at) / on % size(),
                        with exact divisibility required

                        a direction field, when present, must live
                        on `on` and match the frozen state's width

                        the underlying operation must ANSWER on `on`
                        - checked at both apply sites

                        the input_data-absent path was corrected
                        too, not only the branch GMRES exercises

                      *** THE LIMIT, AND IT IS THE POINT ***

                      THIS CITIZEN REMAINS SAME-DOMAIN:

                        L : Q -> Q

                      The tower does NOT earn, and this change does
                      NOT implement:

                        U -> Y with U /= Y      rectangular
                        transpose application
                        reverse action / adjoint
                        bidirectional linearization

                      SEAM B THEREFORE GETS ZERO NEW VOTES. The
                      forward use of difference_linearization inside
                      Newton is not bidirectional-linearization
                      evidence, and calling it that would be
                      counting a tower that did not run the
                      experiment. Seam B stands where it stood: two
                      independent derivative-family towers.

                      REGRESSION: difference_linearization has
                      exactly ONE consumer in the repository,
                      class_graph_newton, whose consumers are
                      graph-marching and graph-minimization. Both
                      were run IMMEDIATELY after this change and
                      passed.

exact caller:         src/class_graph_linearization.f90
                      (derivative_apply, answered_on)

mathematical concept: a linearization inherits the domain of what it
                      linearizes; same-domain is a narrower claim
                      than bidirectional

local necessity:      yes
global necessity:     yes for the SAME-DOMAIN citizen. UNKNOWN for a
                      rectangular one, which does not exist

cross-tower recurrence: third correction of one shape in one tower

graph role:           unchanged; the host is still passed through

comparison:           the reverse review bundles A2 and B as two
                      seams that "one experiment would settle."
                      This tower settles the A2 half for this
                      citizen and deliberately leaves B untouched -
                      which is only possible because the two are
                      genuinely different questions

suspected nucleus implication: applied and complete for the
                      same-domain case. The rectangular contract
                      remains unearned and undesigned.

confidence:           high
action:               ACT - applied, narrow, regression-verified
```

---

## OBSERVATION TI-18

```text
tower:                Time Integration
level:                8  (constitution)
review gate:          C
contextual radius:    1

symptom / fact:       A UNIFORM TIME GRAPH AND A SCALAR STEP ARE
                      SPECIALIZATIONS, NOT YET UNIVERSAL CONTRACTS.

                      Two facts about production, recorded WITHOUT
                      being called defects:

                        marcher regenerates a LINEAR CHAIN from
                        nsteps rather than consuming a supplied
                        G_time

                        marcher carries ONE SCALAR step rather than
                        the field h : E -> reals that Level 5 declared

                      For THIS specimen both are EXACT:

                        the time graph IS a simple chain - proved
                        extensionally in TI-12

                        h IS uniform - and Level 8 checks
                        h(e) = clock % step at every step, so the
                        scalar is an exact specialization of the
                        field rather than a coincidence

                      NO DEFECT IS ESTABLISHED. A specialization
                      that is exact on its specimen is not evidence
                      about the general contract, in either
                      direction.

                      These are the natural experiments for clients
                      that would supply what this one did not:

                        ADAPTIVE TIME     a genuinely varying h
                        COMPOSITE TIME    a branching or nonlinear
                                          time structure

                      One tower cannot decide whether per-edge h or
                      a supplied relational time graph becomes
                      load-bearing, and this one deliberately does
                      not try.

exact caller:         level-8-constitution/test.f90
                      (check_scalar_step_specializes_the_field,
                      check_control_chain_realizes_a1)

mathematical concept: an exact specialization of a general contract
                      on a particular instance

local necessity:      no - nothing here needs the general form
global necessity:     UNKNOWN, and deliberately left so

cross-tower recurrence: first tower to declare a richer time
                      structure than production consumes, and to
                      leave the gap open on purpose

graph role:           the marcher makes its own chain

comparison:           the partitioned tower left "what a genuinely
                      distributed road would need" DERIVED, not
                      implemented, under PIP-8. Same discipline: name
                      the frontier precisely, build nothing

suspected nucleus implication: NONE YET. Recorded so a future
                      Adaptive or Composite Time client knows exactly
                      which two assumptions to attack.

confidence:           high
action:               observe
```

---

## OBSERVATION TI-19

```text
tower:                Time Integration
level:                9  (statement)
review gate:          C
contextual radius:    1

symptom / fact:       THE COMPLETE STATEMENT IS ASKED AND ANSWERED
                      ON Q.

                      One initial-value problem:

                        Q = {x,y},  q(0) = [2,0],  qdot = -S(q),
                        S(q) = [x, y-x],  h = 1/2,  t0 = 0, t4 = 2,
                        bdf-2 with a backward-euler start

                        -> q(t4) = [7/24, 83/144]

                      Both ENDS of the statement are FIELDS ON Q.
                      The initial state arrives as one and its
                      vector is fetched once; the answer is written
                      back into one and its domain checked - the
                      field contract as the nucleus states it, fetch
                      once, work in arrays, write back once.

                      The marcher's raw-array core was left exactly
                      as it is. Nothing was refactored to make a
                      public argument prettier.

                      THE ENDPOINT IS EARNED, not assumed. Not
                      "member 5 is the last integer" but: time(t4)
                      = 2 read off the Level-5 coordinate field,
                      four steps walked, and the control chain's
                      terminal instant reached by FOLLOWING ITS
                      INCIDENCE from t0.

                      The marker carries the COMPUTED field, two
                      tokens, in Q's declaration order - not five,
                      because the answer lives on Q and a five-token
                      marker would undo nine levels of argument at
                      the last line of output.

exact caller:         level-9-statement/test.f90
                      (check_the_statement_s_two_ends,
                      check_the_endpoint_is_earned,
                      check_the_answer, say_the_result)

mathematical concept: a boundary stated in the terms the problem was
                      posed in

local necessity:      yes
global necessity:     yes - a statement that answered somewhere else
                      would not have answered the question asked

cross-tower recurrence: every sealed tower ends with a statement.
                      This is the first whose statement's domain is
                      not any graph's vertex set

graph role:           H_context is passed the whole way and is never
                      the seat of the answer

comparison:           the partitioned tower's L9 returned a field on
                      the global vertex carrier - which WAS a graph's
                      vertex set, correctly, because that problem
                      lived there. This one does not, and the
                      contrast is the tower's summary in one line

suspected nucleus implication: none.

confidence:           high
action:               observe
```

---

# Frontier

```text
THE TOWER IS SEALED. All ten levels are built.

WHAT GATE B ASKED, AND ANSWERED

    When q becomes a field on Q, can the existing temporal
    discretization and minimization stack preserve Q as the state
    domain, independently of whatever graph supplies structural
    context?

    YES - after one narrow correction, and the two halves of the
    stack answered differently:

        graph_minimization   already did          (TI-11)
        class_graph_step     did not; RED         (TI-8)
                             corrected            (TI-9)

WHAT GATE C ASKED, AND ANSWERED

    Can everything the tower earned be constituted into an actual
    multi-step march, and can a complete initial-value problem be
    asked and answered on Q?

    YES. Three narrow production corrections were required, each
    RED-first:

        class_graph_step          Gate B    (TI-8  -> TI-9)
        class_graph_marcher       Gate C    (TI-14 -> TI-15)
        class_graph_linearization Gate C    (TI-16 -> TI-17)

    q(t4) = [7/24, 83/144], a field on Q.

REVERSE EVIDENCE AT SEAL

    SEAM A1   CLOSED. Nothing here reopens it. The host is still
              passed to every action at every level, all the way
              through the march; this action ignores its topology
              only because a triangular 2x2 decay HAS none.

    SEAM A2   3 independent towers - Derivative Action, Adjoint,
              Time Integration - which reaches the strong-evidence
              threshold. Time is the FIRST non-derivative-family
              client, the FIRST full temporal production
              composition, and now carries a genuine CLASS-2
              witness reached naturally through
              marcher -> newton -> difference_linearization.

              Time contributes ONE vote. Three REDs, one client.

              RECOMMENDED NEXT ACTION: a dedicated reverse
              architecture review. NOT performed here, and
              doc/REVERSE-ARCHITECTURE-REVIEW.md is deliberately
              unedited - it is a separate artifact with its own
              process.

    SEAM B    STILL 2 independent derivative-family towers. ZERO
              new votes. difference_linearization remains a
              SAME-DOMAIN citizen L : Q -> Q, and its forward use
              inside Newton is not bidirectional-linearization
              evidence.

    A3        Another successful relational-graph ownership
              pattern - G_time owns two axes and four relations,
              including a carrier no relation names. No production
              change follows automatically. KEEP.

WHAT THIS LEDGER DOES NOT SAY

    that fit / reduction / broadcast must be migrated - none was
        exercised here
    that a rectangular or bidirectional linearization is earned
    that a nonuniform h or a supplied relational time graph is
        needed - see TI-18; both are exact specializations here
    that graph_operation, graph_grammar, graph_fitting,
        graph_calculus or graph_minimization must change - none
        of them did

A NOTE ON THE LABEL, because this ledger uses it both ways:

    SEAM A2   the reverse review's "operations should carry their
              own domain rather than ask a graph for one"
    A2        this tower's two-step reach relation, A1 o A1

They share two characters and nothing else. Every entry above means
the RELATION; the standing caution and this block mean the SEAM.

None of those is established by Levels 0-4, and the tower will be
worth less as evidence if it claims them early.
```
