What is being built

Two files under maths/demo and a set of configuration files. physics_vanderpol.f90 supplies the residual integrand R and the functional integrand F with their state and design partials, for a Van der Pol system generalised to degree N as q^(N) − ν(1−q²)q^(N−1) + q = 0. graph_time_integrator.f90 assembles what src/ provides and prints one table, selected by --config=test_homogeneous and similar. Columns run f, df/dν, … to the configured derivative degree; rows are the homogeneous primitives, then ordered pairs, then ordered triples of BDF, ABM and DIRK, at one order or across orders as the configuration asks. No test suite is kept; a demonstration that runs is the check. Three further demonstrations follow: grid design, coefficient design, and the two together.

Structure

type(graph) is the only permitted structure, so a scheme is built from supplied vertices and edges. The hierarchy is self-similar at every zoom:

  tower      B_1 =====> B_2 =====> B_3        blocks, one family and order each
  block      G_0 --> G_1 --> ... --> G_m      slices, one instant each
  slice      G_k = ( G_k^S , G_k^C )          state first, constraint second
  sub-deck   G_k,1 --> ... --> G_k,s          one slice when not multistage
  degree     [ q , q' , ... , q^(N) ]

with G^S = (primal_state, (adjoint_state, tangent_state)) and G^C the same shape over the equations. That nesting records a real distinction: the primal traversal is nonlinear, while adjoint and tangent are both linear in the operator assembled at the converged primal state. At the top, P^S = (design, state deck) and P^C = (functional, constraint deck), which keeps the count square — N+2 data blocks against N+2 operators — with the design never solved and the functional's multiplier held at one, exactly as Figure 2 of the 2017 paper draws it. The pair order is (data, operator) to match view_epistemic.

A block is homogeneous in family and order; the tower above it is heterogeneous. Boundaries between blocks carry the state forward and the costate backward, and the tower's outer boundaries are the initial and terminal conditions, so those are the degenerate case of a junction rather than a separate mechanism. Only a multistage block may be first. With automatic_order_conservation = .true. a multistage startup block is prepended before a multistep block — 2P slices for BDF of order P, P−1 for ABM — so the order is preserved; .false. is unimplemented and stops. Every family carries a stage sub-deck, degenerate at one stage for BDF and ABM, so that all families are traversed by identical code; that identity is the evidence the abstraction is correct, and it also collapses the per-stage quadrature h_k Σ β_i F_ki to h_k F_k without a family test.

Indexing and the two linear systems
Fortran folds case, so distinct words are used rather than Q against q.
  max_state_degree           degree      0 .. max_state_degree   components, tangents                                                                                                                                                 condition   0 .. max_state_degree   constraints, costates  max_discretization_order   --                                  scheme order  max_derivative_degree      sensitivity 1 ..                    tangent order  max_stages                 stage       1 ..  max_outputs                output      1 ..                    functionals, costate  max_parameters             parameter   1 ..                    design components  max_instants               instant     0 ..  max_blocks                 block       1 ..
degree and condition share a bound because the block is square. The primal system has rows indexed by constraint and columns by component; the adjoint is that same matrix transposed, rows by component and columns by constraint. λ, ψ and φ are dual to R, S and T; q, q̇ and q̈ are dual to tons. The generic (N+1)m block is solved rather than its Schur complement, and thereduced Newton system of the 2017 paper is recovered exactly by eliminating the derived rows. Multiplication by ones and addition of zeros are accepted, the purpose being characterisation rather than speed. operation's existing max_degree is a different quantity from max_state_degree and must not be conflated.
                                                                                                                                                                                                  The scheme as a productEvery weight separates as α(θ)·h_k^(d′−d), the exponent being the derivative degree an it enters, which is determined by the topology alone. The dimensionless factor dependson the steps only through θ_j = (t_k − t_{k−j})/h_k, with α_j = ℓ_j′(0) for BDF and ∫_niform grid is the degeneracy θ_j = j, at which α collapses to the classical tables;this was checked against set_bdf, whose three BDF-2 lines are reproduced exactly and gtio. The product is edgewise over one shared topology, though the factors are notindependent since α is computed from τ. The coefficient builder takes the dt-graph andh and returns the extended one; entries are immutable at fixed design, which gives afree check — incremental and from-scratch construction must agree entry for entry.
Design as the parent notion
State and design are one kind of variable differing in whether the codomain map is square and closable, and in disposition, fixed or free. Solving is tuning, which the tower already encodes: every solver is a minimizer attached to a statement. Two refinements were taken from evidence in the repository: disposition is carried as a map rather than as branches, since it changes during a study and malready carries per-identity status; and state is composition rather than a subtype ofoural differs and the view_directed audit showed what a one-concretion hierarchy costs.The unification is of description and gradient assembly, not of traversal — the primaltuning outer.Carrying ψ and φ as genuine unknowns, rather than assuming them zero, is what makes coall, since ∂S/∂α and ∂T/∂α reach the gradient only through them.
                                                                                                                                                                                                  Consequences flushed out
                                                                                                                                                                                                  Six revisions are required before the design demonstrations will work. The coefficientable operation rather than a producer of data, since grid design needs ∂α/∂h. Thequadrature weight becomes design-dependent, so ∂f/∂h_k carries an explicit F_k term. The instants become design-dependent, so ∂R/∂t must be declarable and chained — invisible for Van der Pol, which has no explicit t. A tower-level constraint is needed, because a fixed duration and the order conditions have no slice to live at while the constraint deck is instant-indexed. The coefficient graph per-entry provenance, computed from τ or free, the two being exclusive. And immutabili.Two limits are scope rather than defect: the number of steps and the assignment of sliy sizes varying; and designing α is a choice between one tableau per block andindependent entries per slice. One caution: joint grid and coefficient design is reduneld by trading α against h, so the problem should be expected ill-conditioned until theduration or the order conditions are imposed.

The acceptance criterion adopted for the three design demonstrations is that they differ only by which entries a configuration file marks free. Any one of them needing its own code path means the abstraction has not earned itself.

Open, and the gaps in src/

Still to settle: whether the startup block takes its slices from the head of the firstnds the duration; the random step distribution, its bounds and its seed; and thedetailed contents of the configuration groups.

src/ supplies none of the traversal yet. There is no traversal over a descriptor step,rajectory-level derivative recursion; the coefficient tables are uniform-step only andstop at order four; and there are no Butcher tableaux or stage assembly. Those are what must be built before either demonstration file can be written.







```mermaid


graph TD
  P["P — problem"]
  P -->|B1| PS["P_S — data"]
  P -->|B2| PC["P_C — operator"]

  PS -->|B1| PSP["primal data"]
  PS -->|B2| PSD["derived data"]
  PSD -->|B1| PSA["adjoint data"]
  PSD -->|B2| PST["tangent data"]

  PSP -->|B1| XI["design — xi"]
  PSP -->|B2| QT["state tower"]
  PSA -->|B1| GR["gradient — df/dxi"]
  PSA -->|B2| LT["costate tower"]
  PST -->|B1| DIR["direction — p"]
  PST -->|B2| MT["tangent tower"]

  PC -->|B1| PCP["primal operator"]
  PC -->|B2| PCD["derived operator"]
  PCD -->|B1| PCA["adjoint operator"]
  PCD -->|B2| PCT["tangent operator"]

  PCP -->|B1| FF["functional — F"]
  PCP -->|B2| RT["constraint tower"]
  PCA -->|B1| SD["seed — e_j"]
  PCA -->|B2| AT["adjoint equation tower"]
  PCT -->|B1| SR["source"]
  PCT -->|B2| TT["tangent equation tower"]

```





```mermaid
graph TD
  X["node — tower, block, slice or stage"]
  X -->|B1| XS["data"]
  X -->|B2| XC["operator"]

  XS -->|B1| XSN["nonlinear"]
  XS -->|B2| XSL["linear"]
  XSN -->|B1| PAR["parameter — this level's own design"]
  XSN -->|B2| VAL["value — the two sub-nodes, or a field at a leaf"]
  XSL -->|B1| TAN["tangent — J forward"]
  XSL -->|B2| ADJ["adjoint — J transposed"]

  XC -->|B1| XCN["nonlinear"]
  XC -->|B2| XCL["linear"]
  XCN -->|B1| OBJ["objective"]
  XCN -->|B2| CON["constraint"]
  XCL -->|B1| TEQ["tangent equation"]
  XCL -->|B2| AEQ["adjoint equation"]
```





```mermaid
graph TD
  T["tower — [0, T]"]
  T -->|B1| A["earlier span"]
  T -->|B2| B["later span"]
  A -->|B1| A1["span"]
  A -->|B2| A2["block — family and order attach here"]
  B -->|B1| B1x["block"]
  B -->|B2| B2x["span"]
  A2 -->|B1| C1["earlier slices"]
  A2 -->|B2| C2["later slices"]
  C1 -->|B1| S1["slice — one instant"]
  C1 -->|B2| S2["slice"]
  S1 -->|B1| G1["earlier stages"]
  S1 -->|B2| G2["later stages"]
  G1 -->|B1| U1["stage"]
  U1 -->|B1| D1["lower degrees"]
  U1 -->|B2| D2["higher degrees"]
  D1 -->|B1| Q0["component — field on the state domain"]
```



instead of splitting `data` as `primal` and `derived`, should we split it as `linear` and `nonlinear`. The tangent and adjoint are two orientations of the same linearized operator. 









```mermaid
graph TD
  TW["tower"] -->|B1| BK["block"]
  TW -->|B2| TWR["rest of tower"]

  BK -->|B1| SL["slice"]
  BK -->|B2| BKR["rest of block"]

  SL -->|B1| ST["stage"]
  SL -->|B2| SLR["rest of slice"]

  ST -->|B1| DG["degree component"]
  ST -->|B2| STR["rest of stage"]

  DG -->|B1| VL["values — field on the state domain"]
  DG -->|B2| NIL["NULL"]
```







```mermaid
graph TD
  P["P — problem"]
  P -->|B1| PS["P_S — data"]
  P -->|B2| PC["P_C — operator"]

  PS -->|B1| PSN["nonlinear data"]
  PS -->|B2| PSL["linear data"]
  PSN -->|B1| XI["design — xi"]
  PSN -->|B2| QT["state tower"]
  PSL -->|B1| PSA["adjoint data"]
  PSL -->|B2| PST["tangent data"]
  PSA -->|B1| GR["gradient — df/dxi"]
  PSA -->|B2| LT["costate tower"]
  PST -->|B1| DIR["direction — p"]
  PST -->|B2| MT["tangent tower"]

  PC -->|B1| PCN["nonlinear operator"]
  PC -->|B2| PCL["linear operator"]
  PCN -->|B1| FF["functional — F"]
  PCN -->|B2| RT["constraint tower"]
  PCL -->|B1| PCA["adjoint operator"]
  PCL -->|B2| PCT["tangent operator"]
  PCA -->|B1| SD["seed — e_j"]
  PCA -->|B2| AT["adjoint equation tower"]
  PCT -->|B1| SR["source"]
  PCT -->|B2| TT["tangent equation tower"]
```







```mermaid
graph TD
  SEQ0["sequence - sensitivity 1..max_derivative_degree"] -->|element| G0

  G0["graph - tower"]
  G0 --> E0["epistemic - data | operator"]
  G0 --> R0["relational - blocks, junction coupling"]
  G0 --> S0["sequence - block 1..max_blocks"]
  S0 -->|element| G1

  G1["graph - block"]
  G1 --> E1["epistemic - data | operator"]
  G1 --> R1["relational - slices, scheme reach"]
  G1 --> S1["sequence - instant 0..max_instants"]
  S1 -->|element| G2

  G2["graph - slice"]
  G2 --> E2["epistemic - data | operator"]
  G2 --> R2["relational - stages, butcher coupling"]
  G2 --> S2["sequence - stage 1..max_stages"]
  S2 -->|element| G3

  G3["graph - stage"]
  G3 --> E3["epistemic - data | operator"]
  G3 --> R3["relational - components x constraints"]
  G3 --> S3["sequence - degree 0..max_state_degree"]
  S3 -->|element| G4

  G4["graph - derivative state"]
  G4 --> E4["epistemic - data | operator"]
  G4 --> R4["relational - spatial coupling"]
  G4 --> S4["sequence - freedom 1..max_freedoms"]
  S4 -->|element| G5

  G5["graph - degree of freedom, leaf"]
```





```mermaid
graph LR
  S1["branch - KNOWN"] --> C1["cell - graph"]
  C1 -->|"B1 KNOWN"| E1["element"]
  C1 -->|"B2"| S2["branch - KNOWN"]
  S2 --> C2["cell - graph"]
  C2 -->|"B1 KNOWN"| E2["element"]
  C2 -->|"B2"| S3["branch - KNOWN"]
  S3 --> C3["cell - graph"]
  C3 -->|"B1 KNOWN"| E3["element"]
  C3 -->|"B2"| N["branch - NULL, the end"]
```







```mermaid
graph TD

  subgraph L0["level 0 — expansion"]
    G0["graph — one identity"]
    G0 --> E0["epistemic"]
    E0 -->|B1| E0D["data — design xi, sweep values"]
    E0 -->|B2| E0O["operator — functional F, constraints"]
    G0 --> R0["relational"]
    R0 -->|B1| R0C["carriers — sweeps s = 0..r"]
    R0 -->|B2| R0R["relations — sweep s reads sweeps below s"]
    G0 --> S0["sequence — sensitivity 0..max_derivative_degree"]
  end
  S0 -->|element| G1

  subgraph L1["level 1 — horizon, the span 0..T"]
    G1["graph"]
    G1 --> E1["epistemic"]
    E1 -->|B1| E1D["data — state over the span"]
    E1 -->|B2| E1O["operator — constraints over the span"]
    G1 --> R1["relational"]
    R1 -->|B1| R1C["carriers — blocks"]
    R1 -->|B2| R1R["relations — junctions, weight 1"]
    G1 --> S1["sequence — block 1..max_blocks"]
  end
  S1 -->|element| G2

  subgraph L2["level 2 — block, one family and order"]
    G2["graph"]
    G2 --> E2["epistemic"]
    E2 -->|B1| E2D["data — slice states"]
    E2 -->|B2| E2O["operator — slice constraints"]
    G2 --> R2["relational"]
    R2 -->|B1| R2C["carriers — slices"]
    R2 -->|B2| R2R["relations — scheme reach, weight tau x alpha"]
    G2 --> S2["sequence — instant 0..max_instants"]
  end
  S2 -->|element| G3

  subgraph L3["level 3 — slice, one instant"]
    G3["graph"]
    G3 --> E3["epistemic"]
    E3 -->|B1| E3D["data — bundle Q_k"]
    E3 -->|B2| E3O["operator — constraints at k, plus recovery"]
    G3 --> R3["relational"]
    R3 -->|B1| R3C["carriers — stages"]
    R3 -->|B2| R3R["relations — butcher a_ij, j <= i"]
    G3 --> S3["sequence — stage 1..max_stages"]
  end
  S3 -->|element| G4

  subgraph L4["level 4 — stage"]
    G4["graph"]
    G4 --> E4["epistemic"]
    E4 -->|B1| E4D["data — components q, q', q''..."]
    E4 -->|B2| E4O["operator — one governing, N derived"]
    G4 --> R4["relational"]
    R4 -->|B1| R4C["carriers — components and constraints"]
    R4 -->|B2| R4R["relations — condition x degree sparsity"]
    G4 --> S4["sequence — degree 0..max_state_degree"]
  end
  S4 -->|element| G5

  subgraph L5["level 5 — component, one derivative order"]
    G5["graph"]
    G5 --> E5["epistemic"]
    E5 -->|B1| E5D["data — values"]
    E5 -->|B2| E5O["operator — spatial constraints, if any"]
    G5 --> R5["relational"]
    R5 -->|B1| R5C["carriers — freedoms"]
    R5 -->|B2| R5R["relations — spatial coupling, empty for an ODE"]
    G5 --> V5["set — freedom, extent in map_set, O(1) objects"]
  end

```





```mermaid
graph TD

  subgraph L2["level 2 — block, one family and order"]
    G2["graph"]
    G2 --> E2["epistemic"]
    E2 -->|B1| E2D["data — the member slices"]
    E2 -->|B2| E2O["operator — coupling and constraints"]
    E2D --> S2["sequence — instant 0..max_instants"]
    E2O --> R2["relational"]
    R2 -->|B1| R2C["carriers — slices, constraint instances"]
  end
  S2 -->|element| G3

  subgraph L3["level 3 — slice, one instant"]
    G3["graph"]
    G3 --> E3["epistemic"]
    E3 -->|B1| E3D["data — the member stages"]
    E3 -->|B2| E3O["operator — constraints at k, plus recovery"]
    E3D --> S3["sequence — stage 1..max_stages"]
    E3O --> R3["relational"]
    R3 -->|B1| R3C["carriers — stages"]
    R3 -->|B2| R3R["relations — butcher a_ij, j <= i"]
  end
  S3 -->|element| G4

  subgraph L4["level 4 — stage (optional for multistage)"]
    G4["graph"]
    G4 --> E4["epistemic"]
    E4 -->|B1| E4D["data — the member components"]
    E4 -->|B2| E4O["operator — one governing, N derived"]
    E4D --> S4["sequence — degree 0..max_state_degree"]
    E4O --> R4["relational"]
    R4 -->|B1| R4C["carriers — components and constraints"]
    R4 -->|B2| R4R["relations — condition x degree sparsity"]
  end
  S4 -->|element| G5

  subgraph L5["level 5 — component, leaf of the state hierarchy"]
    G5["graph"]
    G5 --> E5["epistemic"]
    E5 -->|B1| E5D["data — a field of values, no spine"]
    E5 -->|B2| E5O["operator — spatial constraints, NULL for an ODE"]
    E5D --> V5["set — freedom, extent in map_set, O(1) objects"]
    E5O --> R5["relational"]
    R5 -->|B1| R5C["carriers — freedoms"]
    R5 -->|B2| R5R["relations — spatial coupling, empty for an ODE"]
  end
```

