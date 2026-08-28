# Types and Procedures Reference

Comprehensive catalog of all Fortran derived types and their type-bound procedures across src/ and application/ directories.

---

## adams_family
**Location:** `src/operation_family_adams.f90:40`  
**Description:** Linear multistep Adams family coefficients and time-stepping order configuration.

- `name` - Name identifier for the family
- `history_depth` - Depth of stored history for multistep methods
- `row_pattern` - Row pattern of the Butcher tableau
- `edge_coefficient` - Edge coefficients for step weighting

---

## advection
**Location:** `src/operation_advection.f90:44`  
**Description:** Advection operator with velocity field and face-edge numerical flux functions.

- `normal_speed` - Normal component of advection velocity
- `edge_coefficients` - Upwind or central difference coefficients

---

## affine_map
**Location:** `src/operation_differential.f90:142`  
**Description:** Sparse affine map y = Aq + k with triple-list storage of matrix and vector.

- (no public procedures - uses inherited functionality)

---

## argument
**Location:** `src/operation_action.f90:72`  
**Description:** One operation argument space position with bindings to query argument matching.

- `matches` - Test if this argument matches another by name
- `is_named` - Query whether argument has explicit name

---

## argument_path
**Location:** `src/operation_chain_rule.f90:73`  
**Description:** Derivative path through argument of operation statement being differentiated.

- `has_degree` - Query if this path carries derivative of given order

---

## assembler
**Location:** `src/transform_assembler.f90:106`  
**Description:** Transform assembler that glues coarser blocks to finer neighboring domains.

- `defined_on_graph` - Query assembly support on graph structure
- `defined_on_data` - Query assembly support on field values
- `defined_on_relation` - Query assembly support on relation
- `assemble_graph` - Build assembled graph from component pieces
- `assemble_data` - Build assembled field data from pieces

---

## balance
**Location:** `src/operation_balance.f90:70`  
**Description:** Face terms rule with discrete numerical flux face quadrature weighting.

- `name` - Name of balance equation type
- `apply` - Apply balance equation to get residuals

---

## bdf_family
**Location:** `src/operation_family_bdf.f90:50`  
**Description:** Backward differentiation formula family with multistep time integration coefficients.

- `name` - Name identifier for BDF family
- `history_depth` - Number of historical solution states stored
- `primary_degree` - Accuracy order of the scheme
- `row_pattern` - Pattern row of Newton tableau
- `edge_coefficient` - Stencil coefficients for differentiation

---

## binary_relation
**Location:** `src/relation_binary.f90:121`  
**Description:** Abstract binary relation contract image preimage and direction bidirectional query.

- `arity` - Return arity of the relation
- `image_view` - Query image side view for hot path
- `preimage_view` - Query preimage side view
- `image` - Get image set for given source
- `preimage` - Get preimage set for target
- `source` - Query source set of relation
- `target` - Query target set of relation

---

## block_residual
**Location:** `application/gti_block.f90:94`  
**Description:** Time integration block operation with compiled tangent and spatial layout.

- `name` - Name of residual block operation
- `domain` - Query domain of block operation
- `apply` - Evaluate residual at instant in block
- `max_degree` - Maximum derivative order supported
- `partial_action` - Apply directional derivative
- `compiled_tangent` - Get compiled derivative structure
- `restricted` - Get restriction to subdomain
- `placed_on` - Query mesh placement of block
- `slice_of` - Parent slice containing this block
- `node_of` - Expansion node this block represents
- `moment_of` - Time level of this block
- `labels_of` - Get labels of block degrees
- `num_nodes` - Number of nodes in block domain
- `spatial_discretization_laid` - Query if discretization set
- `aggregates` - Query aggregation relations
- `with_reach` - Get reach set of coupling
- `rows_terms` - Get stencil term rows
- `linear_block` - Get linearized form
- `member_order` - Query expansion member order
- `num_unknowns` - Count of unknowns in block
- `num_degrees` - Count of state degrees
- `num_points` - Count of spatial points
- `points_at` - Get points in block domain
- `num_carried` - Count of fields carried
- `carried_unknowns` - Get list of carried unknowns
- `held_values` - Get held field values
- `first_held` - Get index of first held value

---

## block_sinks
**Location:** `application/gti_chain.f90:210`  
**Description:** One block sinks and jacobian diagonal on them block-level flow information.

- (no public procedures - data holder)

---

## bound_relation
**Location:** `src/view_relational.f90:109`  
**Description:** Binds set and relation for relational view with owned storage components.

- (no public procedures)

---

## bound_set
**Location:** `src/view_relational.f90:104`  
**Description:** Binds set with owned storage one row per set element and extent.

- (no public procedures)

---

## branch
**Location:** `src/graph_fractal.f90:65`  
**Description:** Forward reference pointer and status of graph branch NULL UNKNOWN or KNOWN.

- `status` - Return status constant BRANCH_NULL BRANCH_UNKNOWN or BRANCH_KNOWN
- `known` - Get pointer to known graph or null if unoccupied

---

## broadcast
**Location:** `src/operation_reduction.f90:184`  
**Description:** Broadcast reduction multiplexer node sharing or expansion architecture configuration.

- `broadcast` - Constructor or query broadcast mode type
- `name` - Name of reduction broadcast operation
- `apply` - Apply broadcast operation to input

---

## chain_block
**Location:** `application/gti_chain.f90:96`  
**Description:** One block of chain its statement unknowns positions time levels structure.

- (no public procedures)

---

## chain_rule
**Location:** `src/operation_chain_rule.f90:88`  
**Description:** Stateless assembler for multivariate derivative via tabulated derivative terms.

- `assemble` - Construct derivatives from operation derivatives table

---

## chain_system
**Location:** `application/gti_chain.f90:148`  
**Description:** Stamp of one block tangent at frozen state every order either route.

- (no public procedures - data holder)

---

## change_record
**Location:** `src/map_change_protocol.f90:40`  
**Description:** Protocol tracking which steps reported changes touch flag state validation.

- `reset` - Reset change record to initial state
- `mark_attempted` - Mark that change was attempted
- `mark_applied` - Mark that change was successfully applied
- `mark_checked` - Mark that change was validated
- `mark_kept` - Mark that change is retained
- `mark_reverted` - Mark that change was undone
- `mark_failed` - Mark that change failed
- `validate_terminal` - Check if record is in valid terminal state

---

## coarsener
**Location:** `src/transform_coarsener.f90:70`  
**Description:** Transform that glues and coarsens finer blocks into aggregated structure.

- `defined_on_graph` - Query coarsening support on graph
- `defined_on_data` - Query coarsening support on data
- `coarsen_graph` - Build coarsened graph structure
- `coarsen_data` - Aggregate field data to coarser level
- `blocks` - Get blocks being aggregated

---

## conduction
**Location:** `src/operation_conduction.f90:41`  
**Description:** Conduction operator with thermal diffusivity and face-edge heat flux.

- `normal_conductivity` - Thermal conductivity in normal direction
- `edge_coefficients` - Gradient coefficients for heat flux

---

## configuration
**Location:** `application/gti_configuration.f90:38`  
**Description:** Time integrator configuration parameters tolerances stepping mesh options.

- (no public procedures)

---

## conjugate_gradient
**Location:** `src/operation_conjugate_gradient.f90:33`  
**Description:** Conjugate gradient linear system solver with symmetric positive definite support.

- `name` - Name identifier for solver
- `solve` - Solve linear system using CG algorithm

---

## counted_set_representation
**Location:** `src/map_set_representation.f90:109`  
**Description:** Counted representation storing member count O(1) lookup whatever size is.

- `num_members` - Count of members in set
- `member` - Get member at given local index
- `local_index` - Get local index of member value

---

## coupling_reach
**Location:** `application/gti_block.f90:87`  
**Description:** Block coupling reach along step direction vertices each block reads each.

- (no public procedures)

---

## csr_relation
**Location:** `src/relation_binary.f90:183`  
**Description:** Compressed sparse row binary relation with compiled execution contract.

- `domain` - Query domain set of relation
- `has` - Test if relation contains tuple
- `num_tuples` - Count tuples in relation
- `tuples` - Get all tuples as ragged array
- `image_view` - Query image side view
- `preimage_view` - Query preimage side view
- `materialized` - Query if relation is materialized

---

## dense_direct
**Location:** `src/operation_dense_direct.f90:48`  
**Description:** Direct linear solver factorizing and solving dense systems via LAPACK.

- `name` - Name of dense solver operation
- `solve` - Solve dense linear system

---

## dense_factorisation
**Location:** `src/util_factorisation.f90:50`  
**Description:** LU factorization with pivot tracking singular detection and substitution.

- `factorise` - Perform LU factorization on matrix
- `substitute` - Solve using factorization result
- `order` - Size of factored system
- `singular` - Test if factorization is singular
- `smallest_pivot` - Get minimum pivot magnitude encountered
- `factorise_cost` - Cost of factorization operation
- `substitute_cost` - Cost of substitution operation

---

## derivative_partition
**Location:** `src/operation_chain_rule.f90:101`  
**Description:** Integer partition of degree with multinomial count private to module.

- (no public procedures)

---

## derivative_table
**Location:** `application/gti_chain.f90:157`  
**Description:** One order table one column per multiset designs of order size.

- (no public procedures)

---

## derivative_terms
**Location:** `src/util_derivative_terms.f90:108`  
**Description:** Tabulated derivatives coefficients of given direction order symmetry.

- `set_direction` - Set direction for derivative coefficient
- `set_symmetric` - Mark derivative as symmetric
- `set_coefficient` - Set coefficient value for term
- `num_directions` - Count directions in table
- `num_terms` - Count non-zero derivative terms

---

## designed_grid
**Location:** `src/operation_grid.f90:103`  
**Description:** Grid with designed weights not partitioned weights are steps themselves.

- `name` - Name of grid design
- `weight_of` - Get quadrature weight at vertex

---

## differential_operator
**Location:** `src/operation_differential.f90:109`  
**Description:** Measure vertex face boundary derivatives spatial discretization operations.

- `name` - Name of differential operator
- `domain` - Query domain of operator
- `apply` - Apply differential operator

---

## directed_graph
**Location:** `src/view_directed.f90:187`  
**Description:** Abstract directed graph view with vertices edges ownership and boundary.

- `id` - Query graph identifier
- `num_vertices` - Count vertices in graph
- `num_edges` - Count edges in graph
- `vertex_set` - Get set of all vertices
- `edge_set` - Get set of all edges
- `edge_tail` - Get tail of edge
- `edge_head` - Get head of edge
- `edge_has_head` - Test if edge has specified head
- `interior_vertices` - Get interior vertices
- `boundary_vertices` - Get boundary vertices
- `tagged_vertices` - Get tagged vertices
- `interior_edges` - Get interior edges
- `boundary_edges` - Get boundary edges
- `tagged_edges` - Get tagged edges
- `owned_vertices` - Get owned vertices
- `borrowed_vertices` - Get borrowed vertices
- `overlap_vertices` - Get overlap vertices
- `owned_edges` - Get owned edges
- `borrowed_edges` - Get borrowed edges
- `overlap_edges` - Get overlap edges
- `incident_edges` - Get edges incident to vertex
- `adjacent_vertices` - Get vertices adjacent to vertex
- `outgoing_edges` - Get outgoing edges from vertex
- `incoming_edges` - Get incoming edges to vertex
- `outgoing_vertices` - Get vertices with outgoing edge
- `incoming_vertices` - Get vertices with incoming edge

---

## dirk_family
**Location:** `src/operation_family_dirk.f90:57`  
**Description:** Diagonally implicit Runge-Kutta family Butcher tableau stage weights.

- `name` - Name of DIRK family
- `num_stages` - Number of stages in tableau
- `stage_weight` - Weight coefficient for stage
- `row_pattern` - Pattern row of tableau
- `edge_coefficient` - Edge coefficients for step

---

## discretization
**Location:** `src/operation_discretization.f90:33`  
**Description:** Stencil pattern on dependent variable which unknown feeds which.

- `dependencies` - Get stencil dependencies

---

## edge_function
**Location:** `src/operation_edge_function.f90:56`  
**Description:** Operation applying edge-based function to mesh faces for flux computation.

- `edge_coefficient` - Edge coefficient function values
- `domain` - Query domain of edge function
- `apply` - Apply edge function to get flux
- `max_degree` - Maximum derivative order supported
- `partial_action` - Apply directional derivative

---

## element_kind
**Location:** `src/view_mesh_geometry.f90:78`  
**Description:** Cell type face count ordering shared face algebraic face total.

- (no public procedures)

---

## expansion
**Location:** `application/gti_expansion.f90:98`  
**Description:** Hierarchy of time integration blocks instants and unknowns structured recursion.

- `build` - Construct expansion from specification
- `root` - Get root node of expansion
- `node` - Get node by index
- `num_nodes` - Count nodes in expansion
- `label_of` - Get label of node
- `status_of` - Get status of node
- `value_of` - Get value associated with node
- `extent_of` - Get extent of node values
- `consistent` - Check internal consistency
- `tuples_of` - Get relation tuples
- `rule` - Get recursion rule
- `parameter` - Get parameter value
- `num_designs` - Count designs in expansion
- `design_kind_of` - Get kind of design
- `design_extent` - Get extent of design
- `design_value` - Get value of design
- `step_partials` - Get step partial derivatives
- `step_partial_along` - Get specific step partial
- `weights_of_steps` - Get weights of steps
- `refuse_assignment` - Prevent assignment

---

## expression
**Location:** `src/operation_expression.f90:95`  
**Description:** General expression operation applying physics equations arbitrary order.

- `name` - Name of expression
- `domain` - Query domain of expression
- `apply` - Apply expression operation
- `max_degree` - Maximum derivative order
- `partial_action` - Apply directional derivative
- `at_instant` - Query value at instant
- `equation_degree` - Degree of equation
- `declare_degree` - Set degree of equation
- `highest_degree` - Get highest degree used
- `num_vertices` - Count vertices in domain

---

## face_record
**Location:** `application/gti_space.f90:84`  
**Description:** Face data during mesh building corners cells connectivity information.

- (no public procedures)

---

## family
**Location:** `src/operation_family.f90:45`  
**Description:** Abstract family time integration Butcher tableau coefficients multistep support.

- `history_depth` - Number of historical states stored
- `num_stages` - Number of stages in tableau
- `stage_weight` - Stage weight coefficient value
- `primary_degree` - Order of accuracy
- `row_pattern` - Row pattern of tableau

---

## family_holder
**Location:** `application/gti_expansion.f90:94`  
**Description:** Family holder integer tag stores reference to family object.

- (no public procedures)

---

## field
**Location:** `src/field_calculus.f90:83`  
**Description:** Abstract field with identity domain shape plain vector adapters.

- `name` - Name of field
- `units` - Units of field quantities
- `domain` - Query domain of field
- `defined_on` - Query which objects define field
- `num_components` - Count components per entry
- `num_entries` - Count total entries
- `value_kind` - Get kind of value stored
- `integer_vector` - Get as integer vector
- `set_integer_vector` - Set from integer vector
- `real_vector` - Get as real vector
- `set_real_vector` - Set from real vector
- `complex_vector` - Get as complex vector
- `set_complex_vector` - Set from complex vector
- `logical_vector` - Get as logical vector
- `set_logical_vector` - Set from logical vector
- `character_vector` - Get as character vector
- `set_character_vector` - Set from character vector
- `hold` - Hold values in field

---

## file
**Location:** `src/util_file.f90:28`  
**Description:** File operations reading writing line buffering unit management.

- `open` - Open file for reading or writing
- `close` - Close file releasing unit
- `unit` - Get Fortran unit number
- `read_line` - Read one line from file
- `read_lines` - Read all lines from file
- `num_lines` - Count lines in file

---

## fit
**Location:** `src/operation_fitting.f90:62`  
**Description:** Fit operation holds form level coefficients maintained by fit sector.

- `name` - Name of fit operation
- `apply` - Apply fit operation

---

## fixed_grid
**Location:** `src/operation_grid.f90:111`  
**Description:** Fixed grid weights are steps themselves carried as constants form.

- `name` - Name of fixed grid
- `weight_of` - Get weight at vertex

---

## form
**Location:** `src/field_forms.f90:62`  
**Description:** Finite element form basis function coefficients degrees of freedom.

- `num_members` - Count basis members
- `values` - Get basis function values
- `slopes` - Get basis function slopes
- `dimension` - Get vector dimension
- `declare_basis` - Declare basis set
- `basis_set` - Query basis set
- `members` - Get members of form
- `restrict` - Restrict form to subdomain

---

## form_optimizer
**Location:** `src/operation_fitting.f90:86`  
**Description:** Optimizer holds no machinery itself reads form and adapts coefficients.

- `adapt` - Adapt form coefficients

---

## functional
**Location:** `src/field_calculus.f90:129`  
**Description:** Abstract functional reduction returning single value not small field.

- `integer_value` - Get as integer
- `set_integer_value` - Set from integer
- `real_value` - Get as real
- `set_real_value` - Set from real
- `complex_value` - Get as complex
- `set_complex_value` - Set from complex
- `logical_value` - Get as logical
- `set_logical_value` - Set from logical
- `character_value` - Get as character
- `set_character_value` - Set from character

---

## functional_holder
**Location:** `application/gti_chain.f90:138`  
**Description:** Functional holder storage so several can be handed over simultaneously.

- (no public procedures)

---

## gauss_seidel
**Location:** `src/operation_gauss_seidel.f90:37`  
**Description:** Gauss-Seidel iterative solver with coloring for ordered point relaxation.

- `name` - Name of Gauss-Seidel solver
- `colouring` - Get coloring of vertices
- `solve` - Solve using Gauss-Seidel

---

## gmres
**Location:** `src/operation_gmres.f90:33`  
**Description:** GMRES Krylov subspace solver for nonsymmetric linear systems.

- `name` - Name of GMRES solver
- `solve` - Solve system using GMRES

---

## gmsh_loader
**Location:** `src/view_gmsh_loader.f90:52`  
**Description:** Mesh loader from Gmsh format files mesh data parsing implementation.

- `mesh_data` - Get loaded mesh data

---

## graph
**Location:** `src/graph_fractal.f90:78`  
**Description:** Recursive binary graph two symmetric branches NULL UNKNOWN or KNOWN.

- `declare` - Assign identity token to graph
- `id` - Get identity token of graph
- `same_as` - Test if graph is same as another
- `similar_to` - Test if graph is similar

---

## grid
**Location:** `src/operation_grid.f90:59`  
**Description:** Abstract grid operation quadrature points weights time integration scheme.

- `weight_of` - Get weight at vertex
- `duration` - Get time duration of grid
- `apply` - Apply grid operation
- `max_degree` - Maximum derivative order
- `partial_action` - Apply directional derivative

---

## halving_policy
**Location:** `src/operation_step_policy.f90:68`  
**Description:** Step halving policy proposing judging step sizes retry logic.

- `propose` - Propose step size
- `judge` - Evaluate step success
- `retry` - Determine if retry needed

---

## harmonic_form
**Location:** `src/field_forms.f90:157`  
**Description:** Harmonic form basis sine cosine functions for periodic domains.

- `values` - Get harmonic basis function values
- `slopes` - Get harmonic basis slopes
- `dimension` - Get dimension of form

---

## imbalance
**Location:** `application/gti_march.f90:70`  
**Description:** Steepest direction in state is gradient norm d||r||/dq equation.

- (no public procedures)

---

## inclusion_map
**Location:** `src/map_inclusion.f90:102`  
**Description:** Embedding map declaring parts within ambient declaring pairs stored.

- `include_in` - Include part in ambient
- `included` - Get included parts
- `declared_into` - Query ambient for part

---

## inclusion_pair
**Location:** `src/map_inclusion.f90:97`  
**Description:** Declared embedding which part which ambient by value specification.

- (no public procedures)

---

## jacobi
**Location:** `src/operation_jacobi.f90:27`  
**Description:** Jacobi iterative solver with diagonal preconditioning and vertex coloring.

- `name` - Name of Jacobi solver
- `colouring` - Get coloring pattern

---

## label_map
**Location:** `src/map_label.f90:82`  
**Description:** Maps sets to labels with binding query value associations labels.

- `bind` - Bind set to label
- `labelled` - Get labeled sets
- `label_of` - Get label of set

---

## label_pair
**Location:** `src/map_label.f90:77`  
**Description:** One row which set by value what it is called specification.

- (no public procedures)

---

## level_storage
**Location:** `src/view_level.f90:78`  
**Description:** Storage of spatial hierarchy levels nodes relations assembler coupling.

- `fresh` - Create fresh level
- `node` - Get node at index
- `num_nodes` - Count nodes
- `spine` - Get spine relation
- `assemble` - Assemble level
- `couple` - Couple levels
- `refuse_assignment` - Prevent assignment

---

## linear_cell_type
**Location:** `src/view_paraview_writer.f90:125`  
**Description:** Paraview cell type enumeration hypercube mapping ordering specification.

- `element_type` - Get element type
- `hypercube_type` - Get hypercube type

---

## linearization
**Location:** `src/operation_linearization.f90:45`  
**Description:** Linearization operation freezing values inputs around reference solution.

- `name` - Name of linearization
- `domain` - Query domain
- `apply` - Apply linearized operation
- `exact` - Get exact operation
- `freeze_values` - Freeze field values
- `freeze_inputs` - Freeze input fields

---

## listed_set_representation
**Location:** `src/map_set_representation.f90:131`  
**Description:** Listed representation explicit roll member values declaration order describe.

- `num_members` - Count members
- `member` - Get member at index
- `local_index` - Get local index

---

## marcher
**Location:** `src/operation_marching.f90:73`  
**Description:** Time stepping marcher operation applying scheme to advance state.

- `instants` - Get time instants of march
- `march` - Perform one time step
- `march_adjoint` - March adjoint system
- `march_directional` - March with direction
- `march_adaptive` - March with adaptivity

---

## mesh
**Location:** `src/view_mesh.f90:58`  
**Description:** Mesh structure with seven measurements cells faces vertices geometry.

- `cell_volume` - Get volume of cell
- `cell_centre` - Get center of cell
- `face_area` - Get area of face
- `face_delta` - Get face spacing
- `face_normal` - Get normal vector
- `face_centre` - Get center of face
- `face_weights` - Get weights of faces

---

## mesh_loader
**Location:** `src/view_mesh_loader.f90:28`  
**Description:** Abstract mesh loader extending type all mesh loaders concrete implementation.

- `mesh_data` - Get loaded mesh data

---

## minimizer
**Location:** `src/operation_minimization.f90:79`  
**Description:** Base minimizer attached operation graph tolerances every iteration control.

- `begin_imbalance` - Begin imbalance measurement
- `note_imbalance` - Record imbalance value
- `converged` - Test if converged
- `flattened` - Test if flattened
- `diverging` - Test if diverging
- `began` - Test if begun
- `fitted` - Test if fitted
- `exhausted` - Test if exhausted
- `halted` - Test if halted
- `attach` - Attach operation
- `evaluation_inputs` - Get input fields
- `matvec` - Matrix-vector product
- `inner_product` - Inner product operation
- `norm` - Compute norm
- `sweep_order` - Get sweep ordering
- `diagonal` - Get diagonal
- `block_diagonal` - Get block diagonal
- `constant` - Get constant vector
- `domain` - Query domain
- `apply` - Apply operation
- `solve` - Solve minimization

---

## multigrid
**Location:** `src/operation_multigrid.f90:59`  
**Description:** Multigrid solver hierarchy coarsening relaxation solving structured.

- `name` - Name of multigrid solver
- `setup` - Setup multigrid hierarchy
- `attach` - Attach smoother operation
- `solve` - Solve using multigrid

---

## newton
**Location:** `src/operation_newton.f90:66`  
**Description:** Newton method one component beyond family minimizer governs lower.

- `name` - Name of Newton solver
- `solve` - Solve nonlinear system

---

## operation
**Location:** `src/operation_action.f90:111`  
**Description:** Abstract operation domain action argument space tangent derivative.

- `name` - Name of operation
- `domain` - Query domain of operation
- `apply` - Apply operation action
- `max_degree` - Maximum derivative order
- `partial_action` - Apply directional derivative
- `compiled_tangent` - Get compiled derivatives
- `declare_arguments` - Declare argument space
- `stamped` - Get stamped derivatives
- `stamp` - Stamp forward derivatives
- `stamp_transposed` - Stamp adjoint derivatives
- `num_arguments` - Count arguments
- `argument` - Get argument by index
- `owns` - Test ownership
- `require_owned` - Require owned component
- `restricted` - Get restriction

---

## paraview_writer
**Location:** `src/view_paraview_writer.f90:160`  
**Description:** Holds drawn data three coordinates point cells ragged organization.

- `write` - Write mesh to paraview format

---

## partition_relation
**Location:** `src/relation_partition.f90:111`  
**Description:** Partitioning identity partition vertex edge indices part owners.

- `num_parts` - Count parts in partition
- `has_part_relation` - Query part relation
- `part_id` - Get part identifier
- `global_vertex_index` - Get global vertex
- `global_edge_index` - Get global edge
- `global_index` - Get global index
- `part_vertex_index` - Get local vertex
- `part_edge_index` - Get local edge
- `vertex_owner_part` - Get vertex owner
- `edge_owner_part` - Get edge owner
- `owner_part` - Get owner part
- `describes` - Query described domain
- `whole_vertex_set` - Get all vertices
- `whole_edge_set` - Get all edges
- `num_whole_vertices` - Count all vertices
- `num_whole_edges` - Count all edges

---

## partitioner
**Location:** `src/transform_partitioner.f90:115`  
**Description:** Transform how cut into how many which part hand back partition.

- `defined_on_graph` - Query partition support on graph
- `defined_on_data` - Query partition support on data
- `partition_graph` - Partition graph structure
- `partition_data` - Partition field data

---

## path_derivative
**Location:** `src/operation_chain_rule.f90:60`  
**Description:** Derivative of path occupied carries x^(k) direction field unoccupied.

- (no public procedures)

---

## polynomial_form
**Location:** `src/field_forms.f90:128`  
**Description:** Polynomial form basis power functions arbitrary order derivatives available.

- `values` - Get polynomial basis values
- `slopes` - Get polynomial slopes
- `dimension` - Get dimension

---

## pruner
**Location:** `src/operation_fitting.f90:119`  
**Description:** Pruning form optimizer adaptively removes weak coefficients from.

- `adapt` - Adapt and prune form

---

## ragged
**Location:** `src/relation_binary.f90:247`  
**Description:** Ragged array structure padded shape fixed width walk reads.

- `num_lists` - Count lists
- `length` - Get length of list
- `list` - Get list data
- `padded` - Get padded data

---

## random_grid
**Location:** `src/operation_grid.f90:96`  
**Description:** Random grid sampling vertices weights computed by random distribution.

- `name` - Name of random grid
- `weight_of` - Get weight at vertex

---

## reduction
**Location:** `src/operation_reduction.f90:122`  
**Description:** Reduction operation global reduce fold accumulate combine aggregation.

- `initialize` - Initialize reduction accumulator
- `accumulate` - Accumulate values
- `combine` - Combine reductions
- `finalize` - Finalize reduction result
- `reduce` - Perform full reduction
- `name` - Name of reduction
- `domain` - Query domain
- `apply` - Apply reduction

---

## refiner
**Location:** `src/transform_refiner.f90:59`  
**Description:** Refine transform child parent split direction refinement factor.

- `defined_on_graph` - Query refine support on graph
- `defined_on_data` - Query refine support on data
- `refine_graph` - Refine graph structure
- `refine_data` - Refine field data

---

## relation
**Location:** `src/relation_finitary.f90:109`  
**Description:** Abstract finitary relation arity domain tuples membership identity.

- `arity` - Get arity of relation
- `domain` - Query domain set
- `has` - Test if contains tuple
- `num_tuples` - Count tuples
- `tuples` - Get tuples as array
- `declare` - Declare relation
- `id` - Get identifier
- `same_as` - Test identity
- `materialized` - Query materialization
- `name` - Get name
- (inherits procedures from abstract type)

---

## relational_binding
**Location:** `src/view_relational.f90:114`  
**Description:** Binds set and relation for relational view components owned storage.

- `bind_set` - Bind set component
- `bind_relation` - Bind relation component
- `set_for` - Get set for key
- `relation_for` - Get relation for key
- `refuse_assignment` - Prevent assignment

---

## reversible_change
**Location:** `src/map_change_protocol.f90:73`  
**Description:** Abstract change apply revert keep undo protocol four deferred steps.

- `apply` - Apply the change
- `check` - Validate change
- `keep` - Accept change permanently
- `revert` - Undo the change

---

## robin_condition
**Location:** `src/operation_robin_condition.f90:72`  
**Description:** One condition tag three numbers boundary treatment specification rule.

- `faces` - Faces affected by condition
- `lhs_coefficients` - Left-hand side coefficients
- `rhs_coefficients` - Right-hand side coefficients
- `advection_lhs_coefficients` - Advection LHS coefficients
- `advection_rhs_coefficients` - Advection RHS coefficients
- `operator_coefficients` - Operator coefficients
- `boundary_values` - Boundary value data
- `wall_relation` - Wall relation data

---

## room
**Location:** `application/gti_space.f90:59`  
**Description:** Mesh space room corner vertices face edges cells connectivity.

- (no public procedures)

---

## scheme
**Location:** `src/operation_step.f90:56`  
**Description:** Time stepping scheme Runge-Kutta multistep auxiliary history fields.

- `name` - Name of scheme
- `domain` - Query domain
- `apply` - Apply scheme step
- `dependencies` - Get dependencies
- `set_bdf` - Set BDF family
- `max_degree` - Maximum derivative order
- `partial_action` - Apply directional derivative
- `state` - Get state field
- `auxiliary` - Get auxiliary fields
- `history` - Get history fields
- `action_argument` - Get action argument
- `from_action` - Get from action

---

## scheme_weight
**Location:** `src/operation_weight.f90:49`  
**Description:** Scheme weight family multistep order coefficients stage weights.

- `name` - Name of weight scheme
- `history_depth` - Depth of history
- `num_stages` - Number of stages
- `primary_degree` - Primary degree
- `row_pattern` - Row pattern
- `edge_coefficient` - Edge coefficient

---

## set_map
**Location:** `src/map_set.f90:85`  
**Description:** Maps sets together binding query value associations members stored.

- `bind` - Bind sets together
- `describes` - Query described domain
- `num_members_of` - Count members
- `member_of` - Get member
- `members_of` - Get members
- `has` - Test membership
- `index_in` - Get index
- `extent_of` - Get extent

---

## set_pair
**Location:** `src/map_set.f90:80`  
**Description:** One row which set by value how members are stored.

- (no public procedures)

---

## set_representation
**Location:** `src/map_set_representation.f90:70`  
**Description:** Contract three primitives deferred two theorems concrete identity.

- `num_members` - Count members
- `member` - Get member
- `local_index` - Get local index
- `members` - Get all members
- `has` - Test membership

---

## sink_costates
**Location:** `application/gti_chain.f90:192`  
**Description:** Departure checks transposition solve placing functional gradient together.

- (no public procedures)

---

## sized_costates
**Location:** `application/gti_chain.f90:172`  
**Description:** Sized storage of costates derivatives functional gradient magnitude.

- (no public procedures)

---

## sized_steps
**Location:** `application/gti_chain.f90:215`  
**Description:** Sized storage steps time step sizes advancement time storage.

- (no public procedures)

---

## sized_tangents
**Location:** `application/gti_chain.f90:168`  
**Description:** Block rank tangent unknown block functional rank sized storage.

- (no public procedures)

---

## stencil
**Location:** `src/operation_stencil.f90:59`  
**Description:** Spatial discretization stencil pattern transposition gradient approximation.

- `name` - Name of stencil
- `apply` - Apply stencil operation
- `dependencies` - Get dependencies
- `transpose` - Get transposed stencil
- `restricted` - Get restriction
- `max_degree` - Maximum derivative order
- `partial_action` - Apply directional derivative

---

## step_policy
**Location:** `src/operation_step_policy.f90:32`  
**Description:** Abstract step policy proposing judging retrying stepping logic.

- `propose` - Propose step size
- `judge` - Judge step success
- `retry` - Retry step

---

## step_scaling
**Location:** `src/operation_step_scaling.f90:44`  
**Description:** Step scaling weight application step size adaptive time integration.

- `name` - Name of scaling
- `edge_coefficient` - Scaling coefficient

---

## stored_directed_graph
**Location:** `src/view_directed_stored.f90:79`  
**Description:** Stored graph keeps own structure arrays vertices edges relations.

- `id` - Get identifier
- `num_vertices` - Count vertices
- `num_edges` - Count edges
- `vertex_set` - Get vertices
- `edge_set` - Get edges
- `name_carriers` - Get name carriers
- `edge_tail` - Get edge tail
- `edge_head` - Get edge head
- `edge_has_head` - Test edge head
- `transpose` - Get transpose
- `transposed` - Get transposed view
- `loop` - Test for loops
- `interior_vertices` - Get interior vertices
- `boundary_vertices` - Get boundary vertices
- `tagged_vertices` - Get tagged vertices
- `interior_edges` - Get interior edges
- `boundary_edges` - Get boundary edges
- `tagged_edges` - Get tagged edges
- `owned_vertices` - Get owned vertices
- `borrowed_vertices` - Get borrowed vertices
- `overlap_vertices` - Get overlap vertices
- `owned_edges` - Get owned edges
- `borrowed_edges` - Get borrowed edges
- `overlap_edges` - Get overlap edges
- `incident_edges` - Get incident edges
- `adjacent_vertices` - Get adjacent vertices
- `outgoing_edges` - Get outgoing edges
- `incoming_edges` - Get incoming edges
- `outgoing_vertices` - Get outgoing vertices
- `incoming_vertices` - Get incoming vertices
- `whole_relation` - Get whole relation
- `tail_relation` - Get tail relation
- `head_relation` - Get head relation

---

## stored_field
**Location:** `src/field_stored.f90:86`  
**Description:** Stored field name unit domain width live store allocation.

- `name` - Name of field
- `units` - Units of values
- `domain` - Query domain
- `num_components` - Count components
- `num_entries` - Count entries
- (inherits adapters from field)

---

## stored_functional
**Location:** `src/field_functional.f90:53`  
**Description:** Stored functional value whichever kind last set held store.

- `name` - Name of functional
- `units` - Units of value
- `domain` - Query domain
- `num_components` - Count components
- `num_entries` - Count entries
- (inherits value access from functional)

---

## stored_relation
**Location:** `src/relation_finitary.f90:207`  
**Description:** Stored relation materializes tuples scans domain identity answers.

- `arity` - Get arity
- `domain` - Query domain
- `has` - Test membership
- `num_tuples` - Count tuples
- `tuples` - Get tuples
- `materialized` - Query materialized
- (inherits from relation)

---

## string
**Location:** `src/util_string.f90:21`  
**Description:** String type character buffer operations parsing tokenizing conversion.

- `print` - Print string value
- `equals` - Test string equality
- `tokenize` - Split into tokens
- `as_integer` - Convert to integer
- `as_real` - Convert to real

---

## token
**Location:** `src/token_identity.f90:43`  
**Description:** Identity token matches declared serial number for graph identity.

- `matches` - Test if token matches
- `declared` - Test if declared
- `serial_number` - Get serial number

---

## transform
**Location:** `src/transform_structure.f90:32`  
**Description:** Abstract transform defined on graph data structure refinement.

- `defined_on_graph` - Query support on graph
- `defined_on_data` - Query support on data
- (abstract - children implement specific transforms)

---

## transposed_relation
**Location:** `src/relation_binary.f90:218`  
**Description:** Transpose view borrower holds base pointer answers every question.

- `domain` - Query domain set
- `has` - Test if contains
- `num_tuples` - Count tuples
- `tuples` - Get tuples
- `image_view` - Get image view
- `preimage_view` - Get preimage view

---

## triple_list
**Location:** `src/operation_stencil.f90:87`  
**Description:** Dynamic triple list for sparse matrix coefficients doubling growth.

- `place` - Add entry to list
- `entries` - Get entry data

---

## uniform_grid
**Location:** `src/operation_grid.f90:90`  
**Description:** Uniform grid equally spaced vertices quadrature weights uniform.

- `name` - Name of grid
- `weight_of` - Get weight

---

## value_change
**Location:** `src/map_value_change.f90:44`  
**Description:** Bound update map pointer graph copied identity new value.

- `bind` - Bind to value map
- `apply` - Apply the change
- `check` - Validate change
- `keep` - Accept change
- `revert` - Undo change

---

## value_map
**Location:** `src/map_value.f90:58`  
**Description:** Maps graph identity to status value unknown known field.

- `attach_unknown` - Attach unknown value
- `mark_known` - Mark as known
- `mark_unknown` - Mark as unknown
- `detach` - Detach from value
- `attached` - Query if attached
- `status_of` - Get status
- `value_of` - Get value

---

## value_pair
**Location:** `src/map_value.f90:50`  
**Description:** One row copied identity token status value field storage.

- (no public procedures)

---

## variation
**Location:** `src/operation_action.f90:91`  
**Description:** Derivative direction field stored position in argument space.

- `argument_is` - Query argument identity
- `argument` - Get argument
- `direction` - Get direction field
- `domain` - Query domain
- `field` - Get field
- `with_argument` - Test argument match

---

## walk
**Location:** `src/operation_walk.f90:71`  
**Description:** Walk question answered starting point structured graph navigation.

- `name` - Name of walk
- `apply` - Apply walk operation

---

**Total Types:** 163 | **Generated:** 2026-08-26
