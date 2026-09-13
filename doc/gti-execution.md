# GTI execution ownership

An execution owns the configuration and numerical state for one discrete
trajectory. Two independently initialized executions may advance in any
interleaving without changing either result. Configuration is copied at
initialization; subsequent changes to the caller's configuration affect later
executions only.

`solver_context` owns linear solver choices, row selections, aggregate and
coarse selections, stopping rules, and the retained inner minimizer.
`march_context` adds nonlinear stopping rules, space/time coupling, and the
version sequence. Its `configuration()` method copies the configuration and
starts with no retained minimizer and no issued versions. Changes to settings
that affect a retained minimizer invalidate it; reuse also checks the system
size and component width.

`chain_execution` owns its `march_context`, expansion, trajectory, time grid,
block rules, incremental schedule, and optional streamed Taylor state.
`initialize`, `advance`, and `complete` define its primal execution. Each
`advance` evaluates one actual block, including a startup block when required.
`derivative` evaluates derivatives of a completed retained trajectory using
that execution's solver context. `take_results` transfers ownership of a
completed trajectory and expansion to the caller. It consumes those results.
Reinitialization discards the preceding execution, including partial work.
Taking results also clears remaining rule, schedule and solver storage.

The blocking `march_chain` interface initializes this same type, advances it
to completion, and transfers its results. Existing application, adaptive,
field, and demonstration consumers pass explicit contexts. Numerical helpers
that omit the optional context receive a local default for the duration of
that call. Initial field and adaptive construction receive the configured
stopping criteria before their first consistency solve. There is no mutable module default or shared Taylor execution.

## Dependencies and lifetimes

For every history instant of block b, its predecessor is the latest earlier
block covering that instant. `dependency_incidence` states this relation once
for primal, tangent, and costate scheduling. A startup block belongs to the
same relation. The transpose reverses these dependencies for executable
costate rules. The former combined graph with unevaluated reverse vertices
and the separate descending costate solve loop are removed.

The generic driver constructs its order and release intervals once. It owns
its completed position, and `advance` performs one rule, its writes, and its
final reads. `evaluate` uses that same kernel for a complete traversal.
`advance_with` performs the same step with a rule the caller passes (a
`vertex_rule`), which the driver tells its vertex and the vertices its
arguments read before applying it; the rule is not stored in the pairing.
`released_at` returns only the data whose final reader is that completed
step; `expired_at` adds the data the step's vertex writes that no rule
reads, the one lifetime statement the execution and the derivative passes
read. Pairing with new data starts the position at zero, or at the completed
position the caller names when the branch stores the data `live_after` that
position (written at or before it, read beyond it), which is the driver's
restart state. Replacing or clearing one
rule changes neither the immutable dependency graph nor its cached order.

Primal rules borrow pointers into their execution only during `advance_with`,
which returns without storing the rule, so no rule retains a pointer to
a caller's non-TARGET execution object. A derivative call similarly owns its
Taylor towers and retained costate derivatives. Forward rules reference
these arrays, and release follows the same dependency schedule without a
second numerical store. Each reverse rule computes its own right-hand side,
reads child costates from declared driver bindings, and emits its solved
costate. Accumulation uses descending child block order to preserve the
numerical order of reverse substitution. The driver releases these costate
fields after their final predecessor reads. A separate retained costate
array supplies the higher derivatives and final Lagrangian evaluation.

Streamed Taylor coefficients and each block's functional contribution are
computed immediately after the primal block. The execution uses the driver's
release intervals to deallocate both primal state and tangent towers.
An output with no subsequent reader is consumed by its own block's
functional evaluation and is released at the end of that step, as
`expired_at` states. The independent last-reader calculations formerly
present in Taylor preparation, chain differentiation and the execution's
advance, and the three placement loops (in-neighbourhood, argument count,
set, advance, clear) are removed. The startup family is registered once,
`gti_chain % startup_family`; execution ownership stays with the
application's `chain_execution`, which reads the library's driver for its
order, its step and its lifetimes and states no schedule of its own. Postprocessed forward derivatives use
the same release intervals for tangent towers.

Reverse derivatives are solved block-outer, order-inner in the transposed
order: at a block the costates of every order, multiset and functional are
solved in sequence, order k reading the block's own lower orders and the
child costates the bindings supply, and the block's Lagrangian terms are
evaluated at once from its state, tower and costates. The costates leave
the block as one datum on its unknown domain, released by the transposed
driver after the final predecessor read; the terms are stored per block
and summed ascending afterwards, so the tables equal the order-outer sums
bitwise. The Leibniz parts are per-block sums added ascending, no longer one
running sum over the rows of all blocks.

## Bounded reverse storage

The restart state of a completed position p of the forward traversal is
R_p = { (S_e, W_e) : e in `live_after`(p) }, the states and towers of the
blocks written at or before p and read beyond it; with the immutable rules,
expansion and configuration it recomputes every later block bitwise (Newton
starts from the fixed values, primal rows re-factorise at version zero, a
recomputed tower is solved under a version of its own so that each costate solve
reads the factors of its own tower only where the retained pass does, the
last block). `chain_derivative` and `derivative` read the solver context's
`reverse_entries` (application key `reverse_entries`, default `huge(1)`):
the largest sum of the live entries of five accounts - states, towers,
costates, restart states and per-block Lagrangian terms - returned as the
`derivative_storage` record with the schedule's quantities. With L the
largest restart state, F_max the largest block's forward entries, the window
the largest live costate data during a step with that step's own, and the
terms N nf nd M (2^(top+1) + top + 2) declared for the Lagrangian terms and
Leibniz parts whether or not the parts are requested, retention stores
sum F + window + terms; c stored restart states give
peak <= (c + 1) L + F_max + window + terms, and the limit admits
c = floor((limit - L - F_max - window - terms) / L), c = 0 recomputing every
block from the initial state; a limit at or above retention stores
everything, one below min(retention, the bound at c = 0) is refused with the
accounts before any tower is solved. The schedule is the recursion
`reversed(first, last, c)`: t(n, 0) = n (n + 1) / 2 forward evaluations,
t(n, c) = min over m of m + t(n - m, c - 1) + t(m - 1, c), a restart state
stored after block m carrying that block's own forward data. The bound
excludes the blocks' rows and rules, the expansion, the driver's copies of
the data it passes and the solvers' temporaries. Streamed Taylor storage
pairs count maximum live and total entries of the tower and state accounts.

A reverse derivative over a bounded primal is requested at `initialize`
(`functionals`, `derivative_order`, `pass_kind = reverse_pass`, optionally
`designs`, the leading designs of the tower), as the Taylor mode is. The
block sizes are read from the expansion before any block is solved
(`block_unknowns`), so the schedule is fixed at initialization: the march
solves each block's tower to top = order - 1 with the block, takes the
functional values, stores the restart state - states and towers of the
blocks live after the position - at the positions of the recursion's right
descent, retains the final block (every block under retention) and releases
everything else at `expired_at`. `derivative(reverse_pass, the same order
and functional count)` resumes the recursion with the forward driver paired
at the final position; a block recomputed is the execution's own block rule
applied again, the driver's data branch restored with the live states, so
the Newton solves of recomputation are counted by the tally beside the
tower solves; the costates read the versions of the first evaluation. A
later `derivative` finds no restart state stored and recomputes from the
initial state within the same bound. `take_results` of such an execution
returns no Taylor table. The streamed tables equal the retained post-hoc
pass bitwise (`test/gti-execution`; `test/graph-benchmark` `horizon`'s
`reverse` mode and `scaling.py`'s `reverse` series record the accounts, the
extra solves against retention, the peak resident size and the tables of
every restart count, equal across counts to 17 digits).

Periodic and event closure remain mathematical residual constraints. The
execution graph remains acyclic.

## Copies

An expansion is immutable once built: its hierarchy of level nodes and
its coupling binding are counted cells (`doc/topology-ownership.md`), and
its maps, layout, designs, rule and grid are values. A copy of an
expansion is one more owner of the two cells with its own values, so the
node identities are preserved: the copy is the same mathematical
expansion, and fields, domains and solver objects of the source are
defined on the copy's domains.

An execution is a value: trajectories (`chain(b) % state`), Taylor
coefficients, schedule position, results, block rules, the solver
context and every numerical cache are allocatable or scalar components
without pointers, and the transient rule pointers into the execution are
null between calls. Intrinsic assignment `copy = execution` therefore
produces an independent execution with equal state sharing only the
immutable expansion cells. Copies through a containing object, an array
element, a function result, a block local or an allocated allocatable
are the same copy. Numerical caches (the retained inner minimizer, a
dense factorization, a Gauss-Seidel block diagonal) are copied, not
invalidated: each is a function of the state it was computed from, the
copy has an equal state, and no hook exists that could invalidate them
under the copy mechanisms that bypass defined assignment. A copy and its
source reach exactly equal states, forward and reverse derivative tables
and streamed Taylor coefficients (`test/gti-execution`), and both equal
an isolated execution within the suite's declared tolerance.

The expansion is a nonallocatable component of the execution because
gfortran 15 does not invoke the cells' defined assignment through an
allocatable component (mechanism table in `doc/topology-ownership.md`).
`take_results` assigns the expansion to the caller's tower, which then
owns the cells, and leaves the execution without one.

Copies by `allocate(source=variable)`, a structure constructor or
polymorphic assignment are twins of one binding: their values are
independent and both compute the same results while both live, but
finalizing either, or taking either's results, releases the shared
binding and the other is refused at its next access with
`view_level: this storage's hierarchy has been released`. No storage is
read after release and nothing is freed twice. `allocate(source=)` of an
unevaluated or evaluated execution is therefore not a way to obtain a
second execution; assign instead. Whole-array assignment of arrays of
executions or expansions is a gfortran 15.2 internal compiler error
(`gfc_get_descriptor_field`); assign the elements.

## Limits

Interleaving is supported, and independent concurrent executions on the
threads of one process are certified under the conditions of
"Concurrent executions" below. Simultaneous blocks inside one dependent
trajectory, MPI or coarray images and device kernels are not. Reverse
storage is bounded by the schedule of "Bounded reverse storage" above.
Periodic or event closure remains separate architectural work.

## Accounting

Diagnostic accounting is a value, `tally` (`util_tally`), owned by each
execution's `march_context` as its `account`: the declared levels, the
amounts a(l, k, e) by level, derivative order and event, and the stack of
open levels. `configuration()` copies it restarted (same levels, orders
and recording state, zero amounts), so an execution begins with an empty
account of its caller's shape; `march_chain` adds the execution's
account into the caller's context after taking the results, and
`chain_execution % account()` returns it to a caller that drives the
execution itself. Every scope opened in `advance` is closed before
`advance` returns. A minimizer records through `record_event` into the
tally its caller bound with `bind_account` for the duration of one solve
(`solved` and `swept` bind the context's account before the solve and
null after), composites binding their children, so a stored inner
minimizer copied with an execution references no other execution's
tally; a factorisation is recorded by the minimizer that requests it.
No amount is stored in a module variable.

## Concurrent executions

Certified: independent concurrent CPU executions. Executions E_1, ...,
E_n with independent inputs (`chain_execution` values initialized from
their own `march_context`, copies of one shared template included) run
one per iteration of `!$omp parallel do schedule(dynamic,1)` and give,
per execution, the states, grid, derivative tables, Taylor coefficients,
storage pairs, failure results and counted accounting events of the
serial run at tolerance zero: every reduction inside one execution runs
in the program order of the serial build, and no cross-execution
reduction exists (a sum over executions is formed after the join in
index order). `test/gti-concurrent` checks this for ten heterogeneous
executions (state dimensions 3 and 4; BDF, Adams, DIRK and Newmark
families; direct, GMRES with Gauss-Seidel, GMRES with one multigrid
cycle and numerical elimination; designs and random grids; forward and
reverse derivatives of orders 1 and 2; one streamed Taylor execution;
one execution configured not to converge, whose derivative is a
non-converged result while the other nine equal the serial run and the
process exits normally), five repetitions at 1, 2 and 4 threads, each
execution writing its own file inside the region, every line of
standard output written after the join. The same source without
`-fopenmp` is the serial fallback: the requested thread count is
ignored and the suite passes with one thread.

Conditions: GNU Fortran 15.2 with `libgomp`; `OPENMP=yes ./build.sh`
(`-O3 -g -fbounds-check -fopenmp`, the default `FPETRAP=yes`) and
`-fopenmp` for every compilation unit that links the library (the suite
Makefiles read the same variable); `-fopenmp` implies `-frecursive` and
`-pthread`, so every local array is on its thread's stack and
`OMP_STACKSIZE` must cover the largest automatic array of an execution
(the dense factorisation is allocatable, not automatic); `verbosity`
set before the first parallel region. With the default trap a
non-finite value in any execution raises `SIGFPE` and terminates the
process, because the trap mask is inherited by every thread;
`FPETRAP=no ./build.sh` omits that one flag and the solvers return
`SOLVE_NONFINITE`, which is then a failure result of that execution
alone. `error stop` sites are invalid-input refusals and terminate the
process from any thread.

Shared state and its synchronisation. The version of every coarse
statement of a `multigrid` is that object's own count (`num_statements`),
compared only by its coarse minimizer. The identity serial of
`next_token` is incremented and read in one `!$omp atomic capture`. The
counted-storage registry (`doc/topology-ownership.md`) performs
acquisition, binding and release inside one named critical region, the
owner's `clear` outside it; owner reads take no lock. Unowned and
serialised only by libgfortran: the six write statements to standard
output reachable from an execution (`consistent_states`,
`initial_field`, the three `against_*` checks); start-up only: the gmsh
loader's unit selection by `inquire`; process totals with lost updates:
the malloc counters of the benchmark instrumentation. One location is
generated by the compiler: GNU Fortran stores the length of every
`character(len=:), allocatable` function result in a static variable
at the call site, with every flag combination tried
(`artifacts/remaining-work-2026-09-11/r11/evidence-slices-2-6/slice-6/deferred-length-defect`),
so two threads at one of the 73 such call sites of the library and the
application can exchange the lengths of two labels. No label result is
compared or selects a branch; labels reach diagnostics, refusal
messages and derived labels. The fix, results of
`character(len=this % name_length())`, changes the `name` binding of
every operation and is later work.

Instrumentation and measurement
(`artifacts/remaining-work-2026-09-11/r11/evidence-slices-2-6/slice-6`).
ThreadSanitizer (`-fsanitize=thread -fopenmp -O1 -g`, no trap) over the
suite at 4 threads: 17085 reports, every one classified: 16800
races on the static length above, with two reads past a label's own
buffer into a neighbouring freed block, the over-read an exchanged
length produces (`field_stored % create`); 223 on counted-storage cells
and registry counters accessed inside the critical region, whose libgomp
futex lock ThreadSanitizer does not observe (the ownership suite's
parallel cases pass with the region and fail without it), and
lock-order inversions between libgfortran's unit-table and unit locks;
nothing in the library's own state, and every numerical comparison of
that run exact. Valgrind memcheck of the serial build's suite: no
invalid access and no uninitialised value; the definitely-lost bytes
come from `stencil_term`, `consistent_states` (the R05 origins, solver
internals that predate this work) and `by_aspect` (the same
`class(field), allocatable` result pattern, unchanged since R06).
Throughput of the ten-execution set on four cores shared with two other
agents, medians of five repetitions: 56.9 executions per second on one thread (range 56.3 to 57.4), 92.2 on two (87.5 to 93.1) and 125.1 on four (109.2 to 133.5), speedups 1.62 and 2.20 of the median wall time, with peak resident memory (`VmHWM`) 9.3, 9.9 and 11.2 MB at the fifth repetition (`slice-6/measurement/measurement.json`). Peak resident memory
grows with the threads' simultaneous executions, not with the thread
count itself.

Not certified: simultaneous blocks inside one trajectory (an execution
is one thread for its whole life); MPI and coarray images
(`-fcoarray=single` here; the OpenCoarrays build shares nothing and is
outside this certification); device kernels; any run of the default
trapping build in which a non-finite value occurs; the demonstrations
of the OpenMP build at more than one thread beyond byte-identity of the
comparison set, which run one execution at a time.

## Failure results

A numerical failure inside an execution is a `solve_result` returned to
the execution's caller, never an `error stop`, so one failed execution
leaves the process and every other execution alive:

- a primal block that does not converge leaves the execution's
  `final_imbalance % converged` false, its `outcome` the solver's
  result, and the march continues on the state reached, as before;
- `derivative(..., outcome)` reports the primal result of the first
  block that did not converge (no derivative is solved, the table is
  zero) or the first linear solve of the pass that did not converge
  (`march_context % record_failure`, the first retained, reset at the
  start of every pass); `chain_derivative(..., outcome)` is the same
  contract for a caller that holds the chain;
- a streamed Taylor block whose primal did not converge is recorded as
  the failure of the march and its tower is solved at the state reached,
  so later blocks read defined values and `take_results` reports the
  non-converged `final_imbalance`;
- `adaptive_partition(..., outcome)` returns the step solve that did not
  converge with the steps accepted before it.

Without an `outcome` argument each of these entry points stops the
process with the same reason as before, which `test/gti-contract`
checks.

Output ownership: the paraview path of the main program is
`<export_path>_<label>_r<serial>_<instant>.vtu`, the serial being the
row's position in the run, so two rows of one label never write one
file. The six write statements to standard output reachable from an
execution (`consistent_states`, `initial_field`, the three `against_*`
checks) are serialised by libgfortran per statement and remain
unowned; a concurrent driver prints its own records after the join.

Still shared and therefore not thread-safe: the six write statements
reachable from an execution, the gmsh loader's unit selection by
`inquire` (start-up only), and the malloc counters of the benchmark
instrumentation. `verbosity` must be set before the first parallel
region.

## Solver restriction

A partitioned temporal solve constrains the residual to one member of the
unknowns and solves it with a copy of the inner solver template restricted
to that member. Restriction is each minimizer's own contract: `restrict`
receives the selection `selected(i)`, the whole index of the i-th selected
unknown, as an injective map of the selected domain into the whole. The
family is square, so the selected residual rows are the same indices.

The base discards every quantity evaluated on the whole (statement, affine
part, coupling, stored inputs, block diagonal, residual history), preserves
stopping rules, block width and component width, and refuses an empty,
repeated, out-of-range or block-splitting selection. Each solver with
metadata composes it with the selection and restricts its children by the
selection induced on their domains: an elimination restricts its inner by
the retained positions of the retained selected unknowns; a multigrid
restricts its aggregates, its smoother by the selection and its coarse
solver by the distinct aggregate labels in first-appearance order; Newton,
GMRES and a preconditioned solver pass the selection through; a dense
direct solver discards its retained factorization; the temporal minimizer
restricts its own partition labels, member order and seed transfers.

The operator on the selection is the residual constrained to it with the
exterior fixed; its affine part is the exterior contribution. The temporal
engine dispatches on no solver type, so a new composition, including one
defined only in a test, participates without engine edits.

## Residual boundary

Let Q be the state on the unknown domain U, nu the design on the point
domain P (one value per evaluation point), R(Q, nu) in Y the residual
and A = D_Q R(Q, nu): U -> Y its frozen linearization. The family is
square: row i is the equation of unknown i, so Y = U. `residual_operator`
owns U (`unknown_graph`, the host every consumer passes, and its vertex
set `unknown_domain`) and P (`design_domain`, the vertex set of
`point_domain`). Every consumer reads one frozen tuple (Q, nu) on U x P,
built by `frozen_tuple(x, nu)`, and the same fixed rows F with values h:

| consumer | definition | on F |
|---|---|---|
| value `apply` | stencils on Q plus the physics at the points | R_i = Q_i - h_i |
| state tangent `partial_action` | D_Q R[v], v on U | v_i |
| explicit `explicit_tangent` | J = D_Q R as triples | unit rows |
| design tangent `partial_action` | D_nu R[w], w on P | 0 |
| higher `partial_action` | D^m R[s_1..s_m], m <= the physics' degree, each s_k on U or P | 0 |
| frozen `linearize(rhs, transposed)` | w -> J w - rhs, or J^T w - rhs, on the same U and P, versioned | identity rows of J |
| constrained `constrain(free, h)` | R on the free unknowns with the exterior fixed, on a new U' and P' = `selected_points(free)` | the retained fixed rows |

A state, a direction in the state or a right side defined on a graph
other than U, a design or a direction in the design defined on a graph
other than P or with a value count other than the point count, and a
host graph of another identity are refused, whatever their length.
The residual types its supports from its own graphs: `state_fields`
(U, one value per unknown: state, direction, tangent), `residual_fields`
(Y = U: residual, costate, forcing) and `design_fields` (P, one value
per point). A `typed_field_domain` placed through a directed graph
reads its extent from that graph's vertex count, so no consumer states
an extent read from an array; the extent form of the constructor is
for an owner that stores the declared cardinality of the same identity
(a `discrete_domain`, the residual, a linearization's image). A value
vector that does not fill entries times components is refused at
placement. A `discrete_domain` types the placed law's inputs (state
with the law's component count, design with one value per point) and
a functional's state; the six unused factories for residual, costate,
forcing, direction, tangent and solution on the point set are deleted.
`domain` returns U with one entry per unknown, so a minimizer stated
on a residual reads Y = U from the residual. The linearization emits
its image with the entries and component count of the statement's
result and refuses a result whose entry count is not the domain's.
The temporal partition restricts a stored input of a residual, a field
on P with one value per point, to each member by `selected_points`
onto the constrained residual's P'.

The adjoint of J under the Euclidean pairing on U is the coordinate
transpose J^T, which `linearize(transposed=.true.)` states through the
stencil's orientation reversal. Under <u, v>_M = u^T M v the adjoint is
M^-1 J^T M; the coordinate transpose violates the adjoint identity by
u^T (J^T M - M J^T) v, and the M-weighted sensitivity of F = <g, Q>_M
is obtained from the transposed linearization by solving J^T lambda = M g
(`test/graph-minimization`, `check_residual_boundary`).

Consumers of the boundary (library and `application/module_graph_time_integrator.f90`):

| quantity | producer | consumers |
|---|---|---|
| residual value | `residual_apply` | `minimizer % evaluate` (Newton, temporal partition, `by_aspect`), difference-mode `linearization` |
| constraints | `fixed_rows`, `fixed` (`create`, `constrain`) | `partitioned_solve` (seeds x on F), `at_first_instant`, `sinks_of`, `forcing_of`, `costate_rows`, `lagrangian_term` |
| frozen inputs | `frozen_tuple`; `state_tuple` in `minimizer % evaluate` | `gti_march % frozen_inputs` (`solved`, `swept`, `frozen_at`, `dense_jacobian`, `gti-contract`), Newton (`evaluate(x, y, inputs)` then `freeze`) |
| explicit Jacobian | `residual_explicit_tangent` | Newton (explicit entries), `linearize`, `sinks_of`, `jacobian_of` (`by_aspect`, `dense_jacobian`) |
| tangent action | `residual_partial_action`, one variation | `linearization_apply` exact mode (Newton by the tangent action), `varied`, `design_partial` |
| adjoint action | `linearize(transposed=.true.)`, `stencil % reverse` | `solve_linear` -> `swept` (costate solves), `by_adjoint`, `derivative_rule_apply` |
| versions | `march_context % next_version`, `versioned` | `linearize`, Newton (stamps the explicit stencil), `partitioned_solve` (member versions), `elimination` (complement), `dense_direct` (retained factors by version and transpose) |
| higher partials | `residual_partial_action`, m >= 2 | `halley_correction` (`derivative_of`), `point_terms` (the Taylor towers read the physics expression directly) |

## Family extension

A time family is data of `operation_family`: a geometry tag and the
coefficients that geometry reads (the order and its functional for
Adams-Moulton and BDF, the tableau (a, b) for a diagonally implicit
Runge-Kutta method, the pair (beta, gamma) for Newmark). Every incidence
a discretization states is derived from that data by `row_pattern`,
`block_connectivity` and `stage_connectivity`, every coefficient by
`edge_coefficient`, every quadrature by `step_quadrature` and
`stage_weight`; the application embeds those edges into its tuple layout
and adds none of its own, and neither the temporal engine nor any
solver reads a family name.

Constructor-only, with no library and no engine edit: a tableau through
`dirk_family(a, b)`, a Newmark pair through `newmark_family(beta, gamma)`,
an Adams-Moulton or BDF order through `adams_family(p)` or
`bdf_family(p)`. The application registers the name in `family_named`
(name and order to constructor) and in `family_names`; nothing else
changes. Alexander's two-stage L-stable tableau, gamma = 1 - sqrt(2)/2,
is registered this way (`alexander`, order 2): its accuracy row runs in
the accuracy contract (T14) and `time-integration-tower` level 6 states
its connectivity and weights from the family, assembles its step map on
q' = lambda q from them alone, and checks the stability function, second
order under refinement, the tangent of the step in h from
`weights_terms` and the adjoint identity through the transposed solve.

Not constructor-only under the current design: a new geometry. The four
geometries are the branches of `select case (this % geometry)` in
`history_depth`, `primary_degree`, `row_pattern`, `step_quadrature`,
`stage_weight` and `edge_coefficient`; generalised-alpha, a Nystrom
method or a multistep with another row pattern is a fifth branch in each
of them and a constructor, not a registration by data. **This remains
open**: R12 proved extensibility with a physical law and a second
discretization use, not with a family geometry, so no geometry has been
added through these contracts and the six branches stand as the cost of
adding one. The startup of a multistep family is one
registered family (`gti_chain % startup_family`, Crouzeix's three-stage
tableau) and is not configurable by name.

## Law extension

A physical law is data of `gti_physics`: one Lagrangian over the state
fields and one multiplier each, built from the fields, the design,
constants, the four arithmetic operations, integer and real powers, and
sin, cos, exp, log and sqrt. The residual is its stationarity in the
first multiplier and each functional is the same Lagrangian at a zero
multiplier, so both read one tuple. The radial oscillator
`q'' + q - nu/q**3 = 0` was added that way, with `energy` and
`square_integral` beside it, and `git diff --stat -- src` is empty
across that slice: no library edit, no engine edit, no solver branch.

The zero state is **not** assumed to lie in a statement's domain.
`operation % defined_at_zero` declares it; an expression answers by a
structural scan of its own graph, a residual by the conjunction over its
rules. `minimizer % state` evaluates the constant part `A(0)` only where
it is defined, and `matvec` - with `imbalance`, `block_diagonal` and
`dense_matrix_of`, the operations of a linear solve - refuses a
statement that has none. Before that, every law had to have a value at
the zero state because the minimizer evaluated it there.

One limit the law extension measured and repaired, and one still open:

- The initial tuple supplies the components below the highest and the
  law closes the highest, `q''(0) = nu/q(0)**3 - q(0)` here, which
  depends on the design. `consistent_states` now returns that rate by
  implicit differentiation of the closure, `slope dQ_top/dnu =
  -dRule/dnu`, the same linear system its Newton step already forms;
  `residual_operator` carries it as `fixed_rate` beside `fixed`, since
  a fixed row reads `x(row) - h(row)` and its design partial is
  `-dh/dnu`, not zero. The three places that assembled a design
  derivative - the residual's own `partial_action`, the forward pass's
  `forcing_of` and the reverse pass's `lagrangian_term` - all read it.
  Before the repair the reported design derivative was **not** the
  derivative of the reported functional for a family whose rows read
  the acceleration at the initial instant: at 81 instants the central
  difference of the printed adams2 square integral over `nu` was
  1.1891925 against a printed 1.17836909, and the discrepancy halved
  with the step. It now reads 1.18919269571 against the same
  1.1891925, and the forward and reverse passes agree to 1.1e-15
  where they disagreed by 9.1e-3. Adams-Moulton 2, Adams-Moulton 3 and
  Newmark reach their orders on every quantity (E04, E05, E09);
  Runge-Kutta stages and BDF rows never read that acceleration and are
  unchanged. Van der Pol's rate is identically zero at `q'(0) = 0`, so
  every demonstration is byte-identical across the repair.
- **Open**: only the first design rate of a fixed value is carried. A
  closure nonlinear in the design would need `d^m h/dnu^m`; a repeated
  derivative in the physics design of a closed component is refused
  rather than reported as zero. Both laws in the repository close
  affinely in the design, so no run reaches that refusal.

- A law undefined at a point the solver visits stops the program inside
  a `pure` function several frames below the residual, naming neither
  the rule nor the point.

## Geometry extension

`spatial_geometry = circular` states the disc: a polar mesh whose
angular coordinate is identified across the seam, whose outer ring meets
the curved Neumann boundary `dq/dr = 0`, and whose centre is one
polygonal cell that is its own coarse cell. It is stated through the
same keys as the box and reaches `verify.sh` as the accuracy cases
D01-D08, with `git diff --stat -- src` empty across those slices. What
the disc measures, and does not:

- The balance summed over every cell is zero for every field, to
  1e-15 relative over 129, 513 and 2049 cells: an interior face is
  counted twice with opposite signs and a boundary face has the zero
  Neumann flux.
- Against `kappa` times the laplacian of `(r^2 - a^2)^2` the fitted
  balance at form degree 2 attains order 2 only on the centre cell
  (1.99). The interior converges at 1.74 and the boundary ring at 1.33,
  both declared limitations. At form degree 4 the polar fit is ill
  conditioned, the interior error reading 4.2e+9.
- The marched radial mode `J_0(z_1 r/a) cos(omega t)` nevertheless
  converges at order 2 (2.02) under the fitted balance, one order above
  the operator's own interior rate. Under the jet rows it converges at
  0.38: the compact form at degree 2 fits the second derivatives from
  axis-pure members over a neighbourhood the polar mesh does not align
  with the axes.

The older plan's Newmark, typed-field and continuous/discrete-domain
phases are closed with the source/consumer/test matrix in
`artifacts/remaining-work-2026-09-11/r07/matrix-final.md`: the family's
connectivities are the only source of incidence, the law governs the top
degree at every evaluation point including the arriving instant of a
staged step, fields read their extent from the graph that names their
support, and one continuous law placed on two point graphs keeps two
discrete-domain identities.

## Functional discretization error

`gti_chain % functional_error` estimates F(Q_exact) - F_h(Q_h) for every
functional of a marched chain from the retained coarse state alone: an
enriched chain of blocks (`enriched_family`: the family of order p + 1 on
the same instants, a staged family by its arriving-instant jets, the
identity prolongation of the instant jets; per configured block a block
of the coarse family over its first p + 1 instants, whose costate
transfers the junction sensitivity into the block before, then the
enriched block; the history of each read from the coarse chain by
`transferred`) is built by `block_from` at the coarse state, no primal
solve is performed, the enriched costates lambda+ are solved in descending
block order by `solve_linear` transposed at the frozen state with each
child's costate added on the transfer rows as the reverse pass does, and
the estimate is eta = - lambda+^T R+(P Q_h) + [F+(P Q_h) - F_h(Q_h)], the
residual part per step (the rows between two arriving instants) and the
quadrature part per step (the enriched family's complete rule against the
coarse rule). Class: asymptotically exact, I = eta / (F - F_h) = 1 +
O(h^(min(p+, 2p) - p)), order 1 in |I - 1| for p+ = p + 1; not a bound. The
`functional_error_estimate` carries the estimate, both parts, the scale
S = sum |w_k f(Q_k)|, F_h and F+(P Q_h), the largest fixed-row residual
of the enriched blocks (the transfer identity, zero under a consistent
prolongation) and the indicators by step. The coarse chain's costates are
not read: the enriched rows are differently scaled equations (a BDF row
against an Adams row), so no prolongation of lambda_h onto them exists;
the estimator's costates are its own, solved once per functional. The
application prints the estimate with `check = functional_error` and the
indicators with `indicators`; `test/accuracy-contract` measures the
effectivity and the localization (G01-G12, X13-X18).

On a spatial field the same enrichment in time estimates the temporal
error against the semi-discrete mode energy on the fixed mesh (G13);
enrichment in space (a refined mesh with a prolongation exact to O(H^6))
and a localized spatial error are not implemented and are declared
unsupported.

With `with_derivative` the same object carries the estimate of the error
of the design derivative G_h = dF_h/dnu, the order-1 Lagrangian:
eta_G = - lambda+'^T R+(P Q_h) + [F+_nu - lambda+^T R+_nu](P Q_h) - G_h,
the costate rate lambda+' from J+^T lambda+' = d/dnu[F+_Q - J+^T lambda+]
along the prolonged coarse tangent P w_h (`forward_block` on the coarse
chain, `costate_rows` at the multiset [1]) and the bracket from
`lagrangian_term` at s = [], j = 1 on the enriched chain. It is the
derivative of an asymptotically exact estimate, measured at order 1 in
|I - 1| with its estimate at the order of G_h - G (G14, G15); no per-step
indicator of a derivative functional is produced. The near-zero case is
the functional `mean` = int q dt over one period: the criterion divides
by S, never by |F_h| (G16), and an identically zero integrand gives
S = 0, lambda+ = 0 and eta = 0 exactly (G17).

`functional_error_partition` drives an adaptive grid by the estimate:
accept at |eta| <= tol S, otherwise divide every step with |eta_k| >
tol S / N into ceiling(h_k / h_k') equal steps, h_k' = h_k (tol S /
(N |eta_k|))^(1/(p+1)). It returns an `adaptation_outcome`
(`ADAPTATION_MET`, `ADAPTATION_UNMET` when the next grid would exceed the
instant limit, `ADAPTATION_NONFINITE`), never an accepted grid at
exhaustion; `grid_stationary_partition` takes the same limit and returns
the same outcome (or stops the program without an outcome argument). The
application reports either failure with nonzero status
(`test/gti-contract` modes `functional_error`, `adaptation_unmet`, cases
`dirk_functional_error`, `budget_unmet`; `test/accuracy-contract` A01,
A02, R14).

## Elimination storage

An `elimination` over retained unknowns K and eliminated unknowns E states
the whole system and passes its inner minimizer the Schur complement
S = J_KK - J_KE M with M = (I + N)^-1 D^-1 J_EK, D the eliminated diagonals
and N strictly triangular in the substitution order. Every coefficient of M
and of S is summed over its dependency paths before a dependant row reads
it or the inner minimizer multiplies it; the complement is formed row by
row with one accumulator over the retained columns, so no uncombined
product J_KE M exists. Storage is bounded by `max_entries`, in entries (one
coefficient with its index, 12 bytes), the sum of five accounts recorded in
`elimination % storage`: input (nnz(J) triples read from the stencil and
the retained block split out), substitution (J_KE, J_EK, N, the diagonals
and M), temporary (the position of every unknown and one row accumulator),
schur (nnz(S)) and factorisation (what the inner minimizer declares through
`minimizer % storage_entries` for the retained unknowns: 2 nk^2 for a dense
direct solve, (restart + 1) nk plus the preconditioner's for GMRES, the
block diagonal for Gauss-Seidel, the smoother's and the coarse minimizer's
for multigrid, the inner's for Newton and the temporal minimizer). Every
count is accumulated in 64-bit integers and a total beyond `huge(1)`, the
largest count an index array addresses, is refused whatever the limit.

Each allocation proportional to a count of coefficients follows a check of
the accounts against the limit. A refused statement stores no partition,
records `SOLVE_STORAGE_EXCEEDED` (a failed result) with the accounts at
the refusal, and its solve returns the initial residual with the unknowns
unchanged; Newton and the temporal minimizer propagate it as
`SOLVE_INNER_FAILED`. The limit is metadata: restriction keeps it, and the
application key `elimination_entries` reaches every restricted member
through the solver context's configuration copy. The pattern of the
complement stencil, the minimizer objects and the process's own memory are
outside these counts and are not claimed. No action-based Schur complement
is implemented: the complement is explicit, or the elimination is refused.

## Verification

The repository verification entry point includes context isolation, GTI
execution interleaving, independent concurrent executions
(`test/gti-concurrent`), and generic incremental execution suites. These test
different solver choices, problem dimensions, startup schemes, copied input
configuration, partial reinitialization, forward/reverse derivatives,
streamed Taylor storage, and execution copies before, during and after a
primal and a streamed Taylor march, alternating progress and
reinitialization, either destruction order, containers, array elements,
function results and the `source=` twin refusal. The generic tests also check nontrivial visiting
order, final-read release intervals, argument identity, and refusal cases.
Detailed results for this implementation are recorded in
`artifacts/gti-execution-2026-09-11/implementation.md`.

The final cold library and application builds passed. All 40 repository
suites and five verification demonstrations passed; the three affected GTI
suites and five demonstrations were rerun after the final changes. Five
additional comparisons against commit `13b19c3` produced byte-identical
output: Taylor state, Lagrangian expansion, sensitivity, randomized checks,
and the reduced order study. The forward/reverse dependency demonstration
also passed for all four examples.

The streamed DIRK3 check retained at most 120 tangent entries and 30 state
entries at 20, 40 and 80 instants. Its largest relative difference from the
retained-trajectory calculation was 1.68e-15. The execution tests additionally
verify that every streamed primal array is released after completion.
These are numerical-array counts, not total-memory bounds.
