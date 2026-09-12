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
`released_at` returns only the data whose final reader is that completed
step. Pairing with new data resets the position. Replacing or clearing one
rule changes neither the immutable dependency graph nor its cached order.

Primal rules borrow pointers into their execution only during `advance`.
The stored rule is cleared before returning, so no rule retains a pointer to
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
functional evaluation and is released at the end of that step. The
independent last-reader calculations formerly present in Taylor preparation
and chain differentiation are removed. Postprocessed forward derivatives use
the same release intervals for tangent towers.

Reverse derivatives retain primal state, tangent towers and lower-order
costates because higher derivatives and the final Lagrangian evaluation
still read them. Their storage grows with the horizon. The reported Taylor
storage pairs count maximum live and total allocated numerical entries in
the specified state/tower arrays; they do not measure total process memory,
graph metadata, or temporary solver storage.

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

Interleaving is supported; concurrent threads are not certified. Diagnostic
`tally` accounting is still application-wide, with every accounting scope
closed before `advance` returns. Bounded reverse-memory algorithms remain
separate architectural work.

## Concurrency prerequisites

Concurrent executions are not yet certified: no acceptance case runs two
executions on two threads, and the application-wide state listed below
is still shared. What exists is the build option and the three library
synchronisations that precede such a case.

`OPENMP=yes ./build.sh` and `OPENMP=yes ./application/build.sh` compile
with `-fopenmp` and link `libgomp`; the suite Makefiles read the same
variable from the environment, so a test build must use the setting of
the library it links. The default build omits the flag: every `!$omp`
directive and `!$` sentinel line is then a comment, and the serial
statements are the fallback. The serial build's demonstration output is
byte-identical to that of commit `aebe1f2`, and the OpenMP build under
`OMP_NUM_THREADS=1` and `OMP_NUM_THREADS=2` reproduces the same bytes
for every demonstration of the comparison set
(`artifacts/remaining-work-2026-09-11/r11/evidence-slices-1-4/`).

Three library states are synchronised. The version of every coarse
statement of a `multigrid` is that object's own count (`num_statements`),
compared only by its coarse minimizer, so two objects have independent
sequences. The identity serial of `next_token` is incremented and read
in one `!$omp atomic capture`, so two threads never receive one serial.
The counted-storage registry (`doc/topology-ownership.md`) performs
acquisition, binding and release inside one named critical region, the
owner's `clear` outside it. Owner reads take no lock.

Runtime and compiler requirements of the OpenMP build with GNU Fortran
15: `libgomp` at run time; `-fopenmp` implies `-frecursive` and
`-pthread`, so every local array of every procedure is on its thread's
stack and `OMP_STACKSIZE` must cover the largest automatic array of an
execution (the dense factorisation is allocatable, not automatic). The
library's debug flags include `-ffpe-trap=invalid,overflow,underflow`;
the trap mask is inherited by every thread, so a non-finite value in one
execution raises `SIGFPE` and terminates the process, whereas the solvers
already return `SOLVE_NONFINITE` through `ieee_is_finite` checks. A
concurrent build must either omit the trap for the library or accept
that a trapped execution is not isolated. `error stop` from any thread
terminates the process; three numerical failure paths of the application
still stop rather than return a failed result.

Still shared and therefore not thread-safe: the `util_tally` accounting
stack and amounts, output paths chosen by the caller, the six write
statements reachable from an execution, the gmsh loader's unit selection
by `inquire`, and the malloc counters of the benchmark instrumentation.
`verbosity` must be set before the first parallel region.

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
execution interleaving, and generic incremental execution suites. These test
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
