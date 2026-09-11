# Topology ownership acceptance

`./run.sh` verifies the public read-only topology contract. The suite checks:

- Sparse relation fibres in both directions, duplicate tuple reduction,
  transpose identities and involution, empty fibres, and owning-copy independence.
- Interleaved readers with explicit per-call contexts, each retaining its own
  member count and sum without module-global callback state.
- Directed incidence and endpoint equivalence, boundary degree one, transpose
  orientation, and empty/default graph traversal.
- Compiler refusal of access to every stored graph component, fibre pointer
  storage, pointer-valued relation results, and mutable or `TARGET` readers.
- Constructor refusal of malformed extents, invalid endpoints, invalid partition
  indices or identities, noninjective partition maps, and invalid fibre indices.
- Zero dynamic allocations while reading actual fibre members, fibre arrays,
  explicit fibre contexts and the complete incoming index. An owning-copy measurement must
  allocate, demonstrating that the allocation counter is active.

The counter wraps `malloc`, `calloc` and `realloc` references in the executable
and static library. Allocations inside shared Fortran runtime libraries are
outside this measurement. Construction, owning copies, timing and output are
outside the measured traversal intervals. All traversal checksums are compared
with exact integer sums.

## Repeated performance comparison

Preserve the module files and static library compiled from commit `8d59154` in
a separate directory before building the candidate. Then run:

```sh
python3 test/graph-topology-ownership/benchmark.py \
    --baseline-library /absolute/path/to/baseline/lib \
    --candidate-library /absolute/path/to/candidate/lib
```

The script compiles one source against both interfaces. Baseline traversal uses
the original pointer fibres and directly stored incoming arrays; candidate
traversal uses scalar read methods and read-only array callbacks. It alternates
seven runs of each binary, verifies identical checksums and zero allocations,
and reports every sample, medians and ratios as JSON. `--cpu` optionally fixes
processor affinity. Default data: 20,000 members, degree six, 1,000 repetitions.

The array callbacks are module procedures so their timing measures the reader
interface without an executable-stack closure. Internal procedures that capture
local variables can require compiler-generated trampolines; their timing is
sensitive to stack layout. Application consumers therefore also need separate
measurements. The scalar and bulk measurements are reported separately because
they have different dispatch costs.

The `fibre_context` measurement uses the same explicit sum context in both
executables, passed to a module callback. This is the production traversal
pattern when reading a fibre must update operation-local state.
