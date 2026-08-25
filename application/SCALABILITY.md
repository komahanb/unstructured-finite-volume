# Where the framework stops scaling, on one processor

Measured 2026-08-25 on the graph-time-integrators branch. Every number
below came from a run; nothing here is inferred from the source alone.
`application/jacobian_shape.f90` reproduces the sparsity table,
`--config=uniform` with `--instants` and `--max-derivative-degree` the
timing ones.

The three axes are taken in the order they bite: state variables now,
design variables and functionals later but structurally.

---

## 1. State variables — measured, and the binding constraint

One `bdf2` row, state degree 2, uniform grid, no derivative columns.
Unknowns are `(instants - history) x components`.

```
 unknowns      51    111    231    471    951
 seconds     0.05   0.13   0.46   7.51  98.86
 per doubling   -   2.62   3.58  16.36  13.13
```

Cost grows as **n^3.7 to n^4.0** once the horizon is wide enough for
the elimination to dominate. The last two doublings are the asymptote;
the first two are graph traversal, which is linear and hides it.

### Why: a matrix that is 98 per cent empty is eliminated as if it were full

`jacobian_shape` measures the block jacobian's occupancy. `below` and
`above` are how far a nonzero reaches either side of the diagonal.

```
  scheme      unknowns   filled   below   above   % full
  bdf 1             63      196       8       2     4.94
  bdf 2             63      233      14       2     5.87
  bdf 3             63      258      20       2     6.50
  adams 2           63      223       3       1     5.62
  adams 3           63      253       5       1     6.37
  dirk 2           183      623       9       1     1.86
  bdf 2            183      753      14       2     2.25
  bdf 2            244     1179      27       3     1.98
```

The matrix is **banded**. `below` is the history the family reads,
which is `history_depth x components` and does not grow with the
horizon. `above` never exceeds three: it is the coupling within one
instant, where the row determining a degree reads the degree above it
at the same instant. Occupancy falls as the horizon widens - 5.9 per
cent at 63 unknowns, 2.0 per cent at 244 - because the band is fixed
and the square is not.

The band is bounded end to end. `handed_over` in `gti_chain` gathers
exactly `given` instants, which is the same history depth, so a
junction between two blocks of different layout does not reach outside
the band either.

Against that structure the solver does two things it need not:

- `operation_dense_direct` forms the matrix as `A(:,j) = matvec(e_j)`,
  one matvec per unknown, then runs Gaussian elimination with partial
  pivoting over the whole square. That is `O(n^3)` flops and `O(n^2)`
  storage on a matrix with `O(n b)` entries.
- Partial pivoting may exchange any two rows, so even a banded matrix
  handed to it would be treated as full.

At the largest measured horizon, `n = 951` and `b = 14 + 2`:

```
              dense          banded        ratio
 flops      2.87e+08       2.3e+05        ~1250x
 memory      7.2 MB        0.19 MB           38x
```

The ratio is `O(n^2 / b^2)` and therefore grows without limit. This is
the binding constraint, and it binds on memory before it binds on
time: `O(n^2)` storage per block is what makes a large horizon
impossible rather than merely slow.

### What one processor should do instead

A banded factorisation with partial pivoting inside the band (the
LAPACK `gbtrf` shape), `O(n b^2)` flops and `O(n b)` storage. The band
is known before the solve: `history_depth(equation_degree) x
components` below, and the largest degree gap above.

This is not a preconditioner or an approximation - it is the same
elimination restricted to the entries that can be nonzero, so it reaches
the same solution to the same tolerance.

### The iterative path does not settle this

```
 unknowns      51    111    231    471
 dense       0.05   0.13   0.46   7.51
 krylov      0.04   2.62   4.26   7.65
```

Krylov grows about linearly (1.6x, 1.8x per doubling) because it spends
a fixed iteration budget and each matvec is `O(n b)`. It crosses dense
at roughly 470 unknowns. But it does not converge on a difference
block - at 121 instants it returns a functional of 3.40 against the
dense 11.67 and is reported unconverged - because nothing
preconditions it and a row on the second derivative is weighted by
`1/dt^2`. Iterating buys bounded work, not a solution. A banded direct
solve buys both.

---

## 2. Design variables — latent, and foreclosed by the current shape

`design` is one real number today and there is one functional, so
nothing here is a regression. It is what the architecture forecloses.

`by_tangent` and `by_adjoint` each take one vector and return one
scalar, and each calls `dense_solve`, which builds the matrix and
factorises it from scratch. `chain_by_tangent` and `chain_by_adjoint`
do the same once per block. **No factor is kept anywhere in the
module.**

So `n_d` design variables cost `n_d` full factorisations per block,
where the right algorithm is one factorisation and `n_d` back
substitutions:

```
 now              n_d x blocks x O(n^3)
 available        blocks x O(n b^2)  +  n_d x blocks x O(n b)
```

At `n = 951`, `b = 16`, one block, that is a factor of about `1250 x
n_d` in flops - the banded gain and the reuse gain multiply, because
each is currently paid in full for every design variable.

### What one processor should do instead

A factorisation that outlives the solve which produced it: an object
holding the banded factors, constructed once per block per state, then
applied to as many right-hand sides as there are design variables. The
tangent then costs one factorisation plus `n_d` substitutions.

This is a design decision, not a tuning one. It needs a type that owns
factors, and `dense_solve`'s present signature - matrix in, solution
out, nothing retained - cannot express it.

---

## 3. Functionals — the same shape, transposed

`by_adjoint` has the identical structure: one `g`, one scalar back,
one fresh factorisation. `n_f` functionals cost `n_f` factorisations.

The correct choice between the two modes is by count, and the
framework cannot currently make it because neither mode amortises:

```
 n_d < n_f      tangent:  factorise once, n_d forward substitutions
 n_f < n_d      adjoint:  factorise once, n_f transposed substitutions
```

A banded factorisation serves both - the transpose of a banded matrix
is banded with the two bandwidths exchanged - so one object serves
both modes and the choice becomes a count comparison rather than two
separate code paths.

`functional_of` returns a scalar and `owned_gradient` fills one `g`
per block. Carrying `n_f` functionals means those become a second
dimension, not a loop over independent solves.

---

## 4. The gradient path carries an avoidable order

Measured as the difference between `--max-derivative-degree=1` and
`=0`, which is exactly the sensitivity work:

```
 unknowns      51    111    231    471
 gradient    0.00   0.01   0.08   0.41
 per doubling   -   3.35   6.68   5.22        ~ n^2.5
```

This is additive to the solve, and it is paid once per gradient. Its
cost has a cause that the state path does not share:

`chain_systems` calls `jacobian_of`, which forms the dense jacobian by
`n` partial-action calls. `dense_solve` then wraps that array with
`stencil(a, 'jacobian')`, and `create_dense` in `operation_stencil`
makes **one edge per entry, structural zeros included** - `n^2` edges,
three arrays of `n^2`, and a graph of `n^2` edges. `dense_direct` then
rebuilds the dense array from that graph by `n` matvecs.

So on the gradient path the matrix is materialised three times - dense
array, `n^2`-edge graph, dense array again - before an elimination
flop is done. At `n = 951` the middle representation alone is about
904,000 edges.

This does **not** happen on the Newton hot path: `gti_march` attaches
the `linearization` directly, so the matrix is formed once by matvecs
against it. The waste is confined to sensitivity, where it is paid per
gradient - which is per design variable once axis 2 opens up.

The fix follows axis 1: if the solver takes a banded matrix, the
gradient path hands it the band it already has, and neither the dense
array nor the `n^2`-edge graph is built.

---

## Ranked, by what binds first

1. **Dense storage and elimination of a banded matrix.** `O(n^2)`
   memory is the wall; `n^4` observed time is the symptom. ~1250x
   flops and 38x memory available at the largest size measured, and
   the ratio grows as `n^2/b^2`.
2. **No factorisation reuse.** Multiplies axis 1 by `n_d` and `n_f`.
   Latent today at one design and one functional, and it is the axis
   the user named first, so it should be settled before it is paid.
3. **Triple materialisation on the gradient path.** `n^2` edges to
   carry a matrix already held as an array. Confined to sensitivity,
   not the hot path.
4. **Krylov as the escape hatch.** Linear in `n` but does not converge
   on difference blocks. Superseded by 1, and should be retired as the
   remedy for large horizons once 1 is done.

Nothing here is about parallelism. Every item is a single-processor
algorithm or storage choice, and each one makes the eventual parallel
decomposition easier rather than harder: a banded solve over a chain
of blocks is what a time-parallel method partitions.
