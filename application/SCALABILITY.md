# Where the framework stops scaling, on one processor

Measured 2026-08-25 on the graph-time-integrators branch. Every number
below came from a run; nothing here is inferred from the source alone.
`application/jacobian_shape.f90` reproduces the sparsity table,
`--config=uniform` with `--instants` and `--max-derivative-degree` the
timing ones.

The three axes are taken in the order they bind: state variables now,
design variables and functionals later but structurally.

---

## 0. A convergence tolerance that stops being reachable

This binds before any of the three axes, and it is not an asymptotic
cost at all - it is a threshold, and everything downstream of it was
misread as an asymptotic cost until it was separated out.

`gti_march % solved` sets `tolerance = 1e-12 x max(1, scale_of(...))`.
The iteration counts below are counted by the accounting layer, not
inferred from a ratio of times as they were when this was first
written.

It is one scheme, not every scheme, and it is the step and not the
count of unknowns. At a fixed duration of three, raising the instants:

```
 instants   bdf1 loops   bdf2 loops   bdf2 reported
      106         8            8          converged
      112         8            9          converged
      115         8           37          converged
      121         8           40          converged
      161         8           40          converged
```

bdf1 takes eight iterations at every size and never degrades. bdf2
reaches the iteration limit of forty and stays there, and the row prints as
converged while it does: newton stops at the iteration limit, the residual it
stops at is 1.4e-11, and show_row reports a march as unconverged only
above 1e-8.

Holding the instants at 161 and varying only the duration shows what
the count follows:

```
 dt          0.15    0.075   0.0375   0.01875
 bdf2 loops     3       14       10        40
```

**The step, not the count.** Halving dt from 0.0375 takes bdf2 from ten
iterations to the iteration limit, at 161 instants either way. Refining a
grid is the whole of what scaling up means, so this is met by every
run that refines rather than by every run that grows.

The reference is the wrong one. `scale_of` is the norm of the residual
at the state newton starts from, and that state contains one instant's
components repeated. A difference row of a constant vanishes, by the
same consistency that makes the coefficients sum to zero, so the
starting residual contains none of the `1/dt^2` weight the rows contain. The
rounding-error floor at the solution does contain it. One rises as the
grid refines and the other does not, and where they cross the march
cannot reach the requested tolerance.

A first-order scheme's rows are weighted by `1/dt` and a second-order
scheme's by `1/dt^2`, which is why bdf1 is untouched at every size
measured and bdf2 is not.

**Settled, and not by a tolerance.** The jacobian's diagonal was
proposed here as the reference and it is the wrong quantity. Measured
on a bdf2 block of 244 unknowns the largest entry anywhere is 1.71e5
and the largest on the diagonal is 9.52: the sign convention puts one
on the column a row determines and the `1/dt^d` weights on the
sources, which are off the diagonal. Scaling by the diagonal changes
the tolerance not at all, the starting residual being larger anyway.
The largest row sum does carry the weight - 5.12e5 on the same block -
but scaling by it gives 5.1e-7, looser than the 1e-8 a row is reported
unconverged above.

No tolerance is the control parameter, because the floor is not at a value
a tolerance can be set to. Tracing every iteration of bdf2 at 161 instants:
the residual falls quadratically to 2.7e-11 by the eighth, and then
oscillates between 2.3e-11 and 3.2e-11 for thirty-two more. The tolerance
was 1.25e-11, a factor of two under a floor the iteration cannot reach.

So newton stops when it stops progressing: a residual that has not
improved on the best recorded by a tenth, four iterations running, once it is
already a millionth of where it began. The last clause is what makes
the rule valid - early newton iterates are non-monotone, a residual rising once before it falls
is ordinary, and without the clause the rule stops marches in their first
few iterations and reports them unconverged.

```
 instants   bdf2 loops before   after   wall before   after
      115           37            12       3.2887    1.3472
      121           40            12       3.9669    1.4736
      161           40            12       7.7766    2.9327
      181           40            13      10.5322    4.1743
```

Every functional is identical to all twelve digits printed, the four
configurations are byte-identical, and the two randomized sweeps are
unchanged.

---

## 1. State variables — the cost that remains once the threshold is set aside

One `bdf2` row, state degree 2, uniform grid, no derivative columns.
Unknowns are `(instants - history) x components`.

```
 unknowns      51    111    231    471    951
 seconds     0.05   0.13   0.46   7.51  98.86
```

An earlier account of this row called the growth `n^3.7` to `n^4.0`
and took the last two doublings for an asymptote. That was wrong. The
step between 231 and 471 is the threshold of section 0, and averaging
across it yields a power that nothing in the algorithm produces.

Taken apart, each part matches theory:

```
  unknowns    form(s)   solve(s)      form growth   solve growth
        63      0.002      0.005
       123      0.006      0.034            3.0          6.8
       243      0.023      0.252            3.8          7.4
       483      0.086      1.901            3.7          7.5
       723      0.191      6.453            2.2          3.4
```

Formation is quadratic and low cost - `n` partial-action passes at `O(n)`
each, 0.19s of 6.6s at 723 unknowns, three per cent. **The solve is
cubic** at every step, which is the elimination and nothing
more. The cubic term is the one to reduce; the quartic never
existed.

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

The band is bounded end to end. `transferred` in `gti_chain` gathers
exactly `given` instants, which is the same history depth, so a
junction between two blocks of different layout does not reach outside
the band either.

Against that structure the solver performs two operations it need not:

- `operation_dense_direct` forms the matrix as `A(:,j) = matvec(e_j)`,
  one matvec per unknown, then runs Gaussian elimination with partial
  pivoting over the whole square. That is `O(n^3)` flops and `O(n^2)`
  storage on a matrix with `O(n b)` entries.
- Partial pivoting may exchange any two rows, so even a banded matrix
  passed to it would be treated as full.

At the largest measured horizon, `n = 951` and `b = 14 + 2`:

```
              dense          banded        ratio
 flops      2.87e+08       2.3e+05        ~1250x
 memory      7.2 MB        0.19 MB           38x
```

The ratio is `O(n^2 / b^2)` and therefore grows without limit.

On memory the arithmetic above is right about the order and wrong
about what dominates. Measured, the matrix is not the dominant
memory:

```
 unknowns   peak, factorising   peak, matrix-free
      471          84.7 MB            82.1 MB
      711           185 MB             119 MB
```

The krylov path forms no matrix at all and at 471 unknowns uses
nearly the same memory. So the memory was ascribed to the
representation. It measures otherwise. One part per run, because
within a process the allocator's arena is already grown and a difference
reports nothing:

```
 at 711 unknowns and 9165 edges          peak MB
  a run that builds nothing                 1.54
  a bare vertex set                         1.50
  the same set given its banded edges       3.19
  a field over it                           3.15
  the block statement the march solves      4.07
```

**The representation is low cost.** A vertex set costs nothing
measurable, a field over it nothing, an edge about 184 bytes, and the
whole block statement 2.5 MB above a bare run - one and a half per
cent of the 185 MB that same size reaches while solving.

The memory is transient, and it is the threshold of section 0 measured
in memory. `dense_direct_solve` allocates its `n by n` array inside
the solve, so a new array is allocated every newton iteration and
deallocated; at 471 unknowns that is 1.7 MB an iteration, and the allocator
does not return it. Peak memory over the sweep is 7.5, 8.1, 10.3, 14.6
then 84.7 MB - flat while the march converges in a few iterations, six
times higher as soon as it uses all forty. Per iteration it
is about 2 MB either side of the step, so the step is the iteration
count and nothing else.

Two consequences follow, in order. Settle the tolerance and the memory falls
with the time, by the same factor. Then move the matrix allocation out of the
iteration: it is the same shape every time, and there is no reason to
allocate a new one per step.

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

Krylov grows about linearly (1.6x, 1.8x per doubling) because it runs
to a fixed iteration limit and each matvec is `O(n b)`. It crosses dense
at roughly 470 unknowns. But it does not converge on a difference
block - at 121 instants it returns a functional of 3.40 against the
dense 11.67 and is reported unconverged - because nothing
preconditions it and a row on the second derivative is weighted by
`1/dt^2`. Iterating gives bounded work, not a solution. A banded direct
solve gives both.

---

## 2. Design variables — latent, and foreclosed by the current shape

`design` is one real number today and there is one functional, so
nothing here is a regression. It is what the architecture forecloses.

`by_tangent` and `by_adjoint` each take one vector and return one
scalar, and each calls `dense_solve`, which builds the matrix and
factorises it from scratch. `chain_by_tangent` and `chain_by_adjoint`
do the same once per block. **No factor is retained anywhere in the
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

A factorisation that persists after the solve which produced it: an object
holding the banded factors, constructed once per block per state, then
applied to as many right-hand sides as there are design variables. The
tangent then costs one factorisation plus `n_d` substitutions.

This is a design decision, not a tuning one. It needs a type that owns
factors, and `dense_solve`'s present signature - matrix in, solution
out, nothing retained - cannot express it.

---

## 3. Functionals — the same shape, transposed

`by_adjoint` has the identical structure: one `g`, one scalar returned,
one new factorisation. `n_f` functionals cost `n_f` factorisations.

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
per block. Storing `n_f` functionals means those become a second
dimension, not a loop over independent solves.

---

## 4. The gradient path has an avoidable order

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
gradient - which is per design variable once axis 2 is extended.

The fix follows axis 1: if the solver takes a banded matrix, the
gradient path passes it the band it already has, and neither the dense
array nor the `n^2`-edge graph is built.

---

## Ranked, by what binds first

1. ~~**A tolerance that stops being reachable.**~~ **Settled.** It
   followed the step rather than the count, and no tolerance was the
   control parameter: newton now stops when it stops progressing. Two and a half
   times less time at the sizes where it bound, with every functional
   unchanged. Measurements taken before this are still measuring the
   iteration limit where the step was small.
2. **A matrix allocated again every iteration.** `dense_direct_solve`
   allocates `n by n` inside itself, about 2 MB an iteration at 471
   unknowns, which is the whole of the memory growth. The shape does
   not change between iterations. Follows from 1, and low cost to fix
   on its own.
3. **Dense elimination of a banded matrix.** The solve is
   cubic; ~1250x flops available at the largest size measured, the
   ratio growing as `n^2/b^2`. The representation is not the problem -
   the block statement is 2.5 MB where the run reaches 185.
4. **No factorisation reuse.** Multiplies axis 1 by `n_d` and `n_f`.
   Latent today at one design and one functional, and it is the axis
   the user named first, so it should be settled before it is paid.
5. **Triple materialisation on the gradient path.** `n^2` edges to
   carry a matrix already held as an array. Confined to sensitivity,
   not the hot path.
6. **Krylov as the fallback.** Linear in `n` but does not converge
   on difference blocks. Superseded by 1, and should be retired as the
   remedy for large horizons once 1 is done.

Nothing here is about parallelism. Every item is a single-processor
algorithm or storage choice, and each one makes the eventual parallel
decomposition easier rather than harder: a banded solve over a chain
of blocks is what a time-parallel method partitions.


---

## 5. Measured again with the level below attached (2026-08-26)

Everything above was measured at no optimisation: neither build passed
an -O flag. That does not change any order stated - a cubic is a cubic
at -O0 - but it changes every constant by more than an order of
magnitude, and one of the two "binding" items is removed under it.

```
                                          before        after
 bdf2, 161 instants, the state              2.26 s      0.30 s
 field 8x8, 11 instants, 2112 unknowns     102 s        1.95 s
 fitted operator, 32x32 cells              179 s        0.059 s
```

Three changes did that, in the order they matter:

- `-O3` on the library, `-O2` on the application. The dense elimination
  is 98 per cent of a field march (gprof), and it was 50 times slower
  than it needed to be for lack of a flag.
- The elimination and the substitutions run down columns. The array is
  column-major, and the row-oriented loop it had traversed it against
  the storage order.
- The fitted balance appended one triple at a time to an array, which
  copies the array each time. It grows by doubling now; the assembly
  is linear, and 32x32 cells at any form degree take under a tenth of
  a second where they took three minutes.

None of this is parallelism. All of it is one processor performing the
requested operations in the order the memory is laid out.

**The field block as a statement.** A spatial mesh under every instant
multiplies the unknowns by the cell count, so the dense elimination
binds at a few thousand unknowns - 8x8 cells by 11 instants is 2112,
16x16 by 41 is 31,000, and the latter's factorisation is an hour. The
block now records its own tangent as triples (`explicit_tangent`),
so the matrix is formed at the cost of its nonzeros and not by n
applies, and newton states it as a stencil; that is what a banded
or a multigrid solve needs, and the dense one already gains from it.

**Multigrid over the space-time block does not converge with a point
smoother**, and the reason is structural, not a defect in the
composition. A derived row has one on the degree it determines and
1/dt^d on the same instant's value; a point Gauss-Seidel sweep pivots
on the one and diverges. What the structure requires is a sweep by
instant - the time coupling is lower triangular, so a sweep in
instant order is exact in time - with multigrid inside each instant
over the nodes, on the system reduced to values. That is the
implicit time-stepping solver built as a linear solver for the block,
and it is the next slice; the compiled tangent is its prerequisite
and is done.
