# Four problems from one Lagrangian

Status: implemented, 2026-10-02, but optimization under uncertainty,
which is deferred. The problems are stated in `physics/`; `src/` is the
mathematics and the computing they call.

## The rule

Each condition has a multiplier: the sensitivity of the objective to
that condition. Removing a condition sets its multiplier to zero.

## The sequence

Each problem differs from the one before it by one term of the
Lagrangian.

| Problem | Lagrangian | Computed |
|---|---|---|
| physics | L = <λ, r> + ... | the state |
| sensitivity | L = J + <λ, r> + ... + <κ, s - s0> | the state and its multipliers; κ = -dJ/ds |
| optimization | L = J + <λ, r> + ... | the design s at which dJ/ds = 0 |
| optimization under uncertainty | L = E_ξ[ J + <λ, r> + ... ] | that design, with J averaged over ξ |

The laws, the conditions and their multipliers are written once. The
multipliers supply the gradient in every problem, so that no problem
requires its own solver or its own adjoint.

## On the beam

- **Physics:** `physics/beam/euler_bernoulli_beam.f90`, cases 2, 3, 5,
  6 and 7, which have no objective; cases 6 and 7 by finite elements
  (`doc/finite-element-method.md`).
- **Sensitivity:** the Lagrangian of
  `physics/beam/euler_bernoulli_beam_optimization.f90` with `+ kappa*d`,
  `d = continuous_lagrangian(omega % design(), [s1 - s1_0, s2 - s2_0])`:
  `minimize` returns kappa = -dJ/ds at (s1_0, s2_0).
- **Optimization:** that program as it stands, without `+ kappa*d`.
  `minimize` finds s = (0.51178039, -0.31039058) at 21 x 17 by Newton's
  iteration on the design, the beam discretized by finite elements. By
  the fitted differences of order 6 it found (0.50961827, -0.31071464),
  the digits the program's own iteration gave before it was moved into
  the library.
- **Under uncertainty:** deferred. Its non-intrusive form is a weighted
  sum of Lagrangians, one at each node of the distribution; its
  intrusive form is one Lagrangian of the laws projected onto the
  basis. Each is to be stated before either is built.

## Reading a solution

| Quantity | Read |
|---|---|
| the state | `w_h = solution % fields(['w'])` |
| its derivatives along time and space, as the jet stores them | `w_h % derivative([t])`, `w_h % derivative([x, x])` |
| its derivatives along the design | `w_h % derivative([s1])` |
| a multiplier | `solution % fields(['adjoint'])`, by the name given when it was paired |
| the design | `solution % fields(['s1', 's2'])` |
| the objective and its derivatives along the design | `J % at(solution)`, then `values` |

A derivative along time or space is read, not computed: the component
the discretization stores. One that is not stored, and one asked of a
solution under a checkpoint or snapshot schedule, whose jet keeps the
values alone, are refused.

## The commits

Each verified by every digit of the 92 cases of `physics/` and
`test/graph-scheme`, serially and on 1, 2 and 4 images.

| Commit | What |
|---|---|
| d8c13d6 | the verified refactor: the expansion module, one column axis, the block of the chain, bare conditions, the design multiplier corrected |
| b2db0a8 | `physics/` for the problems, `test/graph-scheme` for the studies, each problem's `run.sh` on one process or, with `IMAGES=n`, on n images |
| cebccf0 | the design read by name |
| db74eb1 | a free design found by `minimize` |
| the commit adding this note | the derivatives of the state read by name |

## Scope

Established for one law per unknown. Systems of laws, such as velocity
and pressure, and gauge conditions require separate treatment. A free
design requires an objective, the equations paired with their
multiplier, and the expansion of order one along each design
coordinate.
