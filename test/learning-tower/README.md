# Learning Tower

A second vertical inhabitation test of the relation-centered tower.
The calculator proved the architecture can compute; this tower asks
whether the very same abstractions can **learn a parameter from
data** — with no ML framework anywhere in production.

\[
\boxed{(x,y)=(2,6),\qquad \hat y = wx,\qquad w_0=0,\qquad w_{\mathrm{learned}}=3}
\]

As a final inference check, the learned model predicts
\(x_\ast=4 \Rightarrow \hat y_\ast=12\). The tower's primary result
is the **learned parameter** \(w=3\).

## The persistent object

```text
value slots V                 operations O            ports P

w    x    ŷ    y    e         predict   error         in₁  in₂  out
│    │
w₀=0 2         6              laws appear only at Level 8
```

Structural flow (six ternary tuples, the one source of truth):

```text
w ──in₁──▶ predict        ŷ ──in₁──▶ error
x ──in₂──▶ predict        y ──in₂──▶ error
predict ──out──▶ ŷ        error ──out──▶ e
```

## The ladder

| Level | Learning interpretation | Required truth |
|---|---|---|
| 0 | symbolic domains | \(V,O,P\) distinct |
| 1 | structural flow | six ternary tuples |
| 2 | operation dependency | `predict → error` |
| 3 | relational graph | one structure |
| 4 | graph calculus | `[predict, error]` |
| 5 | data + parameter fields | \(x=2,\ y=6,\ w_0=0\) |
| 6 | trainable dependency | \(J_\Theta=\{(r,w)\}\) |
| 7 | parameter fitting | \(w\to3\) |
| 8 | model constitution | `predict = w·x`, `error = ŷ−y` |
| 9 | learning statement | learned \(w=3\); then \(x_\ast=4\Rightarrow\hat y_\ast=12\) |

Below Level 8, `predict` and `error` are only members of \(O\).
Training changes \(w\); it never mutates structure. What this tower
deliberately does **not** prove: deep networks, autodiff, gradient
descent, batches, least squares. Each of those must earn its own
tower.

## Status

```text
level 0  carriers .............. PASS
level 1  relation .............. PASS
level 2  relation algebra ...... PASS
level 3  relational graph ...... PASS
level 4  graph calculus ........ PASS
level 5  field calculus ........ PASS
level 6  discretization ........ PASS
level 7  minimization .......... PASS
level 8  constitution .......... PASS
level 9  statement ............. ABSENT
```

Data and parameters differ by domain and role, not by field class;
uncomputed slots have no fabricated values. The relation
`predict → error` existed at Level 2; only Level 4
chooses to read it as directed execution. From topology alone, the
tower derives that residual `r` depends on trainable `w` — through
`ŷ`, transitively, never by adjacency — and the reverse structure
is `J_Θᵀ`, a transpose view of the same dependency, not a
separately stored backprop graph. No derivative value exists yet;
this is the Jacobian's PATTERN, owned by the model, not by any one
data instance. Given the supplied residual `r(w) = 2w − 6`, the
existing residual minimizer changes the trainable state from
`w = 0` to `w = 3` — and, handed `(4,8)` instead of `(2,6)`, fits
`w = 2`: parameter fitting through the ordinary GMRES citizen, not
gradient descent, and the model laws are still unconstituted.
Only at constitution does `predict` become multiplication and
`error` become subtraction. With observed `{y,x}` and trainable
`w`, the structural execution computes `ŷ` and `e` **into** the
computed domain — no field on `U`, no `q(e) − law` — then reads
residual `r` from `e`, reproducing `r(w) = 2w − 6` without running
a solver; and meaning does not bend topology: the trainable slots
the evaluation reads equal Level 6's structural support `{w}`. Each level is added RED-first, one review gate at a time, with the
import gate holding every rung to its own allowlist. Production
changes are expected to be **none** at every level — that is the
experiment.
