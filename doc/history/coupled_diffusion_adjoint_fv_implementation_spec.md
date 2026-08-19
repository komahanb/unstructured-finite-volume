# Implementation Specification: Coupled Diffusion Physics, Functional, Finite-Volume Assembly, and Adjoint Gradients

## 1. Purpose

Implement a small but extensible framework where:

1. `interface_physics` provides **local physics laws** and their partial derivatives.
2. `interface_function` provides **local functional/objective densities** and their partial derivatives.
3. `finite_volume_assembler` integrates these local quantities over the mesh and assembles global residuals, Jacobians, and design partials.
4. `adjoint_solver` solves the transpose residual-Jacobian system.
5. `optimizer` uses the adjoint gradient to update design parameters.

The target demonstration problem is a coupled diffusion system on a box geometry.

The architectural law is:

\[
\boxed{
\text{interfaces define local laws; assemblers define geometric integration}
}
\]

---

## 2. Mathematical Problem

### 2.1 Domain

Use a rectangular box:

\[
\Omega = [0,L_x]\times[0,L_y]\times[0,L_z]
\]

with boundary:

\[
\partial\Omega =
\Gamma_{x0}\cup \Gamma_{xL}\cup \Gamma_{y0}\cup \Gamma_{yL}\cup \Gamma_{z0}\cup \Gamma_{zL}.
\]

---

## 3. State Variables

Use two coupled scalar fields:

\[
q =
\begin{bmatrix}
u\\
v
\end{bmatrix}
\]

where:

| Symbol | Meaning |
|---|---|
| \(u\) | first diffusing quantity |
| \(v\) | second diffusing quantity |
| \(q\) | vector state field |
| \(n_q = 2\) | number of state components |

---

## 4. Design Parameters

Use three design parameters:

\[
p =
\begin{bmatrix}
D\\
\alpha\\
\beta
\end{bmatrix}
\]

where:

| Design parameter | Meaning |
|---|---|
| \(D\) | common diffusion coefficient |
| \(\alpha\) | coupling rate from \(v\) to \(u\) |
| \(\beta\) | coupling rate from \(u\) to \(v\) |

The coupling matrix is:

\[
C(\alpha,\beta)
=
\begin{bmatrix}
-\alpha & \alpha\\
\beta & -\beta
\end{bmatrix}.
\]

The source/coupling term is:

\[
S(q,p)=Cq.
\]

Explicitly:

\[
S_u = \alpha(v-u)
\]

\[
S_v = \beta(u-v).
\]

---

## 5. Governing Residual

### 5.1 Strong Differential Form

The coupled diffusion equation is:

\[
q_t = D\nabla^2 q + Cq.
\]

Equivalently, the residual is:

\[
R(q,q_t,\nabla^2 q;p)
=
q_t - D\nabla^2q - Cq = 0.
\]

In component form:

\[
R_u = u_t - D\nabla^2u - \alpha(v-u)
\]

\[
R_v = v_t - D\nabla^2v - \beta(u-v).
\]

---

## 6. Finite-Volume Form

The finite-volume residual for cell \(i\) is:

\[
R_i =
\int_{\Omega_i} q_t\,d\Omega
-
\int_{\partial\Omega_i} F(q,\nabla q,p)\cdot n\,dS
-
\int_{\Omega_i} S(q,p)\,d\Omega.
\]

For diffusion:

\[
F = D\nabla q.
\]

Thus the discrete residual is:

\[
R_i =
V_i \dot q_i
-
\sum_{f\in\partial\Omega_i} A_f F_f\cdot n_f
-
V_i S_i.
\]

For a two-cell face \(f=(i,j)\), use a centered gradient approximation:

\[
\nabla q_f\cdot n_f
\approx
\frac{q_j-q_i}{d_{ij}}.
\]

Then the face flux is:

\[
F_f\cdot n_f
=
D\frac{q_j-q_i}{d_{ij}}.
\]

The contribution to cell \(i\) is:

\[
- A_f D\frac{q_j-q_i}{d_{ij}}.
\]

---

## 7. Boundary Conditions

Use the following initial demonstration boundary conditions.

### 7.1 Dirichlet Boundaries

On the left face:

\[
q = q_L
\quad\text{on } x=0.
\]

On the right face:

\[
q = q_R
\quad\text{on } x=L_x.
\]

where:

\[
q_L =
\begin{bmatrix}
u_L\\
v_L
\end{bmatrix},
\qquad
q_R =
\begin{bmatrix}
u_R\\
v_R
\end{bmatrix}.
\]

### 7.2 No-Flux Boundaries

On the other faces:

\[
\nabla q\cdot n = 0.
\]

That is:

\[
\nabla u\cdot n = 0,
\qquad
\nabla v\cdot n = 0.
\]

---

## 8. Functional of Interest

Use a tracking objective:

\[
J(q,p)
=
\int_0^T\int_\Omega
\phi(q,p)\,d\Omega\,dt.
\]

The local functional density is:

\[
\phi(q,p)
=
\frac12(q-q_d)^T W(q-q_d)
+
\frac{\gamma_D}{2}(D-D_0)^2
+
\frac{\gamma_\alpha}{2}(\alpha-\alpha_0)^2
+
\frac{\gamma_\beta}{2}(\beta-\beta_0)^2.
\]

Here:

| Symbol | Meaning |
|---|---|
| \(q_d\) | desired target field |
| \(W\) | positive semidefinite state-weight matrix |
| \(D_0,\alpha_0,\beta_0\) | reference design values |
| \(\gamma_D,\gamma_\alpha,\gamma_\beta\) | design regularization weights |

For steady-state implementation, drop the time integral:

\[
J(q,p)
=
\int_\Omega \phi(q,p)\,d\Omega.
\]

---

## 9. Functional Partials

The functional state partial is:

\[
J_q = W(q-q_d).
\]

The design partials are:

\[
J_D = \gamma_D(D-D_0)
\]

\[
J_\alpha = \gamma_\alpha(\alpha-\alpha_0)
\]

\[
J_\beta = \gamma_\beta(\beta-\beta_0).
\]

In local-density form:

\[
\phi_q = W(q-q_d)
\]

\[
\phi_D = \gamma_D(D-D_0)
\]

\[
\phi_\alpha = \gamma_\alpha(\alpha-\alpha_0)
\]

\[
\phi_\beta = \gamma_\beta(\beta-\beta_0).
\]

The finite-volume assembler integrates these local partials:

\[
J_q^{global} = \sum_i V_i \phi_{q,i}
\]

\[
J_p^{global} = \sum_i V_i \phi_{p,i}.
\]

If design parameters are global scalars, the integrated design partials are scalar sums.

---

## 10. Physics Partials

The physics interface should expose local partials for:

1. time term,
2. diffusive flux term,
3. source/coupling term,
4. design parameters.

---

### 10.1 Time Partial

For the transient residual:

\[
R_i^{time}=V_i\dot q_i.
\]

The local time partial is:

\[
\frac{\partial R_i}{\partial \dot q_i}=V_i I.
\]

At the physics-interface level, before finite-volume integration, expose:

\[
R_{\dot q}=I.
\]

The assembler multiplies by the cell volume.

---

### 10.2 Source Term Partials

The source term is:

\[
S(q,p)=Cq.
\]

The residual contains:

\[
-R_i^{source}=-V_iS_i.
\]

The local source state partial is:

\[
S_q=C.
\]

The residual source contribution to the state Jacobian is:

\[
\frac{\partial(-V_iS_i)}{\partial q_i}
=
-V_iC.
\]

At the physics-interface level, expose:

\[
S_q=C.
\]

The assembler applies the residual sign and the volume.

---

### 10.3 Coupling Matrix Partials

The coupling matrix is:

\[
C(\alpha,\beta)
=
\begin{bmatrix}
-\alpha & \alpha\\
\beta & -\beta
\end{bmatrix}.
\]

Its partial derivatives are:

\[
C_\alpha =
\frac{\partial C}{\partial \alpha}
=
\begin{bmatrix}
-1 & 1\\
0 & 0
\end{bmatrix}
\]

\[
C_\beta =
\frac{\partial C}{\partial \beta}
=
\begin{bmatrix}
0 & 0\\
1 & -1
\end{bmatrix}.
\]

Therefore:

\[
S_\alpha = C_\alpha q
=
\begin{bmatrix}
-v_u\\
0
\end{bmatrix}
\]

More explicitly:

\[
S_\alpha =
\begin{bmatrix}
v-u\\
0
\end{bmatrix}.
\]

and:

\[
S_\beta = C_\beta q
=
\begin{bmatrix}
0\\
u-v
\end{bmatrix}.
\]

The residual contains \(-S\), so the residual design partials from the source term are:

\[
R_\alpha^{source}
=
-
S_\alpha
=
\begin{bmatrix}
u-v\\
0
\end{bmatrix}
\]

\[
R_\beta^{source}
=
-
S_\beta
=
\begin{bmatrix}
0\\
v-u
\end{bmatrix}.
\]

After integration over a cell:

\[
R_{\alpha,i}^{source}
=
-V_i S_{\alpha,i}
\]

\[
R_{\beta,i}^{source}
=
-V_i S_{\beta,i}.
\]

---

### 10.4 Diffusive Flux Partials

For an internal face \(f=(i,j)\):

\[
F_f\cdot n_f =
D\frac{q_j-q_i}{d_{ij}}.
\]

The cell residual receives:

\[
R_i^{face}
=
-A_fD\frac{q_j-q_i}{d_{ij}}.
\]

Therefore:

\[
\frac{\partial R_i^{face}}{\partial q_i}
=
+\frac{A_fD}{d_{ij}}I
\]

\[
\frac{\partial R_i^{face}}{\partial q_j}
=
-\frac{A_fD}{d_{ij}}I.
\]

The design partial with respect to \(D\) is:

\[
\frac{\partial R_i^{face}}{\partial D}
=
-A_f\frac{q_j-q_i}{d_{ij}}.
\]

The neighboring cell \(j\) receives the opposite contribution.

---

## 11. Global Discrete Residual

After assembly, the discrete residual has the form:

\[
R_h(q,\dot q,p)
=
M\dot q - K(D)q - S_h(q,\alpha,\beta)
\]

where:

| Symbol | Meaning |
|---|---|
| \(M\) | finite-volume mass matrix, usually diagonal cell-volume matrix |
| \(K(D)\) | diffusion operator assembled from face fluxes |
| \(S_h\) | integrated source/coupling term |

For a steady-state problem:

\[
R_h(q,p) = -K(D)q - S_h(q,\alpha,\beta).
\]

Depending on sign convention, one may also write:

\[
R_h(q,p)=K(D)q - S_h(q,\alpha,\beta).
\]

The implementation must keep one residual sign convention consistently.

This specification uses:

\[
R_h =
M\dot q
-
K_{flux}(D)q
-
S_h(q,p).
\]

---

## 12. Global Residual Jacobians

The global state Jacobian is:

\[
R_q =
\frac{\partial R_h}{\partial q}.
\]

For the steady-state diffusion-coupling problem:

\[
R_q
=
-K_q - S_q.
\]

For the transient problem:

\[
R_q
=
-K_q - S_q
\]

and:

\[
R_{\dot q}=M.
\]

If using implicit time integration with time step \(\Delta t\), the solve Jacobian is:

\[
A =
\frac{1}{\Delta t}M + R_q.
\]

For backward Euler:

\[
R^{n+1}
=
M\frac{q^{n+1}-q^n}{\Delta t}
-
Kq^{n+1}
-
S(q^{n+1},p).
\]

Then:

\[
\frac{\partial R^{n+1}}{\partial q^{n+1}}
=
\frac{1}{\Delta t}M
-
K_q
-
S_q.
\]

---

## 13. Residual Design Partials

The global residual design partial is:

\[
R_p =
\frac{\partial R_h}{\partial p}.
\]

For:

\[
p=[D,\alpha,\beta]^T
\]

the columns are:

\[
R_p =
\begin{bmatrix}
R_D & R_\alpha & R_\beta
\end{bmatrix}.
\]

### 13.1 Diffusion Design Partial

From the diffusive face term:

\[
R_D^{face}
=
-A_f\frac{q_j-q_i}{d_{ij}}.
\]

Assemble over all faces.

### 13.2 Coupling Design Partials

From the source term:

\[
R_\alpha^{cell}
=
-V_i
\begin{bmatrix}
v_i-u_i\\
0
\end{bmatrix}
=
V_i
\begin{bmatrix}
u_i-v_i\\
0
\end{bmatrix}
\]

\[
R_\beta^{cell}
=
-V_i
\begin{bmatrix}
0\\
u_i-v_i
\end{bmatrix}
=
V_i
\begin{bmatrix}
0\\
v_i-u_i
\end{bmatrix}.
\]

---

## 14. Adjoint Method

The constrained optimization problem is:

\[
\min_p J(q,p)
\]

subject to:

\[
R(q,p)=0.
\]

Define the Lagrangian:

\[
\mathcal{L}(q,p,\lambda)
=
J(q,p)+\lambda^TR(q,p).
\]

The adjoint equation is:

\[
R_q^T\lambda = -J_q^T.
\]

After solving for \(\lambda\), compute the total gradient:

\[
\frac{dJ}{dp}
=
J_p + R_p^T\lambda.
\]

Component-wise:

\[
\frac{dJ}{dD}
=
J_D + R_D^T\lambda
\]

\[
\frac{dJ}{d\alpha}
=
J_\alpha + R_\alpha^T\lambda
\]

\[
\frac{dJ}{d\beta}
=
J_\beta + R_\beta^T\lambda.
\]

---

## 15. Required Components

Implement the following components.

---

# 15.1 `interface_physics`

## Responsibility

`interface_physics` provides local physics laws and their local partial derivatives.

It must not know:

- global sparse matrix structure,
- mesh connectivity beyond local cell/face data passed to it,
- optimizer,
- adjoint solve,
- global integration rules.

## Required Methods

### `num_state_components()`

Returns:

```text
2
```

### `num_design_parameters()`

Returns:

```text
3
```

### `source(q, p)`

Computes:

\[
S(q,p)=Cq.
\]

Input:

```text
q: local state vector of size 2
p: design vector [D, alpha, beta]
```

Output:

```text
S: source vector of size 2
```

Expected result:

```text
S[0] = alpha * (v - u)
S[1] = beta  * (u - v)
```

---

### `source_state_partial(q, p)`

Computes:

\[
S_q = C.
\]

Output:

```text
2 x 2 matrix
```

Expected result:

```text
[[-alpha,  alpha],
 [ beta,  -beta ]]
```

---

### `source_design_partial(q, p)`

Computes:

\[
S_p =
\begin{bmatrix}
S_D & S_\alpha & S_\beta
\end{bmatrix}.
\]

Since the source does not depend on \(D\):

\[
S_D=0.
\]

Expected output shape:

```text
2 x 3
```

Expected result:

```text
S_p[:,0] = [0, 0]
S_p[:,1] = [v - u, 0]
S_p[:,2] = [0, u - v]
```

---

### `diffusive_flux(q_left, q_right, face_geometry, p)`

Computes the diffusive normal flux:

\[
F_f\cdot n_f = D\frac{q_R-q_L}{d}.
\]

Input:

```text
q_left: state on owner cell
q_right: state on neighbor or boundary ghost value
face_geometry:
    area
    distance
    normal
p:
    [D, alpha, beta]
```

Output:

```text
normal_flux: vector of size 2
```

Expected result:

```text
normal_flux = D * (q_right - q_left) / distance
```

The area multiplication is owned by the assembler, not the physics interface.

---

### `diffusive_flux_state_partials(q_left, q_right, face_geometry, p)`

Computes:

\[
\frac{\partial F}{\partial q_L}
=
-\frac{D}{d}I
\]

\[
\frac{\partial F}{\partial q_R}
=
+\frac{D}{d}I.
\]

Output:

```text
dflux_dq_left:  2 x 2
dflux_dq_right: 2 x 2
```

---

### `diffusive_flux_design_partial(q_left, q_right, face_geometry, p)`

Computes:

\[
\frac{\partial F}{\partial D}
=
\frac{q_R-q_L}{d}.
\]

Since the flux does not depend on \(\alpha\) or \(\beta\):

\[
F_\alpha=0,
\qquad
F_\beta=0.
\]

Expected output shape:

```text
2 x 3
```

Expected result:

```text
F_p[:,0] = (q_right - q_left) / distance
F_p[:,1] = [0, 0]
F_p[:,2] = [0, 0]
```

---

# 15.2 `interface_function`

## Responsibility

`interface_function` provides local functional density and partial derivatives.

It must not know:

- global sparse matrices,
- adjoint vectors,
- optimization algorithms,
- finite-volume connectivity.

## Required Methods

### `density(q, p, q_desired, weights)`

Computes:

\[
\phi(q,p)
=
\frac12(q-q_d)^TW(q-q_d)
+
\frac{\gamma_D}{2}(D-D_0)^2
+
\frac{\gamma_\alpha}{2}(\alpha-\alpha_0)^2
+
\frac{\gamma_\beta}{2}(\beta-\beta_0)^2.
\]

Output:

```text
scalar phi
```

---

### `state_partial(q, p, q_desired, weights)`

Computes:

\[
\phi_q = W(q-q_d).
\]

Output:

```text
vector of size 2
```

---

### `design_partial(q, p, q_desired, weights)`

Computes:

\[
\phi_p =
\begin{bmatrix}
\gamma_D(D-D_0)\\
\gamma_\alpha(\alpha-\alpha_0)\\
\gamma_\beta(\beta-\beta_0)
\end{bmatrix}.
\]

Output:

```text
vector of size 3
```

---

# 15.3 `finite_volume_assembler`

## Responsibility

The finite-volume assembler owns:

- mesh traversal,
- cell-volume integration,
- face-area integration,
- boundary condition treatment,
- sparse matrix assembly,
- global residual assembly,
- global residual Jacobian assembly,
- global residual design partial assembly,
- global functional partial assembly.

It must not define the physics law or objective law.

---

## Required Inputs

```text
mesh
physics: interface_physics
function: interface_function
state q
state_time_derivative qdot
design p
boundary_conditions
functional_target q_desired
functional_weights
```

---

## Required Global Outputs

The assembler must compute:

```text
R      : global residual vector
Rq     : global residual state Jacobian
Rqdot  : global residual time Jacobian
Rp     : global residual design partial matrix
J      : scalar functional value
Jq     : global functional state partial vector
Jp     : global functional design partial vector
```

Shapes:

```text
num_cells = Nc
num_state_components = ns = 2
num_design_parameters = np = 3

R     shape: (Nc*ns)
Rq    shape: (Nc*ns, Nc*ns)
Rqdot shape: (Nc*ns, Nc*ns)
Rp    shape: (Nc*ns, np)
J     shape: scalar
Jq    shape: (Nc*ns)
Jp    shape: (np)
```

---

## 16. Assembly Rules

### 16.1 Cell Time Term

For each cell \(i\):

\[
R_i \mathrel{+}= V_i \dot q_i.
\]

Jacobian contribution:

\[
R_{\dot q,ii} \mathrel{+}= V_i I.
\]

---

### 16.2 Cell Source Term

For each cell \(i\):

\[
R_i \mathrel{-}= V_i S(q_i,p).
\]

State Jacobian:

\[
R_{q,ii} \mathrel{-}= V_i S_q(q_i,p).
\]

Design partial:

\[
R_{p,i} \mathrel{-}= V_i S_p(q_i,p).
\]

---

### 16.3 Internal Face Diffusion Term

For an internal face \(f=(i,j)\):

1. Evaluate flux from owner \(i\) to neighbor \(j\):

\[
F_f = D\frac{q_j-q_i}{d_{ij}}.
\]

2. Residual contribution:

\[
R_i \mathrel{-}= A_fF_f
\]

\[
R_j \mathrel{+}= A_fF_f.
\]

3. State Jacobian contribution:

The physics interface provides:

\[
F_{q_i}=-\frac{D}{d_{ij}}I
\]

\[
F_{q_j}=+\frac{D}{d_{ij}}I.
\]

Assembler applies signs:

For owner cell \(i\):

\[
R_{q,ii} \mathrel{-}= A_f F_{q_i}
\]

\[
R_{q,ij} \mathrel{-}= A_f F_{q_j}
\]

For neighbor cell \(j\):

\[
R_{q,ji} \mathrel{+}= A_f F_{q_i}
\]

\[
R_{q,jj} \mathrel{+}= A_f F_{q_j}.
\]

4. Design partial:

The physics interface provides:

\[
F_p =
\begin{bmatrix}
(q_j-q_i)/d_{ij} & 0 & 0
\end{bmatrix}.
\]

Assembler applies:

\[
R_{p,i} \mathrel{-}= A_f F_p
\]

\[
R_{p,j} \mathrel{+}= A_f F_p.
\]

---

### 16.4 Dirichlet Boundary Face

For a boundary face with prescribed value \(q_b\):

\[
F_b = D\frac{q_b-q_i}{d_{ib}}.
\]

Residual contribution:

\[
R_i \mathrel{-}= A_bF_b.
\]

State Jacobian:

\[
F_{q_i}=-\frac{D}{d_{ib}}I.
\]

So:

\[
R_{q,ii} \mathrel{-}= A_bF_{q_i}.
\]

Design partial:

\[
F_D=\frac{q_b-q_i}{d_{ib}}.
\]

So:

\[
R_{p,i} \mathrel{-}= A_bF_p.
\]

If the boundary value itself depends on design, add:

\[
F_p \mathrel{+}= \frac{\partial F}{\partial q_b}\frac{\partial q_b}{\partial p}.
\]

For the first implementation, assume boundary values are design-independent.

---

### 16.5 No-Flux Boundary Face

For no-flux boundary:

\[
F_b=0.
\]

No residual, Jacobian, or design partial contribution is added.

---

## 17. Functional Assembly Rules

For each cell \(i\):

1. Compute local density:

\[
\phi_i=\phi(q_i,p).
\]

2. Add to global objective:

\[
J \mathrel{+}= V_i\phi_i.
\]

3. Add state partial:

\[
J_{q,i} \mathrel{+}= V_i\phi_{q,i}.
\]

4. Add design partial:

\[
J_p \mathrel{+}= V_i\phi_{p,i}.
\]

Important:

If \(p\) is global, then \(J_p\) is accumulated globally.

If \(p\) is cell-local in a future extension, then \(J_p\) should be assembled into a distributed parameter vector.

---

## 18. Solver Workflow

### 18.1 Primal Solve

Given \(p\), solve:

\[
R(q,p)=0.
\]

For steady state:

```text
while not converged:
    assemble R, Rq
    solve Rq * dq = -R
    q = q + dq
```

For transient backward Euler:

```text
for each time step:
    while not converged:
        assemble R, A = (1/dt)M + Rq
        solve A * dq = -R
        q_new = q_new + dq
```

---

### 18.2 Functional and Partials

After primal convergence:

```text
assemble J, Jq, Jp, Rq, Rp
```

---

### 18.3 Adjoint Solve

Solve:

\[
R_q^T\lambda = -J_q^T.
\]

Implementation:

```text
solve transpose(Rq) * lambda = -Jq
```

For implicit transient problems, use the correct time-discrete residual Jacobian.

---

### 18.4 Gradient Evaluation

Compute:

\[
g = J_p + R_p^T\lambda.
\]

Implementation:

```text
gradient = Jp + transpose(Rp) * lambda
```

Expected shape:

```text
gradient shape: (num_design_parameters)
```

---

### 18.5 Optimization Loop

```text
initialize design p

for optimization_iteration:
    solve primal R(q,p)=0
    assemble J, Jq, Jp, Rq, Rp
    solve adjoint transpose(Rq) lambda = -Jq
    gradient = Jp + transpose(Rp) lambda
    update p using optimizer
```

The optimizer can initially be gradient descent:

```text
p = p - step_size * gradient
```

Later extensions may use L-BFGS, SQP, IPOPT, or custom constrained optimizers.

---

## 19. Suggested File Structure

```text
src/
  interfaces/
    interface_physics.*
    interface_function.*

  physics/
    coupled_diffusion_physics.*

  functions/
    tracking_function.*

  mesh/
    mesh.*
    box_mesh.*
    face.*
    cell.*

  assembly/
    finite_volume_assembler.*

  solvers/
    primal_solver.*
    adjoint_solver.*
    linear_solver.*

  optimization/
    optimizer.*
    gradient_descent_optimizer.*

  drivers/
    run_coupled_diffusion_optimization.*

tests/
  test_physics_partials.*
  test_function_partials.*
  test_fv_residual_assembly.*
  test_fv_jacobian_fd_check.*
  test_adjoint_gradient_fd_check.*
```

Use the language-specific extension appropriate for the implementation target.

---

## 20. Data Structures

### 20.1 State Storage

Use a cell-major layout:

```text
q[cell_id, component_id]
```

Flattening rule:

```text
global_id = cell_id * num_state_components + component_id
```

For two fields:

```text
u_i = q[i,0]
v_i = q[i,1]
```

---

### 20.2 Design Storage

Use:

```text
p[0] = D
p[1] = alpha
p[2] = beta
```

---

### 20.3 Face Geometry

Each face should provide:

```text
owner_cell
neighbor_cell
area
distance
normal
boundary_marker
```

For boundary faces:

```text
neighbor_cell = none
boundary_marker = x_min, x_max, y_min, y_max, z_min, z_max
```

---

## 21. Boundary Condition Objects

Boundary condition type:

```text
Dirichlet
NoFlux
```

Boundary condition data:

```text
boundary_marker
type
value(component)
```

Example:

```text
x_min:
    type: Dirichlet
    value: [u_L, v_L]

x_max:
    type: Dirichlet
    value: [u_R, v_R]

y_min:
    type: NoFlux

y_max:
    type: NoFlux

z_min:
    type: NoFlux

z_max:
    type: NoFlux
```

---

## 22. Minimal Numerical Example

Use a 1D box first:

\[
\Omega=[0,1].
\]

Parameters:

```text
Lx = 1.0
num_cells = 50

D     = 0.05
alpha = 1.0
beta  = 0.5

u_left  = 1.0
v_left  = 0.0

u_right = 0.0
v_right = 1.0
```

Functional target:

```text
q_desired[:,0] = 0.5
q_desired[:,1] = 0.5
```

Weights:

```text
W = identity(2)

D0     = 0.05
alpha0 = 1.0
beta0  = 0.5

gamma_D     = 1.0e-4
gamma_alpha = 1.0e-4
gamma_beta  = 1.0e-4
```

---

## 23. Acceptance Tests

### 23.1 Physics Source Test

Given:

```text
q = [2.0, 5.0]
p = [0.1, 3.0, 4.0]
```

Expected:

```text
S[0] = 3.0 * (5.0 - 2.0) = 9.0
S[1] = 4.0 * (2.0 - 5.0) = -12.0
```

---

### 23.2 Source State Partial Test

Expected:

```text
S_q =
[[-3.0,  3.0],
 [ 4.0, -4.0]]
```

Verify against finite differences.

---

### 23.3 Source Design Partial Test

Expected:

```text
S_D     = [0.0, 0.0]
S_alpha = [3.0, 0.0]
S_beta  = [0.0, -3.0]
```

Verify against finite differences.

---

### 23.4 Flux Test

Given:

```text
q_left  = [1.0, 2.0]
q_right = [3.0, 5.0]
D = 0.1
distance = 0.5
```

Expected:

```text
flux = 0.1 * ([3.0,5.0] - [1.0,2.0]) / 0.5
     = [0.4, 0.6]
```

---

### 23.5 Flux State Partial Test

Expected:

```text
dflux_dq_left  = -(0.1/0.5) * I = -0.2 I
dflux_dq_right = +(0.1/0.5) * I = +0.2 I
```

Verify against finite differences.

---

### 23.6 Flux Design Partial Test

Expected:

```text
dflux_dD = ([3.0,5.0] - [1.0,2.0]) / 0.5
         = [4.0, 6.0]
```

and:

```text
dflux_dalpha = [0.0, 0.0]
dflux_dbeta  = [0.0, 0.0]
```

---

### 23.7 Functional State Partial Test

Given:

```text
q = [2.0, 5.0]
q_desired = [1.0, 1.0]
W = identity(2)
```

Expected:

```text
phi_q = [1.0, 4.0]
```

Verify against finite differences.

---

### 23.8 Functional Design Partial Test

Given:

```text
D = 0.2
D0 = 0.1
gamma_D = 10.0
```

Expected:

```text
phi_D = 10.0 * (0.2 - 0.1) = 1.0
```

Similarly check \(\alpha\) and \(\beta\).

---

### 23.9 Global Jacobian Finite-Difference Test

After assembling \(R(q,p)\) and \(R_q\), choose a random perturbation \(\delta q\).

Check:

\[
R_q\delta q
\approx
\frac{R(q+\epsilon\delta q,p)-R(q,p)}{\epsilon}
\]

with:

```text
epsilon = 1.0e-6
```

Expected relative error:

```text
< 1.0e-6 to 1.0e-5
```

depending on floating-point precision and solver details.

---

### 23.10 Global Design Partial Finite-Difference Test

Choose random perturbation \(\delta p\).

Check:

\[
R_p\delta p
\approx
\frac{R(q,p+\epsilon\delta p)-R(q,p)}{\epsilon}
\]

Expected relative error:

```text
< 1.0e-6 to 1.0e-5
```

---

### 23.11 Adjoint Gradient Test

After solving the primal and adjoint equations, compute:

\[
g_{adj}=J_p+R_p^T\lambda.
\]

Compare with finite-difference total derivative:

\[
g_k
\approx
\frac{J(q(p+\epsilon e_k),p+\epsilon e_k)-J(q(p),p)}{\epsilon}.
\]

Important:

For each perturbed design parameter, re-solve the primal equation.

Expected relative error:

```text
< 1.0e-5
```

---

## 24. Sign Convention Warning

The implementation must choose and preserve a single residual sign convention.

This document uses:

\[
R_i =
V_i\dot q_i
-
\sum_f A_fF_f
-
V_iS_i.
\]

With:

\[
F_f = D\frac{q_j-q_i}{d_{ij}}.
\]

Therefore:

\[
R_q^T\lambda = -J_q^T
\]

and:

\[
g = J_p + R_p^T\lambda.
\]

If the residual sign is changed, the adjoint equation and gradient expression must be checked accordingly.

---

## 25. Extension Points

The design should allow future extensions:

### 25.1 More State Variables

Allow:

\[
q\in\mathbb{R}^{n_q}
\]

instead of hard-coding \(n_q=2\).

---

### 25.2 More Design Parameters

Allow:

\[
p\in\mathbb{R}^{n_p}
\]

instead of hard-coding \(n_p=3\).

---

### 25.3 Tensor Diffusion

Replace scalar \(D\) with tensor diffusion:

\[
F = \mathbf{D}\nabla q.
\]

---

### 25.4 Nonlinear Sources

Allow:

\[
S=S(q,p)
\]

with nonlinear state and design dependence.

---

### 25.5 Boundary Functionals

Allow:

\[
J =
\int_\Omega \phi(q,p)d\Omega
+
\int_{\Gamma} \psi(q,p)d\Gamma.
\]

Then `interface_function` should provide:

```text
boundary_density()
boundary_state_partial()
boundary_design_partial()
```

and the finite-volume assembler should integrate these using face areas.

---

### 25.6 Time-Dependent Adjoint

For transient problems, the adjoint must be derived from the fully discrete time-marching residual.

For backward Euler, each time step residual is:

\[
R^{n+1}
=
M\frac{q^{n+1}-q^n}{\Delta t}
-
Kq^{n+1}
-
S(q^{n+1},p).
\]

The global time-coupled adjoint system is block lower/upper triangular depending on ordering.

The initial implementation may restrict itself to steady-state optimization.

---

## 26. First Implementation Target

The first working version should implement:

1. 1D uniform finite-volume mesh.
2. Two-state coupled diffusion.
3. Scalar global design vector \(p=[D,\alpha,\beta]\).
4. Dirichlet boundary at left and right.
5. Steady-state residual.
6. Tracking functional.
7. Newton primal solve.
8. Adjoint solve.
9. Adjoint gradient.
10. Finite-difference gradient verification.

After this passes, extend to:

1. 2D box,
2. 3D box,
3. transient solve,
4. nonuniform mesh,
5. nonlinear sources,
6. tensor diffusion,
7. boundary functionals.

---

## 27. Implementation Checklist

### Physics Interface

- [ ] Implement `num_state_components`.
- [ ] Implement `num_design_parameters`.
- [ ] Implement `source`.
- [ ] Implement `source_state_partial`.
- [ ] Implement `source_design_partial`.
- [ ] Implement `diffusive_flux`.
- [ ] Implement `diffusive_flux_state_partials`.
- [ ] Implement `diffusive_flux_design_partial`.

### Function Interface

- [ ] Implement `density`.
- [ ] Implement `state_partial`.
- [ ] Implement `design_partial`.

### Mesh

- [ ] Implement cells.
- [ ] Implement faces.
- [ ] Implement volumes.
- [ ] Implement areas.
- [ ] Implement distances.
- [ ] Implement boundary markers.

### Assembler

- [ ] Assemble residual.
- [ ] Assemble state Jacobian.
- [ ] Assemble time Jacobian.
- [ ] Assemble residual design partial.
- [ ] Assemble functional value.
- [ ] Assemble functional state partial.
- [ ] Assemble functional design partial.

### Solver

- [ ] Implement primal Newton solve.
- [ ] Implement linear solve.
- [ ] Implement adjoint solve.
- [ ] Implement gradient evaluation.

### Tests

- [ ] Test local physics partials.
- [ ] Test local functional partials.
- [ ] Test global residual Jacobian.
- [ ] Test global residual design partial.
- [ ] Test adjoint gradient against finite differences.

---

## 28. Core Principle

The framework must preserve this ownership law:

\[
\boxed{
\texttt{interface\_physics}
\;\text{owns local residual laws}
}
\]

\[
\boxed{
\texttt{interface\_function}
\;\text{owns local objective laws}
}
\]

\[
\boxed{
\texttt{finite\_volume\_assembler}
\;\text{owns spatial integration and sparse accumulation}
}
\]

\[
\boxed{
\texttt{adjoint\_solver}
\;\text{owns transpose sensitivity propagation}
}
\]

\[
\boxed{
\texttt{optimizer}
\;\text{owns design updates}
}
\]

This is the main design rule the implementation agent should follow.
