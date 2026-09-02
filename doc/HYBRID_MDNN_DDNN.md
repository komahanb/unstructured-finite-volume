# Hybrid MDNN/DDNN Architecture: Unifying Model-Driven and Data-Driven Learning

## Overview

This framework unifies **Model-Driven Neural Networks (MDNN)** and **Data-Driven Neural Networks (DDNN)** into a single, coherent hybrid system. The central observation: the fractal graph's binary tree structure encodes this duality.

```
graph (any level: expansion, sweep, block, slice, component)
│
├─ branch(1): DDNN
│             Learned coefficients, corrections, adaptive weights
│             Data-driven adaptation
│
└─ branch(2): MDNN
              Differential operators, stencils, time integration
              Model structure (guaranteed to be satisfied)
```

---

## The Duality

### MDNN (Model-Driven): The Basis

**What it provides**: Mathematical structure that must be preserved

- Differential equations (∂u/∂t, ∇²u, etc.)
- Finite volume stencils (discretization)
- Time integration schemes (BDF, Adams, DIRK)
- Boundary and initial conditions
- Conservation laws (mass, momentum, energy)
- Physical constraints (positivity, monotonicity, etc.)

**Guarantees**: 
- Interpretable (every parameter is a stencil weight, solver tolerance, time step)
- Convergent (numerical analysis bounds apply)
- Conservative (laws built in)
- Transferable (structure transfers across meshes, parameters don't)

**Limitation**: A fixed model cannot adapt to data outside its model

### DDNN (Data-Driven): The Coefficients

**What it provides**: Adaptive learning from measurements or simulations

- Learned parameter values (diffusion coefficient k, damping μ, etc.)
- Correction terms (model error compensation)
- Basis function weights (which modes matter)
- Adaptive refinement (where model fails)
- Uncertainty quantification (parameter distributions)

**Advantages**: 
- Learns from data (no manual tuning)
- Adapts to measured data (unmodeled phenomena)
- Flexible (no fixed assumptions)
- Data-efficient (only learns unknowns)

**Limitation**: Black-box NN weights are uninterpretable; needs large data sets for accuracy

---

## Unified Architecture

### Binary Structure of Fractal Graph

The `graph` type from `graph_fractal.f90` has exactly two branches:

```fortran
type :: graph
   type(branch) :: branch(2)  ! ← Always binary
   type(token), private :: identity
end type graph
```

This binary structure is **not accidental**: it is the representation of the MDNN ↔ DDNN pair:

```
Expansion (root graph)
├─ branch(1): DDNN side
│             └─ Sweep (derivative order)
│                └─ Horizon (time blocks)
│                   └─ Block (one integrator)
│                      └─ Slice (one time step)
│                         └─ Component (stencil weights to learn)
│
└─ branch(2): MDNN side
              └─ Sweep (derivative structure)
                 └─ Horizon (time integration scheme)
                    └─ Block (one solver step)
                       └─ Slice (forward evaluation)
                          └─ Component (stencil application)
```

### Recursive Composition

At every level, the same pattern repeats:

| Level | MDNN (branch 2) | DDNN (branch 1) |
|-------|---|---|
| **Expansion** | Full PDE system | Global learned corrections |
| **Sweep** | Derivative order structure | Learned sensitivities |
| **Horizon** | Time block sequencing | Learned block-specific params |
| **Block** | Time integrator scheme | Learned step sizes |
| **Slice** | One step's stencil | Learned stencil weights |
| **Component** | Gradient/divergence/Laplacian | Learned coefficients per component |

---

## Evaluation: Forward Pass

**Pattern**: Combine MDNN structure with DDNN coefficients at each level

```fortran
recursive function evaluate_hybrid(g) result(output)
  type(graph), intent(in) :: g
  real(dp), allocatable :: output(:)
  
  real(dp), allocatable :: mdnn_output(:), ddnn_output(:)
  
  ! MDNN (branch 2): Evaluate the model structure
  ! Applies stencil, time integrator, differential operator
  mdnn_output = evaluate_mdnn(g % branch(2) % known_())
  
  ! DDNN (branch 1): Evaluate learned modifications
  ! Produces corrections, learned coefficients, basis weights
  ddnn_output = evaluate_ddnn(g % branch(1) % known_())
  
  ! Combine: structure + learned parameters
  output = mdnn_output + ddnn_output
  
end function evaluate_hybrid
```

**Semantics**:
- **MDNN computes**: The exact discretization of the PDE
- **DDNN computes**: Corrections, parameter adjustments, learned features
- **Result**: Model-respecting solution with data-informed parameters

---

## Backward Pass: Hybrid Adjoint

**Pattern**: Backpropagate through both branches simultaneously

```fortran
recursive subroutine adjoint_hybrid(g, adjoint_in, gradient_mdnn, gradient_ddnn)
  type(graph), intent(in) :: g
  real(dp), intent(in) :: adjoint_in(:)
  real(dp), intent(out) :: gradient_mdnn(:), gradient_ddnn(:)
  
  ! Backprop through MDNN (branch 2)
  ! Computes: ∂loss/∂(stencil weights, time steps, solver params)
  call adjoint_mdnn(g % branch(2) % known_(), adjoint_in, gradient_mdnn)
  
  ! Backprop through DDNN (branch 1)
  ! Computes: ∂loss/∂(learned coefficients, corrections, basis)
  call adjoint_ddnn(g % branch(1) % known_(), adjoint_in, gradient_ddnn)
  
end subroutine adjoint_hybrid
```

**Information flow**:
- The same adjoint vector is propagated through both branches
- Each computes its own gradient
- Gradients are independent and additive

---

## Three Levels of Collaboration

### Level 1: Parameter Learning (Simplest)

**Use case**: Known PDE structure, unknown parameters

```
MDNN provides: ∂u/∂t = k·∇²u (structure, but k unknown)
DDNN learns:  k from data

Example: Diffusion with unknown diffusion coefficient
  Train: minimize ||u_model(k_learned) - u_measured||²
  Result: Interpretable learned k
```

### Level 2: Correction Terms (Hybrid)

**Use case**: Accurate model, but incomplete

```
MDNN provides: ∂u/∂t = k·∇²u + f(x,t) (base model + unknown correction)
DDNN learns:  f(x,t) (model error correction)

Example: Reaction-diffusion with unknown reaction rate
  Train: minimize ||∂u/∂t - k·∇²u - f(u,v)||²
  Result: Model + learned correction term
```

### Level 3: Adaptive Basis (Most Flexible)

**Use case**: Structure known, discretization uncertain

```
MDNN provides: Stencil operators, time integration (structure)
DDNN learns:  Basis functions, refinement strategy, local metrics

Example: Multi-scale problem
  MDNN: Global solver structure (time stepping, boundary conditions)
  DDNN: Local basis selection (which stencil per region, which timestep per block)
  Result: Adaptive discretization that respects structure
```

---

## Implementation Pattern

### Attached to Any Graph Node

```fortran
type :: hybrid_node
   type(graph) :: structure_branch      ! MDNN (branch 2)
   class(neural_network) :: adapt_branch  ! DDNN (branch 1)
end type hybrid_node
```

### Evaluation at Each Level

```fortran
! At block level: compute solution using both MDNN + DDNN
call net_mdnn % forward(physics, schemes, instants, dt)  ! Structure
call net_ddnn % forward(measurements)                    ! Learned params

solution = net_mdnn % functional() + net_ddnn % functional()
```

### Training Loop

```fortran
do iteration = 1, num_iterations
   ! 1. Forward: evaluate both MDNN and DDNN
   call net_mdnn % forward(...)
   call net_ddnn % forward(...)
   
   ! 2. Compute loss
   loss = compute_loss(net_mdnn, net_ddnn, target)
   
   ! 3. Backward: gradients flow through both
   call net_mdnn % backward(...)
   call net_ddnn % backward(...)
   
   ! 4. Update: only DDNN weights change; MDNN structure stays
   grad_mdnn = net_mdnn % gradient()      ! For analysis
   grad_ddnn = net_ddnn % gradient()      ! For training
   
   ddnn_params = ddnn_params - learning_rate * grad_ddnn
end do
```

---

## Concrete Example: Van der Pol with Uncertain Damping

**Problem**: d²u/dt² + μ·(1-u²)·du/dt + u = 0, where μ is unknown

### MDNN Side (Structure)

```fortran
! Define the exact model structure
mdnn_eq = d('u', t=2) + variable('mu') * (1.0_dp - variable('u')**2) * d('u', t=1) + variable('u')

! Time integration: BDF or Adams
! Stencils: gradient, Laplacian
! Boundary conditions: fixed
! → All deterministic, all interpretable
```

### DDNN Side (Learning)

```fortran
! From measurements with measurement error, learn mu
measurements = read_experimental_data()

! Neural network learns: mu(t) or just constant mu
net_ddnn = neural_network(
   input=measurements,
   output=learned_mu,
   structure='learn_scalar'
)

call net_ddnn % train(mdnn_eq, measurements, learning_rate=0.01_dp)
```

### Evaluation (Hybrid)

```fortran
! Forward pass
call net_mdnn % forward(vdp_physics, schemes, instants, dt)
   ! Uses MDNN structure with nominal μ

call net_ddnn % forward(measurements)
   ! Produces learned correction to μ

! Combine
mu_effective = mu_nominal + net_ddnn % learned_mu
solution = solve_vdp(mu_effective)

! Backward pass
call net_mdnn % backward()
   ! ∂loss/∂(time step, solver tolerance, mesh)

call net_ddnn % backward()
   ! ∂loss/∂(learned μ, correction term)
```

### Result

- **Interpretable**: μ_learned is a single number, directly verifiable against physics
- **Data-efficient**: Only μ is learned; everything else is structure
- **Robust**: If measurements fail, μ defaults to physics-based prior
- **Transferable**: μ learned on one mesh transfers to another

---

## Example: Multi-Fidelity Learning

**Scenario**: Coarse simulation and fine measurement data

```fortran
! MDNN coarse: Fast model on coarse mesh
model_coarse_mdnn = define_diffusion_pde(mesh_coarse)

! DDNN learns: Correction from coarse to fine
net_ddnn = neural_network(
   input=model_coarse_mdnn % solution(),
   output=correction_to_fine,
   structure='learn_correction'
)

! MDNN fine: Exact model on fine mesh (reference)
model_fine_mdnn = define_diffusion_pde(mesh_fine)

! Training: minimize ||coarse + learned_correction - fine||²
do iteration = 1, num_iterations
   call model_coarse_mdnn % forward(...)
   call net_ddnn % forward(model_coarse_mdnn)
   call model_fine_mdnn % forward(...)
   
   loss = ||model_coarse_mdnn + net_ddnn - model_fine_mdnn||²
   
   call net_ddnn % backward()
   learned_params = learned_params - alpha * net_ddnn % gradient()
end do

! Deployment: use the coarse model + low-cost correction
prediction = model_coarse_mdnn + net_ddnn
```

---

## Advantages Over Pure MDNN or Pure DDNN

| Property | Pure MDNN | Pure DDNN | MDNN + DDNN |
|----------|----------|----------|-----------|
| **Interpretability** | ✓✓✓ (exact) | ✗ (black-box) | ✓✓ (structure visible) |
| **Data required** | ✗ (none) | ✓✓ (large) | ✓ (only unknowns) |
| **Robustness** | ✓ (guaranteed) | ✗ (extrapolation fails) | ✓✓ (structure + learning) |
| **Accuracy** | Limited by model | ✓✓ (if data adequate) | ✓✓ (both) |
| **Computational cost** | Low (exact) | High (NN overhead) | Low (minimal NN) |
| **Physics enforcement** | ✓ (by design) | ✗ (learned) | ✓ (by design) |
| **Adaptability** | ✗ (fixed) | ✓ (flexible) | ✓ (guided flexibility) |
| **Transferability** | ✓ (structure) | ✗ (data-dependent) | ✓ (structure transfers) |

---

## Theoretical Foundation

### Why This Works

1. **MDNN computes the modelled part**: Physics that is well understood, modeled accurately, discretizable
2. **DDNN computes the unmodelled part**: Unmodeled dynamics, parameter uncertainty, systematic errors
3. **Together**: Model structure and learning flexibility, with neither removed

### Automatic Differentiation Through Both

Since both MDNN (stencils, solvers) and DDNN (neural network) are differentiable:

```
∂loss/∂design = ∂loss/∂(MDNN) + ∂loss/∂(DDNN)
```

No additional code is needed; the chain rule applies.

### Data Efficiency

Pure DDNN needs O(n^d) data to learn in d dimensions (curse of dimensionality).
MDNN provides structure that reduces data requirement to O(m) where m = number of unknowns.

**Example**: 
- Pure DDNN on 2D diffusion: need ~10,000 simulations
- MDNN + DDNN: need ~10 data points (to learn k and one correction term)

---

## Computational Graph Integration

### Each Node Can Be Hybrid

```fortran
! At expansion level
exp_graph = construct_expansion(verts, edges)
! branch(2) = MDNN (PDE structure)
! branch(1) = DDNN (global learned corrections)

! At block level
block_graph = construct_block(...)
! branch(2) = MDNN (one time step, one solver)
! branch(1) = DDNN (learned block-specific params)

! At component level
component_graph = construct_stencil(...)
! branch(2) = MDNN (the stencil itself)
! branch(1) = DDNN (learned stencil weight corrections)
```

### No Separate Systems

One graph, two branches. No parallel infrastructure needed.

---

## Next Steps

1. **Implement hybrid evaluation and adjoint**
   - Add `evaluate_hybrid()` method to graph type
   - Add `adjoint_hybrid()` for backpropagation

2. **Test on concrete problems**
   - Van der Pol with learned damping
   - Diffusion with learned diffusivity
   - Multi-fidelity coarse+fine

3. **Benchmark against baselines**
   - Pure MDNN (fixed parameters)
   - Pure DDNN (no structure)
   - Hybrid (this work)

4. **Extend to hierarchical learning**
   - Learn different parameters at different levels
   - Adaptive refinement using DDNN gradients

5. **Uncertainty quantification**
   - Represent learned parameters as distributions
   - Propagate uncertainty through MDNN + DDNN

---

## Conclusion

The hybrid MDNN/DDNN architecture unifies model-driven and data-driven learning into a single, coherent framework. The fractal graph's binary structure was not accidental: it is the representation of this duality.

**Result**: interpretability and adaptability are both retained; no choice between them is required.

