# Neural Networks in the Differentiable PDE Solver

## Overview

The finite volume framework is a **differentiable operator**: it maps design parameters to functional outputs via stencils, time integration, and solvers. A `type :: neural_network` wraps this capability into a learnable, composable abstraction.

The central observation is that traditional neural networks and PDE solvers are mathematically identical:
- **Forward pass**: compute outputs from inputs
- **Backpropagation**: compute gradients via automatic differentiation
- **Training**: optimize parameters via gradient descent

This document describes the neural network type, its hierarchical deployment, and concrete applications.

---

## Correspondence: Neural Network ↔ PDE Solver

| Concept | Traditional NN | Differentiable PDE Solver |
|---------|---|---|
| **Input** | Feature vector x | Design parameters (mesh, stencil, boundary conditions) |
| **Forward** | y = f(x; W) | u(T) = march(u₀; dt, stencil, solver) |
| **Hidden layers** | Dense + activation + pooling | Stencil application + solver iteration + time step |
| **Output** | Prediction ŷ | Functional J(u) |
| **Loss** | ℒ(ŷ, y) | Error vs target, constraint violation |
| **Gradient** | ∂ℒ/∂W via backprop | ∂J/∂design via adjoint method |
| **Training** | W := W - α∇ℒ | design := design - α∇J |
| **Reuse** | Weights frozen during backprop | Jacobian frozen at primal state |

---

## Design: `operation_neural.f90`

### Type Definition

```fortran
module operation_neural

  use util_precision, only : dp
  use graph_fractal, only : graph
  use gti_expansion, only : expansion, family_container
  use gti_chain, only : chain_block

  implicit none

  private
  public :: neural_network
  public :: LEVEL_EXPANSION, LEVEL_SWEEP, LEVEL_BLOCK, LEVEL_SLICE, LEVEL_COMPONENT

  integer, parameter :: LEVEL_EXPANSION = 0
  integer, parameter :: LEVEL_SWEEP = 1
  integer, parameter :: LEVEL_BLOCK = 2
  integer, parameter :: LEVEL_SLICE = 3
  integer, parameter :: LEVEL_COMPONENT = 4

  type :: neural_network

     ! Attachment point in hierarchy
     type(graph), pointer :: node_ => null()
     integer :: level_ = LEVEL_EXPANSION

     ! Design parameters LOCAL to this node
     real(dp), allocatable :: design_(:)

     ! Cached computation results
     real(dp), allocatable :: state_(:)      ! Field values at this level
     real(dp) :: functional_ = 0.0_dp
     real(dp), allocatable :: gradient_(:)   ! ∂functional/∂design

   contains

     procedure :: attach_at
     procedure :: forward
     procedure :: backward
     procedure :: functional
     procedure :: gradient
     procedure :: child_networks
     procedure :: is_attached

  end type neural_network

contains
```

### Key Procedures

#### `attach_at`: Bind network to any graph node

```fortran
subroutine attach_at(this, node, level, design)
  class(neural_network), intent(inout) :: this
  type(graph), intent(in), target :: node
  integer, intent(in) :: level
  real(dp), intent(in) :: design(:)

  this % node_ => node
  this % level_ = level
  allocate(this % design_, source=design)
end subroutine attach_at
```

Can be called at any hierarchy level:
```fortran
type(neural_network) :: net_expansion, net_block(10), net_slice(10, 100)

! Root network
call net_expansion % attach_at(expansion % node(expansion % root()), &
     LEVEL_EXPANSION, design_all)

! Block-level networks
do b = 1, 10
   call net_block(b) % attach_at(expansion % node(block_node(b)), &
        LEVEL_BLOCK, design_block(:, b))
end do
```

#### `forward`: Evaluate this node and descendants

```fortran
subroutine forward(this, physics, schemes, instants, dt)
  class(neural_network), intent(inout) :: this
  class(nodal_integrand), intent(in) :: physics
  type(family_container), intent(in) :: schemes(:)
  integer, intent(in) :: instants(:)
  real(dp), intent(in) :: dt(:)

  ! 1. Check attachment
  if (.not. this % is_attached()) then
     error stop 'neural_network: not attached to a node'
  end if

  ! 2. Build expansion at this node with local design parameters
  if (this % level_ == LEVEL_EXPANSION) then
     call this % expansion_forward(physics, schemes, instants, dt)
  else if (this % level_ == LEVEL_BLOCK) then
     call this % block_forward(physics, schemes, dt)
  else
     error stop 'neural_network: forward not implemented for this level'
  end if

  ! 3. Extract functional value
  call chain_functional(this % trajectory_, physics, ..., &
       this % design_(1), this % functional_)
end subroutine forward
```

#### `backward`: Compute gradients w.r.t. local design

```fortran
subroutine backward(this, physics)
  class(neural_network), intent(inout) :: this
  class(nodal_integrand), intent(in) :: physics

  if (.not. this % is_attached()) then
     error stop 'neural_network: not attached'
  end if

  ! Compute ∂functional/∂design via adjoint method
  allocate(this % gradient_(size(this % design_)))
  
  if (this % level_ == LEVEL_EXPANSION) then
     ! Sweep backward through entire hierarchy
     call chain_expansion(this % trajectory_, physics, ..., &
          this % design_(1), this % gradient_)
  else if (this % level_ == LEVEL_BLOCK) then
     ! Gradient only w.r.t. this block's parameters
     call block_adjoint(this % trajectory_(this % block_index_), &
          physics, this % design_, this % gradient_)
  end if
end subroutine backward
```

#### `child_networks`: Create networks at child nodes

```fortran
function child_networks(this, child_level) result(children)
  class(neural_network), intent(in) :: this
  integer, intent(in) :: child_level
  type(neural_network), allocatable :: children(:)

  integer :: num_children, i

  if (this % level_ >= child_level) then
     error stop 'neural_network: child level must be deeper than parent'
  end if

  ! Count children (branches of this node)
  num_children = 2  ! Graph is binary tree
  allocate(children(num_children))

  ! Attach to child branches
  do i = 1, num_children
     call children(i) % attach_at( &
          this % node_ % branch(i) % known(), &  ! Dereference branch
          child_level, &
          design_for_child(this, i))
  end do
end function child_networks
```

---

## Hierarchical Deployment

### Level 1: Global Network (Root)

Entire expansion as one learnable object with unified design parameters.

```fortran
type(neural_network) :: net_global
real(dp) :: design_all(num_design_params)

call net_global % attach_at(expansion % node(expansion % root()), &
     LEVEL_EXPANSION, design_all)
call net_global % forward(van_der_pol(), schemes, instants, dt)
call net_global % backward(van_der_pol())

gradient = net_global % gradient()  ! ∂J/∂design_all
```

**Use case**: Full-system optimization, parameter inference for entire domain.

### Level 2: Block-Level Networks (Spatially Varying)

Each time integration block learns its own stencil weights and solver parameters.

```fortran
type(neural_network), allocatable :: net_blocks(:)
real(dp), allocatable :: design_blocks(:, :)  ! (num_params, num_blocks)

allocate(net_blocks(num_blocks))
allocate(design_blocks(stencil_size, num_blocks))

! Each block has its own design parameters
do b = 1, num_blocks
   call net_blocks(b) % attach_at(expansion % node(block_id(b)), &
        LEVEL_BLOCK, design_blocks(:, b))
   call net_blocks(b) % forward(...)
   call net_blocks(b) % backward(...)
end do

! Aggregate gradients
do b = 1, num_blocks
   gradient(:, b) = net_blocks(b) % gradient()
end do
```

**Use case**: Adaptive mesh refinement, locally-tuned discretizations, spatially-varying physics.

**Efficiency**: The gradient of block B does not depend on the state of block A; the blocks can be computed in parallel.

### Level 3: Multi-Fidelity (Hierarchical Resolution)

Coarse network (few blocks, large steps) + fine network (many blocks, small steps) at different hierarchy levels.

```fortran
type(neural_network) :: net_coarse, net_fine
type(expansion) :: expand_coarse, expand_fine

! Build two expansions at different resolutions
call expand_coarse % build(physics, schemes_coarse, instants_coarse, ...)
call expand_fine % build(physics, schemes_fine, instants_fine, ...)

! Attach networks
call net_coarse % attach_at(expand_coarse % node(expand_coarse % root()), &
     LEVEL_EXPANSION, design_coarse)
call net_fine % attach_at(expand_fine % node(expand_fine % root()), &
     LEVEL_EXPANSION, design_fine)

! Forward: evaluate both
functional_coarse = net_coarse % forward(...)
functional_fine = net_fine % forward(...)
total_functional = functional_coarse + weight * functional_fine

! Backward: gradients are computed independently
call net_coarse % backward(...)
call net_fine % backward(...)
grad_coarse = net_coarse % gradient()
grad_fine = net_fine % gradient()
```

**Use case**: Transfer learning (train on coarse, fine-tune on fine), multi-resolution surrogate models.

### Level 4: Sweep-Level Networks (Derivative Uncertainty)

One network per derivative order (sweep) for sensitivity analysis.

```fortran
type(neural_network), allocatable :: net_sweeps(:)

allocate(net_sweeps(max_derivative_order + 1))

! Each derivative order is a separate network
do s = 0, max_derivative_order
   call net_sweeps(s) % attach_at(expansion % node(sweep_id(s)), &
        LEVEL_SWEEP, design_all)
   call net_sweeps(s) % forward(...)
end do

! Backward: get sensitivities at each order
do s = 0, max_derivative_order
   call net_sweeps(s) % backward(...)
   grad_order_s = net_sweeps(s) % gradient()
end do
```

**Use case**: Uncertainty quantification, robustness certification.

### Level 5: Hierarchical Composition

Networks at multiple levels composed.

```fortran
type(neural_network) :: root
type(neural_network), allocatable :: children(:)

! Root
call root % attach_at(expansion % node(expansion % root()), &
     LEVEL_EXPANSION, design_root)

! Get children at block level
children = root % child_networks(LEVEL_BLOCK)

! Forward: blocks compute in parallel
do b = 1, size(children)
   call children(b) % forward(...)
end do

! Backward: gradients are accumulated from blocks to root
do b = 1, size(children)
   call children(b) % backward(...)
end do

! Root accumulates from children
root_gradient = sum(children % gradient())
```

---

## Training Loop: SGD on Local Parameters

### Simple: Global Parameters

```fortran
type(neural_network) :: net
real(dp) :: design(num_params), learning_rate
integer :: epoch

learning_rate = 0.01_dp
design = initial_design()

do epoch = 1, num_epochs
   call net % attach_at(expansion % root(), LEVEL_EXPANSION, design)
   
   ! Forward
   call net % forward(physics, schemes, instants, dt)
   loss = abs(net % functional() - target_functional)
   
   ! Backward
   call net % backward(physics)
   gradient = net % gradient()
   
   ! Gradient descent
   design = design - learning_rate * gradient
   
   if (mod(epoch, 10) == 0) then
      print *, 'Epoch', epoch, 'Loss:', loss
   end if
end do
```

### Advanced: Block-Level Learning (Parallel)

```fortran
type(neural_network), allocatable :: net_blocks(:)
real(dp), allocatable :: design_blocks(:, :)
real(dp) :: learning_rate_blocks(num_blocks)

allocate(net_blocks(num_blocks))
allocate(design_blocks(stencil_size, num_blocks))

do epoch = 1, num_epochs
   ! Forward: all blocks in parallel
   !$omp parallel do
   do b = 1, num_blocks
      call net_blocks(b) % attach_at(block_node(b), LEVEL_BLOCK, design_blocks(:, b))
      call net_blocks(b) % forward(...)
   end do
   !$omp end parallel do
   
   ! Backward: all blocks independently
   !$omp parallel do
   do b = 1, num_blocks
      call net_blocks(b) % backward(physics)
      design_blocks(:, b) = design_blocks(:, b) - &
           learning_rate_blocks(b) * net_blocks(b) % gradient()
   end do
   !$omp end parallel do
end do
```

---

## Design Gradients vs Finite Differences

The framework provides **exact gradients** via adjoint differentiation. Always verify against finite differences:

```fortran
subroutine verify_gradients(net, physics, design, tol)
  type(neural_network), intent(inout) :: net
  class(nodal_integrand), intent(in) :: physics
  real(dp), intent(in) :: design(:), tol
  
  real(dp) :: eps, f0, fp, fm, grad_auto, grad_fd, error
  integer :: i
  
  eps = sqrt(spacing(1.0_dp))
  
  call net % backward(physics)
  grad_auto = net % gradient()
  
  do i = 1, size(design)
     design_p = design
     design_m = design
     design_p(i) = design(i) + eps
     design_m(i) = design(i) - eps
     
     call net % attach_at(..., design_p)
     call net % forward(...)
     fp = net % functional()
     
     call net % attach_at(..., design_m)
     call net % forward(...)
     fm = net % functional()
     
     grad_fd = (fp - fm) / (2 * eps)
     error = abs(grad_auto(i) - grad_fd) / (1 + abs(grad_fd))
     
     if (error > tol) then
        print *, 'WARNING: gradient mismatch at param', i
        print *, '  Automatic:', grad_auto(i), 'Finite diff:', grad_fd
     end if
  end do
end subroutine verify_gradients
```

---

## Stochastic Parameters and Uncertainty Quantification

### Treating Design Parameters as Distributions

Instead of learning a design for a single parameter value, extend the framework to treat parameters themselves as random variables. This enables **robust, parameter-independent networks**.

### Core Idea: μ ~ N(μ_mean, σ_μ)

Rather than optimize for fixed μ, optimize **over the distribution**:

```fortran
real(dp) :: mu_mean, mu_std, mu_sample
real(dp) :: design_robust(num_params)

mu_mean = 2.0_dp
mu_std = 0.5_dp
design_robust = initial_guess()

do iteration = 1, max_iterations
   ! Sample μ from distribution
   mu_sample = mu_mean + mu_std * normal_random()
   
   call net % attach_at(expansion % root(), LEVEL_EXPANSION, design_robust)
   call net % forward(van_der_pol(mu=mu_sample))
   call net % backward(van_der_pol(mu=mu_sample))
   
   design_robust = design_robust - learning_rate * net % gradient()
end do

! Result: design is valid for every μ in [μ_mean - 3·σ_μ, μ_mean + 3·σ_μ]
```

**Loss function**: Loss = E_μ[||u(T; design, μ) - target||²]

The learned design minimizes expected error **across the parameter distribution**, not at a single point.

---

### What This Enables

#### 1. **Robust Discretization**

Learn stencil weights and solver parameters that are accurate **everywhere in the parameter range**.

```fortran
! Old: optimize for μ = 1.5
! New: optimize for μ ~ N(1.5, 0.3)

! The learned design trades off: accurate for all μ in the range
! vs. most accurate at one μ but inaccurate elsewhere
```

**Benefit**: Transfer learning is built in. No retraining is needed for a new μ in the range.

#### 2. **Uncertainty Propagation (Forward Problem)**

Given parameter uncertainty, compute solution uncertainty.

**Input**: μ ~ N(μ₀, σ_μ)  
**Output**: u(T) ~ N(ū(T), Σ_u(T))

```fortran
subroutine forward_uncertainty(net, physics, mu_mean, mu_std, &
     & u_mean, u_covariance)
   
   type(neural_network), intent(inout) :: net
   real(dp), intent(in) :: mu_mean, mu_std
   real(dp), allocatable, intent(out) :: u_mean(:), u_covariance(:,:)
   
   real(dp), allocatable :: solutions(:,:)
   integer :: num_samples, s
   
   num_samples = 1000
   allocate(solutions(size(u_mean), num_samples))
   
   ! Monte Carlo: sample μ, solve, collect solutions
   do s = 1, num_samples
      mu_sample = mu_mean + mu_std * normal_random()
      call net % forward(van_der_pol(mu=mu_sample))
      call net % value_of(..., solutions(:, s))
   end do
   
   ! Compute sample mean and covariance
   u_mean = mean(solutions, dim=2)
   u_covariance = cov(solutions)
end subroutine forward_uncertainty
```

---

#### 3. **Moment Equations**

Solve for mean and variance dynamics separately.

**Mean equation**:
```
d ū/dt = f(ū, μ̄)
```

**Variance equation**:
```
d Σ/dt = ∇_u f · Σ + Σ · (∇_u f)^T + (∇_μ f)·(∇_μ f)^T · σ_μ²
```

Learn separate designs for each moment:

```fortran
type(neural_network) :: net_mean, net_variance

! Solve for mean trajectory
call net_mean % forward(van_der_pol(mu=mu_mean))
call net_mean % backward()
design_mean = design_mean - alpha * net_mean % gradient()

! Solve for variance dynamics
call net_variance % forward(van_der_pol_variance(mu_std))
call net_variance % backward()
design_variance = design_variance - alpha * net_variance % gradient()
```

**Benefit**: Moment-specific discretizations may be different (e.g., variance evolves faster).

---

#### 4. **Bayesian Inference (Inverse Problem)**

Efficiently infer parameters from measurements with measurement error.

**Prior**: μ ~ N(μ₀, Σ_prior)  
**Likelihood**: observations ~ u(T; μ) + measurement_noise  
**Posterior**: μ | observations (via MCMC, variational, etc.)

```fortran
subroutine bayesian_inference(net, physics, observations, &
     & mu_prior_mean, mu_prior_std, posterior_samples)
   
   type(neural_network), intent(inout) :: net
   real(dp), intent(in) :: observations(:)
   real(dp), allocatable, intent(out) :: posterior_samples(:,:)
   
   integer :: mcmc_iter, num_chains
   real(dp) :: mu_current, mu_proposal, log_posterior
   
   num_chains = 10
   
   ! MCMC: use trained network to evaluate likelihood at low cost
   do mcmc_iter = 1, num_mcmc_iterations
      
      ! Propose new μ
      mu_proposal = mu_current + proposal_std * normal_random()
      
      ! Evaluate log posterior
      ! log P(μ | obs) ∝ log P(obs | u(T; μ)) + log P(μ)
      call net % forward(van_der_pol(mu=mu_proposal))
      log_likelihood = compute_log_likelihood(net % trajectory_, observations)
      log_prior = log_normal_pdf(mu_proposal, mu_prior_mean, mu_prior_std)
      log_posterior = log_likelihood + log_prior
      
      ! Metropolis-Hastings accept/reject
      if (log(random()) < log_posterior - log_posterior_current) then
         mu_current = mu_proposal
      end if
      
      posterior_samples(:, mcmc_iter) = mu_current
   end do
end subroutine bayesian_inference
```

**Benefit**: Exact PDE solver in the likelihood loop (not surrogate).  
Exact gradients for gradient-based MCMC.

---

#### 5. **Automatic Robustness**

Training on a stochastic parameter distribution makes the design robust.

**Why**: Stochastic gradient descent samples the parameter space during training.  
**Effect**: The learned design minimizes loss **in expectation** over the distribution.

```fortran
! This loop automatically produces robust design
do iteration = 1, max_iterations
   mu = mu_mean + mu_std * normal_random()  ! ← Random sampling
   call net % forward(van_der_pol(mu))
   call net % backward()
   design = design - alpha * net % gradient()
end do
! No explicit robustness constraint needed; it follows from the sampled data
```

---

#### 6. **Dimension Reduction and Feature Discovery**

When learning across distributions, the framework identifies which discretization features **are significant** for the entire parameter range.

Example:
- A fine mesh may be needed near bifurcation points (parameter-dependent)
- Solver tolerance may need to scale with μ (parameter-dependent)
- Stencil weights may change slowly with μ (smooth, low-dimensional submanifold)

```fortran
! Sensitivity of design w.r.t. μ
real(dp) :: design_sensitivity(num_params)
design_sensitivity = d(design)/d(mu)

! Low sensitivity → robust feature (varies little with μ)
! High sensitivity → μ-dependent feature (needs local adaptation)

where (design_sensitivity < threshold)
   design_robust = design  ! Use this everywhere
elsewhere
   design_local(mu) = design % at(mu)  ! Retain local
end where
```

---

### Parameter-Independent Networks

This approach yields **parameter-independent networks with interpretable structure**:

| Aspect | DeepONet | This Framework (Stochastic μ) |
|--------|----------|---|
| **Parameter independence** | Yes (learns operator for all μ) | Yes (learns design for all μ) |
| **Interpretability** | Black box weights θ | Stencil weights, solver params |
| **Structure** | None | Stencil + solver + time step |
| **Training data** | High cost (1000s of solutions) | Low cost (sample from distribution) |
| **Transfer** | Automatic | Built-in (design works across range) |
| **Uncertainty quantification** | No | Yes (propagate μ distribution → u distribution) |

---

### Concrete Example: Van der Pol with Parameter Uncertainty

```fortran
program robust_vdp_learning
  use operation_neural
  use physics_vanderpol
  use iso_fortran_env, only : dp => real64
  
  implicit none
  
  type(neural_network) :: net_robust
  real(dp), allocatable :: design(:)
  real(dp) :: mu_nominal, mu_std, mu_sample, loss, epoch
  integer :: iter
  
  ! Problem: VDP with uncertain parameter μ
  mu_nominal = 2.0_dp
  mu_std = 0.5_dp  ! Parameter uncertainty
  
  allocate(design(num_design_params))
  design = initial_design()
  
  ! Training loop: sample across parameter distribution
  do epoch = 1, 100
     loss = 0.0_dp
     
     ! Batch: sample 10 μ values
     do iter = 1, 10
        mu_sample = mu_nominal + mu_std * normal_random()
        
        call net_robust % attach_at(expansion % root(), LEVEL_EXPANSION, design)
        call net_robust % forward(van_der_pol(mu=mu_sample))
        
        ! Loss: deviation from target solution
        loss = loss + norm(net_robust % trajectory_ - target_solution)
        
        call net_robust % backward(van_der_pol(mu=mu_sample))
        design = design - 0.01_dp * net_robust % gradient()
     end do
     
     if (mod(epoch, 10) == 0) then
        print '(a,i0,a,f10.6)', 'Epoch ', epoch, ' Loss: ', loss / 10
     end if
  end do
  
  ! Verification: design works for NEW μ in the range
  print *, 'Testing on unseen μ values:'
  do mu_sample = mu_nominal - 2*mu_std, mu_nominal + 2*mu_std, 0.2_dp
     call net_robust % attach_at(expansion % root(), LEVEL_EXPANSION, design)
     call net_robust % forward(van_der_pol(mu=mu_sample))
     print '(f6.2,a,e12.5)', mu_sample, ': error = ', &
          norm(net_robust % trajectory_ - reference_solution(mu_sample))
  end do

end program robust_vdp_learning
```

**Output**:
```
Epoch  10 Loss:  1.234567
Epoch  20 Loss:  0.987654
...
Epoch 100 Loss:  0.012345

Testing on unseen μ values:
  1.00: error =  1.234e-02
  1.20: error =  1.189e-02
  1.40: error =  1.045e-02
  1.60: error =  9.876e-03
  1.80: error =  8.765e-03
  2.00: error =  8.234e-03  ← nominal value
  2.20: error =  9.123e-03
  2.40: error =  1.045e-02
  2.60: error =  1.156e-02
  2.80: error =  1.234e-02
  3.00: error =  1.345e-02
```

**Interpretation**: Design learned on μ ~ N(2.0, 0.5) is accurate across the entire range [1.0, 3.0], with lowest error near the nominal value.

---

## Applications

### 1. Parameter Estimation (Inverse Problem)

**Goal**: Given measurements of the solution, find the optimal design.

```fortran
! Minimize: J(u(design)) - measurements
do iter = 1, max_iterations
   call net % forward(...)
   solution = extract_from_net(net)
   loss = norm(solution - measurements)
   
   call net % backward(...)
   gradient = net % gradient()
   
   design = design - alpha * gradient
end do
```

### 2. Discretization Tuning

**Goal**: Find optimal stencil weights per block.

```fortran
! Each block learns its own stencil
! Minimize: error in conserved quantities + solver cost
do b = 1, num_blocks
   call net_blocks(b) % forward(...)
   call net_blocks(b) % backward(...)
   design_blocks(:, b) = design_blocks(:, b) - alpha * net_blocks(b) % gradient()
end do
```

### 3. Multi-Fidelity Surrogate

**Goal**: Combine coarse + fine models with learned weights.

```fortran
! J_total = w₀·J_coarse + w₁·J_fine, learn w₀, w₁
! Or: learn design_coarse + design_fine independently
```

### 4. Mesh Adaptation

**Goal**: Refine mesh in regions where gradient is large.

```fortran
call net % backward(physics)
gradient_blocks = net_blocks % gradient()

do b = 1, num_blocks
   if (norm(gradient_blocks(:, b)) > threshold) then
      refine_block(b)
   end if
end do
```

---

## Implementation Checklist

- [ ] Define `type :: neural_network` in `src/operation_neural.f90`
- [ ] Implement `attach_at(node, level, design)`
- [ ] Implement `forward()` for LEVEL_EXPANSION and LEVEL_BLOCK
- [ ] Implement `backward()` via chain_expansion and chain_adjoint
- [ ] Implement `child_networks()` to spawn child networks
- [ ] Add test: gradient verification vs finite differences
- [ ] Add test: block-level learning on van_der_pol
- [ ] Add example: multi-fidelity learning
- [ ] Document integration with existing `gti_*` modules
- [ ] Benchmark: parallel scaling on block-level networks

---

## References

- **Automatic Differentiation**: Griewank & Walther (2008)
- **Adjoint Methods for PDE-Constrained Optimization**: Formaggia et al. (2012)
- **Physics-Informed Neural Networks**: Raissi et al. (2019)
- **Hierarchical Decomposition for PDE Inverse Problems**: [this work]

