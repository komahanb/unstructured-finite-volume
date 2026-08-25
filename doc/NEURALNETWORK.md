# Neural Networks in the Differentiable PDE Solver

## Overview

The finite volume framework is fundamentally a **differentiable operator**: it maps design parameters to functional outputs via stencils, time integration, and solvers. A `type :: neural_network` wraps this capability into a learnable, composable abstraction.

The key insight is that traditional neural networks and PDE solvers are mathematically identical:
- **Forward pass**: compute outputs from inputs
- **Backpropagation**: compute gradients via automatic differentiation
- **Training**: optimize parameters via gradient descent

This document describes the neural network type, its hierarchical deployment, and concrete applications.

---

## Analogy: Neural Network ↔ PDE Solver

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
  use gti_expansion, only : expansion, family_holder
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
  type(family_holder), intent(in) :: schemes(:)
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

#### `child_networks`: Spawn networks at child nodes

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

**Efficiency**: Gradient of block B doesn't depend on state of block A—can be parallelized.

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

! Backward: gradients flow independently
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

Networks at multiple levels working together.

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

! Backward: gradients flow from blocks to root
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

## Applications

### 1. Parameter Estimation (Inverse Problem)

**Goal**: Given measurements of solution, find best design.

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
- **Hierarchical Decomposition for PDE Inverse Problems**: [Your work]

