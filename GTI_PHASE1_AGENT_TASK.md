# Agent task: GTI Phase 1 — differentiable form public interface

Repository: `komahanb/unstructured-finite-volume`  
Base: current `master` unless the user gives another branch  
Reference architecture: `GTI.md`  
Task scope: Phase 1 only

## Objective

Install the public GTI differentiable-form skeleton.

This PR must make it possible for a user to define both a residual form `R` and a functional form `F` through one shared abstract interface:

```text
value
partial_action
max_degree
input_signature
output_signature
```

No time integration is required in this PR.

## Hard boundaries

Do not edit:

```text
src/fractal_graph.f90
```

Do not implement:

```text
BDF
DIRK
ABM
adaptive stepping
time graph construction
forward traversal
adjoint traversal
tangent traversal
Newton coupling
linear solver coupling
higher-order chain-rule assembly
FV coupling
mesh changes
```

Do not add hidden semantic state to:

```text
graph
fields
set representations
relation representations
forms
transforms
partition/assembly relation objects
```

## Files to add

Add these modules:

```text
src/gti_value_buffer.f90
src/gti_state_bundle.f90
src/gti_design_bundle.f90
src/gti_evaluation_point.f90
src/gti_form_interface.f90
```

Add one test directory, name may follow existing repo style:

```text
test/gti-form-interface/test.f90
```

Update build/test manifests only as needed for the new modules and test.

## Required module responsibilities

### `gti_value_buffer.f90`

Provide a small real-valued buffer for Phase 1 tests.

Minimum shape:

```fortran
type :: gti_value_buffer
   integer :: value_kind = GTI_VALUE_REAL
   integer :: nentries = 0
   integer :: ncomp = 1
   real(dp), allocatable :: rvals(:)
contains
   procedure :: clear
   procedure :: set_real
   procedure :: get_real
end type
```

It does not need to support integer, complex, logical, or character values in this PR.

### `gti_state_bundle.f90`

Provide a state bundle with polymorphic field slots.

Do not declare `type(graph_field)`, because `graph_field` is abstract.

Use:

```fortran
type :: gti_field_slot
   class(graph_field), allocatable :: value
end type
```

Minimum state bundle:

```fortran
type :: gti_state_bundle
   type(gti_field_slot), allocatable :: component(:)
contains
   procedure :: differential_degree
   procedure :: has_component
end type
```

Fortran storage convention:

```text
component(1) = q
component(2) = qdot
component(3) = qddot
```

### `gti_design_bundle.f90`

Provide a design bundle using the same `gti_field_slot` concept or a compatible local slot.

Minimum:

```fortran
type :: gti_design_bundle
   type(gti_field_slot), allocatable :: component(:)
contains
   procedure :: size
   procedure :: has_entries
end type
```

### `gti_evaluation_point.f90`

Provide the local point passed to forms.

Minimum:

```fortran
type :: gti_evaluation_point
   type(gti_state_bundle)  :: state
   type(gti_design_bundle) :: design
   real(dp)                :: time = 0.0_dp
   integer                 :: window_id = 0
   integer                 :: step_id = 0
   integer                 :: stage_id = 0
end type
```

Do not add mesh or FV dependencies in this PR.

### `gti_form_interface.f90`

Provide the public abstract form interface.

Minimum:

```fortran
type, abstract :: gti_differentiable_form
contains
   procedure(form_name_interface)            , deferred :: name
   procedure(form_input_signature_interface) , deferred :: input_signature
   procedure(form_output_signature_interface), deferred :: output_signature
   procedure(form_max_degree_interface)      , deferred :: max_degree
   procedure(form_value_interface)           , deferred :: value
   procedure(form_partial_action_interface)  , deferred :: partial_action
end type
```

Also define:

```fortran
type :: gti_partial_request
   integer :: order = 0
   integer, allocatable :: argument_kind(:)
   integer, allocatable :: state_component(:)
end type
```

and:

```fortran
type :: gti_direction_bundle
   integer :: argument_kind = GTI_ARG_STATE
   integer :: state_component = GTI_STATE_Q
   type(gti_value_buffer) :: values
end type
```

Public constants:

```fortran
GTI_ARG_STATE
GTI_ARG_DESIGN
GTI_ARG_TIME
GTI_ARG_GEOM

GTI_STATE_Q
GTI_STATE_QDOT
GTI_STATE_QDDOT

GTI_VALUE_REAL
```

## Required tests

Create a toy residual form and toy functional form that both extend `gti_differentiable_form`.

Example mathematical content is allowed to be simple:

```text
R(q, xi) = q^2 + xi
F(q, xi) = 0.5*q^2 + xi
```

Test:

```text
R and F share the same abstract interface
value works for R
value works for F
partial_action order 0 works
partial_action order 1 works
partial_action order 2 works
max_degree rejects unsupported degree
state bundle reports components correctly
design bundle reports entries correctly
```

No numerical solver test is required.

## Acceptance criteria

The PR is acceptable only if:

```text
new modules compile
new test passes
existing tests still pass
src/fractal_graph.f90 unchanged
no scheme modules added
no traversal modules added
no solver adapter added
no mesh/FV dependency added
no hidden parent/ambient state introduced
```

## Commit/PR title suggestion

```text
Add GTI differentiable form interface skeleton
```

## Summary sentence for PR description

```text
This PR installs the public GTI form contract shared by residual and functional objects, including value and partial-action calls up to toy second order, without adding time schemes, traversals, solvers, or kernel changes.
```
