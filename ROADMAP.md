# roadmap: euler -> rans -> navier-stokes

where the solver is headed: from the scalar advection-diffusion it does today to
compressible euler, then rans, then full navier-stokes. the focus of this doc is
INFRASTRUCTURE READINESS - what the existing seams give us for free, and the
concrete gaps each milestone must close. (companion to PLAN.md, the
general-purpose refactor, now largely done.)

## what we already have (the reusable spine)

- law-agnostic flux/source/objective seam (`interface_physics`): pointwise
  F(q,grad q), S, f as `type(scalar)`, with BLOCK state jacobians
  `dflux_dq (3,nv,nv)` and `dflux_dgradq (3,nv,3,nv)`. the (nv,nv) blocks are the
  systems hook - designed for coupled variables, used only diagonally so far.
- `point_state` carries q(nv) + grad q(3,nv); `type(scalar)` (scalar.fpp) so a
  complex-step jacobian drops in.
- the fv assembler integrates ∮F.n over faces + ∫S over the volume; `flux_keff`
  (diffusion, from dF/dgradq) and `flux_vn` (advection, from dF/dq); central or
  upwind faces (`convection_scheme`).
- the IMPLICIT stack - exactly what compressible cfd wants: matrix-free newton
  (`class_newton_solver`) + bdf 1-6 (`class_bdf`) + krylov (`class_conjugate_gradient`,
  `class_gmres`, `class_normal_cg`) + preconditioners (`class_algebraic_multigrid`,
  block-jacobi) + steady & transient adjoint (adjoint trajectory in `class_bdf`).
- assembled csr (`get_operator_csr`) + GMRES for nonsymmetric operators; coarray
  domain decomposition (`class_distributed_cg`) + RCB partitioner (`graph % partition_rcb`).
- mesh: cell-centred fv, centroids / volumes / face areas+normals, 2d & 3d, named
  boundary groups.
- clean interface/class architecture: every concrete class extends `interface_*.f90`;
  stateless solvers (`linear_solver % solve(system, x, mode)`); graph owns partitioning;
  flux projections on the flux base; convergence monitor on the solver base.

## cross-cutting infrastructure gaps (needed by ALL three milestones)

do these first - they unlock everything downstream:

1. **multi-field BLOCK assembly.** euler is nv = 4 (2d) / 5 (3d) coupled conserved
   variables. the flux blocks (3,nv,nv) are ready; the assembler writes only the
   diagonal. need: block `get_operator_csr` (nv x nv dense block per graph edge ->
   block-csr, the deferred BCSRMat path), block matvec, block-jacobi / block-ILU
   smoother. `graph % dof` already interleaves variables.
   files: `class_csr` (block_size), `class_assembler` (ivar x jvar block loop),
   `class_amg` / smoothers (nodal/block).

2. **state-dependent (nonlinear) assembly + newton.** today the jacobian is built
   at st%q = 0 (the operator is linear). euler/ns are nonlinear:
   R(q) = ∮F(q).n - ∫S, with dR/dq evaluated at the CURRENT q. need: the assembler
   to reconstruct the real face state and evaluate F + dF there, and feed R and
   dR/dq to the existing newton. `flux % value` already exists; this is the key
   assembler change.
   files: `class_assembler` (residual + state-dependent jacobian), `nonlinear_marching`.

3. **numerical-flux (riemann) seam.** hyperbolic systems need a two-state face
   flux F^(q_L, q_R, n), not F(q_face).n. new abstract operator: rusanov /
   local-lax-friedrichs (start) -> hll -> hllc -> roe, with left/right jacobians
   for the implicit path. `face_state` already gives q_L (owner) and q_R (neighbour).
   files: `interface_physics` (numerical_flux type), `class_assembler` (face loop
   calls F^), new `class_riemann_*`.

4. **gradient reconstruction + LIMITERS (2nd order).** 1st-order upwind is stable
   but diffusive; 2nd order needs cell-gradient reconstruction (green-gauss /
   least-squares) + a slope limiter (barth-jespersen, venkatakrishnan) for
   monotonicity at shocks. generalises the normal-gradient face reconstruction
   already used for diffusion.
   files: `class_fvm_field` (limited reconstruction), `class_mesh` (lsq stencils).

5. **thermodynamics / EOS.** ideal gas p = (g-1)(rE - 1/2 r|u|^2); conserved <->
   primitive conversion; sutherland mu(T) for the viscous milestones.
   files: new `class_thermodynamics`.

6. **system boundary conditions.** scalar robin won't do: characteristic /
   farfield (riemann invariants), slip wall (euler), no-slip wall (ns),
   sub/supersonic in/outflow - boundaries set the boundary numerical flux.
   files: `class_boundary_condition` (system bc kinds) or a new bc layer.

## milestone 1 - euler (inviscid compressible)

dq/dt + div F_inv(q) = 0, q = (r, ru, rv, rw, rE), F_inv carrying pressure (EOS).
the first nonlinear hyperbolic system. needs gaps 1, 2, 3, 5, 6 (+ 4 for 2nd
order). time: explicit ssp-rk (cfl-limited) to start, or the implicit newton+gmres
stack (well-suited - reuse).
verify: 1d sod shock tube (exact riemann), 2d isentropic vortex (smooth -> order
study), supersonic wedge (oblique-shock angle), naca0012 inviscid (cp).

## milestone 2 - rans (reynolds-averaged navier-stokes)

rans = mean-flow navier-stokes (euler + VISCOUS flux) + a turbulence model, so
this milestone lands the viscous terms first, then turbulence transport:
- **viscous flux** F_visc(q, grad q): stress tau = mu(grad u + grad u^T - 2/3 div u I)
  and heat flux -k grad T. gradient-based -> reuses the dF/dgradq path (diffusion
  generalised to a tensor); needs face gradients (gap 4).
- **turbulence transport**: spalart-allmaras (1 equation, robust) first, then
  k-omega sst (2 equations). stiff production/destruction SOURCES (the source seam)
  + eddy viscosity mu_t feeding the viscous flux. each model adds fields -> rides
  on block assembly (gap 1).
files: new `class_viscous_flux`, `class_turbulence_sa` (source + mu_t), no-slip /
wall-function bc, `class_thermodynamics` (mu(T)).
verify: flat plate (blasius / skin-friction cf), backward-facing step
(reattachment length), naca0012 (cl/cd vs experiment).

## milestone 3 - full navier-stokes (the capstone)

the complete compressible ns, integrated and validated: unsteady (urans/des/
dns-ready), all bc types, 2nd-order + limiters throughout, robust newton+gmres with
block-amg / schur preconditioning across the coarray decomposition, and adjoint for
shape/design (the steady+transient adjoint exists - extend to the system + node/
shape sensitivities). more HARDENING + SCALING than new physics.
verify: manufactured-solution order study for the full ns system; unsteady cylinder
vortex shedding (strouhal number); 3d wing; adjoint dJ/dshape vs finite difference.

## suggested order of work

infra gaps 1 (block) + 2 (newton / state-dependent) + 5 (eos) are the foundation -
prove them first on a SCALAR nonlinear test, burgers (F = 1/2 u^2, dF/dq = u), so
the nonlinear machinery is solid before euler's 4-5 fields. then 3 (riemann) + 6
(bcs) -> euler 1st-order (sod). then 4 (limiters) -> euler 2nd-order. then viscous
flux + spalart-allmaras -> rans. then harden + scale -> full ns. the adjoint,
coarray DD, amg, gmres and bdf stacks already exist and ride along throughout.

## tracked deferrals — each with its retirement condition

deferred is acceptable; deferred-but-untracked is not. every standing
deferral in the source lives here as the index (the source comments
stay); an entry retires when its condition is met, not before.

- **the product trio** (get_jacobian_vector_product, the forward
  routing wrapper; add_jacobian_vector_product_transpose; and
  add_design_residual_transpose_product) consolidates into the one
  product — retires when the steady and transient adjoints migrate onto
  the digraph's accumulate_adjoint.
- **cg's linearized march override** (the newton/bdf inner path holding
  lin_coeff and an external right-hand side on the solver) — retires at
  the linearization commit, when that frozen-operator state moves onto
  the system.
- **bdf's hand-rolled march_backwards** — retires when the transient
  adjoint tenants the digraph's reverse accumulation; march_backwards
  becomes an alias, then dissolves.
- **the genuine non-symmetric transpose inside class_assembler** —
  gated on the reserved system-implementation decision; until then
  transpose seats refuse on undeclared non-symmetric instances, and
  normal_cg on such operators runs only where a genuine transpose
  exists (the csr-backed test fixture).
- **direction-threading through the shared linear iteration** — a
  REVERSE march through converge is refused on anything but a
  declared-symmetric operator (the loop is direction-blind); retires
  when the residual query and the sweep carry the mode.
