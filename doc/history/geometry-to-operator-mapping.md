# geometry to operator mapping

rung 1 of the migration. how every geometry number the old assembler
reads becomes a constructor argument of a tower operator. nothing here
is code; this is the dictionary rung 2 builds against.

## the law

geometry enters at construction, once. an operator receives
coefficients, spacings, measures, boundary values — plain per-edge and
per-cell numbers — and never reads a mesh. the mesh answers typed
questions (`mesh % cell_volume()`, `mesh % face_area()`); no string
names any internal thing. a tag string enters exactly once, at setup,
when a boundary condition resolves `tagged_edges(tag)` to a member
set; nothing downstream holds the string.

## the seats

the differential operator carries four parameter families, each one
number for the uniform case with a per-entity array that wins when
allocated:

    coefficient   c_e   applied at the innermost step
    spacing       h_e   the length of an edge, in G
    measure       m_v   the size of a cell, in D
    boundary      b_e   the value standing in for the missing end
                        of an edge with no head

## the dictionary

old assembler read            tower seat
--------------------------    ------------------------------------------
face_areas    (farea)         folded into the per-edge coefficient
face_deltas   (fdelta)        spacings h_e
cell_volumes                  measures m_v
keff = n^T K n                per-edge diffusion coefficient, with area:
                              c_e = keff_e * farea_e
adv weights   (wp, wn)        the one-sided average step: upwind is the
                              sign of the coefficient picking the end
cell/face centers, lvec,      construction-time inputs to the skewness
tangents t1 t2                correction (below); never stored on an
                              operator
vertex/face cell weights      interpolation coefficients, order 0

## the interior stencil

the old two-point orthogonal jacobian rows

    diag  = -farea*( keff/fdelta + wp)
    neigh = -farea*(-keff/fdelta + wn)

decompose as two balance terms:

    diffusion   the order-2 operator with c_e = keff_e*farea_e,
                h_e = fdelta_e, m_v = volume_v
    advection   the order-0 average step with signed coefficient
                c_e = farea_e*vn_e, one-sided by upwind — wp and wn
                are the split of vn the sign already encodes

the balance sums both onto the cells through incidence, once per face.

## boundary conditions

every condition is stored as robin, a*phi + b*dphi/dn = c, per tag.
the face value is eliminated one-sidedly:

    phi_b = (c + (b/delta)*phi_p) / (a + b/delta)

which yields, per boundary face (verified against
class_boundary_condition.f90):

    diffusive diagonal   -kappa*area*a / (delta*(a + b/delta))
    diffusive rhs        -kappa*area*c / (delta*(a + b/delta))
    advective            area*vn*phi_b, split the same way into a
                         phi_p coefficient and a constant

in the tower: one operator instance per tag, aimed at
`tagged_edges(tag)`, its diagonal part riding as edge coefficients on
the edges without heads, its constant part a source term in the
balance. the boundary seat b_e carries the eliminated face value where
a chain needs one.

## skewness

the old assembler adds a deferred correction on skewed faces: a
minimum-correction vector plus a least-squares tangential gradient,
diffusion only (get_skew_source). in the tower this is one more edge
term in the balance, its coefficients computed at construction from
centers and tangents, re-evaluated per outer iteration the same way
the old source was. it does not change the stencil; it adds to the
right hand side.

## what rung 2 must build

the mesh, extending stored_graph, answering the typed questions the
dictionary needs:

    cell_volume()    real field on cells
    cell_center()    real field on cells, 3 wide
    face_area()      real field on faces
    face_delta()     real field on faces
    face_normal()    real field on faces, 3 wide
    face_center()    real field on faces, 3 wide
    face_weights()   real field on faces (interpolation)

loaders fill these at construction. every operator above is then two
lines: fetch the fields, hand their arrays to a constructor.

## the acceptance checks

- the order-2 operator with the dictionary's coefficients reproduces
  the old two-point row on a hand mesh, entry for entry
- the robin coefficients above, per tag, match the old
  boundary_condition functions to machine precision
- the balance conserves: interior terms sum to zero over the mesh,
  boundary terms sum to the boundary flux
