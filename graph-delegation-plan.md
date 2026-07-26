# graph-delegation-plan: move the repeated edge-building loop into the graph interface

status: implemented. all test suites pass (graph 90 checks, orbit,
amg, mesh, assembly, solve, krylov, advection, coupled_diffusion,
regression, unsteady, integrator, adjoint, spatial-discretization,
and the distributed parallel suite). actual size below.

## the problem

five places in the code build an edge list the same way: loop over
vertices, use some rule to find which other vertices each one connects
to, skip duplicates with a marker array, count the edges on a first
pass, allocate, then fill on a second pass. the rule is different in
each place, but the surrounding loop is copied five times. two of the
copies are in interface_graph itself.

## where the copies are

- `interface_graph.f90:1025` quotient_edges - for each partition part,
  loop over its owned vertices and their neighbours; record an edge
  between two parts when a neighbour lies in a higher-numbered part
- `interface_graph.f90:2160` condensation_edges - same loop over
  strongly-connected components and their out-neighbours; record an
  edge when a neighbour lies in a different component
- `class_mesh.f90:1780` node_graph - two cells get an edge when they
  share a mesh point
- `class_csr.f90:602` simple_graph - two rows get an edge when either
  matrix entry (i,j) or (j,i) exists, recorded once per pair
- `class_algebraic_multigrid.f90:68` coarsen - two rows get an edge
  when the matrix entry between them is strong

each site is about 20-30 lines of identical loop machinery wrapped
around a 3-6 line rule.

## the fix

add one subroutine to the graph interface that owns the loop, and let
each caller pass in its rule as a procedure argument. the interface
already does this in two places: orbit takes a successor function, and
accumulate_adjoint takes an edge_apply subroutine. same pattern here:

    subroutine harvest_edges(this, rule, tails, heads, group, directed)
      ! loop over the vertices; rule(v) returns the candidate
      ! neighbours of v; keep each new edge once, using a marker
      ! array to skip repeats. if `group` is given, edges are
      ! recorded between groups instead of vertices (this is what
      ! quotient_edges and condensation_edges need). `directed`
      ! controls whether an edge is kept once per pair or once per
      ! direction.

each caller then keeps only its rule, written as a small contained
procedure (a contained procedure can see the caller's variables, so
the mesh rule can read the mesh and the csr rule can read the matrix):

- quotient_edges: rule = neighbours, group = the partition
- condensation_edges: rule = out_neighbours, group = the components
- node_graph: rule = cells sharing a mesh point
- simple_graph: rule = out-edges filtered by the symmetry check
- amg coarsen: rule = out-edges filtered by the strength test

no new interface file is needed; this is one new method on the
existing graph type.

## a second, smaller cleanup

three places contain the same five lines for "given a face and one of
its two cells, find the cell on the other side":

- `class_fvm_field.f90:111`
- `class_assembler.f90:498`
- `class_assembler.f90:614`

replace all three with one pure mesh function:

    ncell = grid % across(icell, gface)

## checked and found already clean

these looked like candidates during the scan but need no change:

- the stationary solvers (jacobi, gauss-seidel, sor) - all matrix
  access already goes through the assembler interface
- the partitioned assembler - already calls owned, ghosts, gather,
  scatter, and dofs_of on the graph instead of doing its own loops
- multigrid prolongation - tentative_prolongation is already shared
  by both multigrid classes
- bfs, transpose_adjacency, counting_sort, power_iteration - each has
  one implementation and all callers use it
- the two-pass loop in mesh agglomerated - it sizes polygons, not
  edges, so it does not fit this kernel
- the two-pass loop in ghost_copies - it collects (key, image, index)
  triples, not edges; forcing it into the kernel would complicate the
  kernel for one caller
- csr add and matmat - two similar count-then-fill row merges in one
  file; they could share a kernel, but merging them trades two
  readable functions for one harder one. noted, not proposed.

## actual size of the change

198 lines added, 160 removed: net +38, not the deletion the estimate
promised. where the estimate went wrong:

- the kernel carries a 14-line banner and the interface a new list
  comment; the two converted callers inside interface_graph also
  gained short delegation notes. comments are most of the growth.
- each contained rule costs its declaration, contains, and end lines
  (about 8 per caller). for quotient, condensation, and node_graph
  the removed walk was much bigger than that, so they shrink. for
  simple_graph and amg the rule keeps nearly the whole old filter
  loop, so those two files grow slightly.

what the change buys despite the count: the count-allocate-fill walk
and the mark-stamp dedup now exist in exactly one place. a bug in
that machinery is now one fix; a new graph builder is now just a
rule. if the line count matters more than that, reverting only the
csr and amg conversions recovers most of the growth while keeping
the three big wins.

## order of work

1. add harvest_edges and convert the two callers inside
   interface_graph (quotient_edges, condensation_edges); graph tests
   verify
2. convert node_graph, simple_graph, and the amg strength loop; mesh,
   orbit, and amg tests verify
3. add across to the mesh and convert the three call sites; assembly,
   regression, and parallel tests verify

each step compiles and passes on its own, so work can stop after any
of them.
