# graph-benchmark is baseline-red, and its carrier migration has not compiled

**Status:** open. Not fixed by the carrier cutover, deliberately.

## The failure

    Fatal Error: Cannot open module file 'class_graph_support.mod'
                 for reading at (1): No such file or directory
    make: *** [Makefile:35: bench.o] Error 1

`bench.f90` imports `class_graph_support`, deleted in `a3817e3` ("the field
finds its true home, and the last edgeless ghost leaves the tower"). The suite
has failed on exactly this since that commit. Nothing in the carrier cutover
touched it, and nothing in the cutover fixed it.

## The part that matters, and is easy to miss

`bench.f90` **was** migrated off `graph_carrier` — in `133d29c`, and on purpose
*before* the module was deleted, so that its failure kept its shape instead of
becoming a fresh dangling `use`. Four lines changed:

    use fractal_graph          , only : set_graph => graph
    use graph_set_representation, only : counted_set_representation
    use graph_set_map          , only : set_map

    type(set_graph)                 :: vcarrier, ecarrier
    type(set_map)                   :: sets

    call vcarrier % declare()
    call ecarrier % declare()
    call sets % bind(vcarrier, counted_set_representation(nv))
    call sets % bind(ecarrier, counted_set_representation(ne))

    rel = csr_relation('tail', vcarrier, ecarrier, tab, sets)

**None of it has been through the compiler.** The `class_graph_support` import
sits at line 33, above every changed line, so compilation stops before reaching
any of them. The signatures were checked by inspection against

    create_csr(name, source, target, table, sets)
        type(set_graph), intent(in) :: source, target
        type(set_map)  , intent(in) :: sets

and they agree — but inspection is not compilation, and this file should be
treated as unverified until `class_graph_support` is resolved.

## What closing this requires

Decide what `support` becomes in the new ontology — a declared identity with a
listed representation, most likely, which is what `graph-characterization`'s
pinned checks already say the destination is. Then rebuild `bench.f90` and
find out whether the four lines above are right.

Until then this suite is red for a reason that predates the cutover, and the
32-suite sweep reports `31 / 32` with this as the one failure.
