# the next strokes

what a six-lens review of src/ and test/ found, and the order i would
mend it in. sixteen findings survived adversarial cross-examination;
four were refuted, and two of those refutations change the rulings we
had already made, so they are written up here beside the defects.

nothing below is built. this is the map.

## the three defects that are not style

    class_diffusion_statement.f90:69 ── a wall's a and b never leave
                                        the condition. the statement
                                        reads faces() and
                                        boundary_values() and nothing
                                        else, so
                                          neumann('w', q)  builds the
                                          operator of dirichlet('w', 0)
                                        a zero-FLUX wall is compiled as
                                        a zero-VALUE wall, and
                                          robin('w', 2, 3, 5)
                                        is byte-identical to
                                          dirichlet('w', 2.5)
                                        the suite never catches it
                                        because no neumann has ever
                                        reached a statement: constitution
                                        and minimization both test two
                                        dirichlets and stop.

    graph_minimization.f90:154, :318 ── level 2 holds a level-1 field
                                        that knows its own width and
                                        reads it without asking
                                        num_components(). the word
                                        appears nowhere in the whole
                                        family. so an implicit march of
                                        a multi-component state is
                                        broken TODAY - the marcher hands
                                        its step operator to
                                        inner % attach, and everything
                                        past that point counts one value
                                        per cell. no test catches it: the
                                        multi-component fixtures all
                                        march explicitly.

    class_graph_marcher.f90:174 ─────── march_adjoint never reads
                                        this % rule. the forward walk
                                        solves an implicit step per edge;
                                        the reverse walk always does
                                        explicit euler. the pairing the
                                        banner promises across every step
                                        holds for one rule of three.
                                        latent: both call sites march
                                        forward. do NOT simply delete the
                                        verb - two suite checks stand on
                                        it, and they are the ones that
                                        retired the complex step.

and one contract that lies:

    defined_on_data ─── the data gate on every transform. four
                        implementations, ZERO callers. two of the four
                        bodies are tautologies - the dummy is already
                        class(graph_field) and the guard asks
                        class is (graph_field), so the class default arm
                        is unreachable and num_entries() >= 0 is
                        vacuous. meanwhile the roads it is supposed to
                        guard (assemble_data, partition_data) narrow to
                        the one concrete with no fallback, so any other
                        field returns an UNALLOCATED intent(out) to the
                        caller. the gate says yes to everything and the
                        road serves one thing.

## the law the grammar states and the code breaks

graph_grammar.f90:131 asks "CAN A GRAPH CHANGE?" and answers "No.
Everything a graph holds - structure, tags, its relation to the whole
it came from - goes in at construction". two citizens break it, and
both were FORCED to, because no door exists:

    class_graph_partitioner.f90:262 ── writes seven components of a
                                       stored_graph after building it:
                                       cut, nparts, me, vglobal, vowner,
                                       eglobal, eowner. the constructor
                                       takes (nv, tails, heads, vtags,
                                       etags, number) and nothing else.
                                       the same object answers
                                       has_part_relation() differently
                                       before and after line 264 - the
                                       exact repeatability the grammar
                                       gives as its reason.
                                       (me is not missing: it duplicates
                                       number, and can die.)

    class_form_pruner.f90:78 ───────── writes a form's inherited support
                                       parent by name. form's contract
                                       has size_of, values, slopes and
                                       nothing that admits a membership
                                       change, so the pruner reached for
                                       the component. a fit holding that
                                       form is silently re-typed by a
                                       governor it never spoke to.

## the cull, counted

roughly 460 lines in src, all consolidation, none of it design:

    -111  defined_on_data: the deferred symbol, its interface block,
          four bindings and four bodies
    -106  class_graph_support: twelve byte-identical degenerate bodies
          that three procedures serve (six neighbourhood answers, three
          edge-set answers, three part-set answers). the review found
          two further pairs; call it sixteen across five groups
    -113  class_robin_condition: six published coefficient arrays that
          all generate from two numbers, w = a/denom and v = c/denom.
          five of the six have no consumer in src at all
     -56  ten hand-rolled identity-index fill loops where the house
          one-liner support(SIDE, [(v, v = 1, nv)]) already exists
     -52  four of functional's scalar getters, called from nowhere
     -24  conduction and advection carry byte-identical
          edge_coefficients bodies

and the censuses in both banners are wrong - they will be wrong again
after the cull, so they get recounted once at the end rather than
patched twice.

## what the refutations changed

TWO RULINGS NEED AMENDING. i had them wrong; the review caught it.

1. SUBSTITUTION CANNOT INHERIT WHAT WE SAID IT WOULD - as stated.
every inherited word of minimizer is sized to the attached graph's
vertex set: attach does nv = on % num_vertices(), and matvec, norm and
inner_product all build their seats from that support. a marcher whose
x is one instant wide cannot express a trajectory, and a marcher that
overrides attach to dodge it inherits nothing - the reparenting buys a
name, not a contract.

    the amendment: substitution attaches the CHAIN, not the space graph.

        vertices = instants        one per time level
        ncomp    = the whole space state at that instant
        x        = the trajectory
        solve    = walk the chain forward, one governed block per edge
        adjoint  = the same walk backward, a direction component

    this is what block forward substitution actually is, it makes the
    inheritance real, and it dissolves the march_adjoint bug at the
    root: once the trajectory is in scope there is a state to linearize
    against. it also makes the ncomp fix a hard prerequisite rather
    than a nicety - the chain vertex IS a wide entry.

2. THE ONE-SYMBOL LEVEL-3 ROLE FAILS ARITHMETICALLY - as stated.
a role whose concretes feed each other breaks on units: robin
multiplies by area itself, conduction's edge_coefficients already
returns keff*area, so composing them gives area squared.

    the amendment: area leaves the laws entirely. a law answers a
    FIELD and the concrete chooses its width -

        conduction  ── ncomp 1   the bare normal projection n.K.n
        advection   ── ncomp 1   the bare normal speed v.n
        robin       ── ncomp 2   the eliminated face relation,
                                 phi_b = (1-w) phi_p + v

    the statement multiplies by area, because area is geometry and
    geometry is level 1. one role, no new symbol, no unit clash - and
    the SAME change fixes the neumann bug, because a two-number answer
    is exactly what the wall needed and could not send.

## the order

    stroke 1  the defects ········ the wall's a and b; ncomp through
              (+ their checks)     the minimizer family; march_adjoint
                                   honest about its rule; the data
                                   gates and their unallocated road.
                                   each lands with the check that would
                                   have caught it - a neumann statement,
                                   a two-component implicit march, an
                                   adjoint under bdf2.

    stroke 2  the doors ·········· stored_graph's constructor gains the
                                   frame (vglobal, vowner, eglobal,
                                   eowner, nparts; me dies as a
                                   duplicate of number). form gains
                                   restrict(). the two write-after-
                                   construction sites become
                                   constructor calls. the grammar's own
                                   law becomes true.

    stroke 3  the cull ··········· ~460 lines, no design change, one
                                   recount of both censuses at the end.

    stroke 4  level 2 ············ two minimizers, fit and form; the
                                   marcher becomes substitution on the
                                   chain; form_minimizer reparents to
                                   graph_transform and needs stroke 2's
                                   restrict; multigrid stops select-
                                   typing and asks dependencies(), the
                                   symbol that today has zero callers.

    stroke 5  level 3 ············ the role, as amended: laws answer
                                   fields, statements own the geometry.
                                   robin's six arrays collapse to two
                                   numbers. the two class_-prefixed
                                   function modules become honest.

    stroke 6  the funeral ········ the native reader; 6,680 legacy
                                   lines die; gfortran-16 comes out of
                                   quarantine, since the legacy mesh is
                                   the only code it miscompiles.

strokes 1 through 3 need no ruling - they are defects, doors and
deletions. strokes 4 and 5 carry the amendments above and want your
word before anyone writes code.
