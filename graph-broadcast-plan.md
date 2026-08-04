# Plan: graph_broadcast

Written 2026-08-03. Implementation follows approval.

## The contract addition

`graph_broadcast` is the transpose of `graph_reduction`, and the
pair is visible in the signatures: `reduce` takes a field and
yields a functional; `broadcast` takes a functional and fills a
field. One abstract type, one deferred procedure, three arguments -
the operation itself and the two data of the map. The contract
becomes 92 deferred procedures, 23 types.

    type, abstract :: graph_broadcast
     contains
       procedure(broadcast_interface), deferred :: broadcast
    end type graph_broadcast

    pure subroutine broadcast_interface(this, functional, field)
      import :: graph_broadcast, graph_functional, graph_field
      class(graph_broadcast) , intent(in)    :: this
      class(graph_functional), intent(in)    :: functional
      class(graph_field)     , intent(inout) :: field
    end subroutine broadcast_interface

No graph argument and no support argument: a field is constructed
on its support and knows its own entry count and side, so the
caller's field brings everything. No measure argument: the
transpose of the weighted sum is the order-0 operator applied to
the copy, a composition of two existing pieces.

## The concrete

It lives in `class_graph_reduction.f90`, beside its pair. One type,
two rules, following the reduction's one-type-with-a-rule pattern.

    integer, parameter :: BROADCAST_COPY  = 1   ! transpose of sum
    integer, parameter :: BROADCAST_SHARE = 2   ! transpose of average

    type, extends(graph_broadcast) :: broadcast
       integer :: rule = BROADCAST_COPY
     contains
       procedure :: broadcast => broadcast_functional
    end type broadcast

The type name and the verb coincide, so the specific procedure
carries the descriptive name `broadcast_functional`; the binding
keeps the plain verb. The rule component is public, so the
intrinsic structure constructor `broadcast(BROADCAST_SHARE)`
suffices and no constructor interface is written.

The body, stated as its cases:

    n = field % num_entries()
    J = the functional's value, by its live kind
    COPY  : every entry receives J
    SHARE : every entry receives J / n
    a real J fills a real field; a complex J fills a complex field,
    which carries a complex-step seed; any other kind fills zeros -
    the value-kind rule the fields already state

## The reduce correction

The same audit applied to the existing pair finds two defects in

    reduce(this, input_graph, field, support, functional, measure)

`input_graph` is never consulted. `support` is accepted but not
honored: accumulate loops over every entry of the field, so a
subset support is silently ignored - an unkept promise. Subsets are
served by constructing the field on the subset support, which is
what the partitioned tests already do. `measure` is load-bearing -
it is the second field of the inner product - and stays.

The corrected interfaces, amended in the same opening of the
contract:

    reduce(this, field, functional, measure)
    accumulate(this, field, state, measure)

`initialize`, `combine`, and `finalize` are already clean. The
concrete and every call site in the suite change accordingly. No
procedure is added or removed by this correction, so the counts
remain 92 and 23 with the broadcast included.

## The checks

Four, in a subroutine `check_broadcast` after
`check_inner_products`, on a four-entry support with J = 6:

    round trip, copy    reduce_average(copy of 6)   = 6
    round trip, share   reduce_sum(shares of 6)     = 6
    pairing             with u = [1,2,3,4]: the inner product of
                        the copy of 6 with u is 60 = J * sum(u),
                        the defining property of a transpose
    complex seed        the copy of 6 + 0.001i, summed, returns
                        24 + 0.004i, both parts intact

## Verification

A standalone stub binds all 92. Rebuild, install, rebuild the
contract suite, run: 169 checks pass (165 committed plus four).
Nothing is committed or pushed without explicit instruction.

## Alternative recorded and set aside

The transposes are expressible without a type, as constant fields
built by callers. Set aside because it leaves the reduction without
its pair in the taxonomy and repeats the seed construction at every
adjoint call site.
