# Independent graph execution

The driver owns its position in the dependency order. `advance` executes one
scheduled vertex and releases paired data after that position's final reads.
`next_rule` returns the next vertex identity; `num_completed` returns the number
of completed positions. These differ when vertex identities are permuted.

`set_rule` and `clear_rule` modify one stored operation without copying the
pairing or changing its incidence, order or release intervals. A caller can
therefore bind transient execution references immediately before an advance
and remove those references immediately afterward.

`advance_with` executes the next scheduled vertex with a rule the caller
passes instead of one stored in the pairing: the driver states the rule's
vertex and the vertices its arguments read (`vertex_rule`), declares one
argument per read when the rule declares another count, and performs the
same writes and final reads as `advance`. `expired_at` returns `released_at`
followed by the data the step's vertex writes that no rule reads, so a
caller retaining outputs reads one lifetime for both.

`evaluate` retains full traversal semantics: it resets the position and uses
the same `advance` implementation until complete. It does not reconstruct data
released during an earlier traversal. `pair_with` supplies a new data branch
and resets the position for a new execution.

The suite verifies independent interleaved and copied executions, nontrivial
vertex order, incremental rule replacement, complete replay, unbound argument
identities, absent-rule positions, fresh output per operation, and exact final
read intervals for a combined forward/reverse graph. The initial source must
remain available until the final reverse rule consumes it. Temporal minimizers
delegate these operations to their own drivers and report completion only once
the entire order has been evaluated.

Run `./run.sh`, or use the root `verify.sh` after building the library.
