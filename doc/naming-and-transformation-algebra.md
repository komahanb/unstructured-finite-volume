# Naming algebra, and the transformations classified by law
Analysis only. Nothing is renamed and nothing is migrated. Read after the member-set and
relation maps.

## 1. The naming algebra
```
graph_<mathematical noun>_<role>
```
Roles are closed: **view**, **map**, **representation**, **algorithm**, **algebra**, **calculus**.
A name is admitted only when both halves are stated and neither needs a paragraph.

| role | answers | keyed on | may be large |
|---|---|---|---|
| *(graph)* | which object? | — | no |
| **view** | what do its parts mean here? | — | no |
| **map** | what is associated with this identity? | identity | yes |
| **representation** | how is this stored and enumerated here? | — | **yes** |
| **algorithm** | what is computed from it? | — | — |
| **algebra / calculus** | what laws hold, and what operations respect them? | — | — |

The one architectural rule that follows: **only maps and representations may be O(n).** Views
and graphs are the semantic control plane and stay small. This is what makes 10⁹ cells possible
without 10⁹ graph objects, and it is a law, not a preference.

Rejected nouns, and why each is not a mathematical object: `profile` (an interpretation — say
which), `ordinary` (an adjective), `carrier` (a metaphor for *set*), `support` (retired already),
`kernel` (overloaded; `fractal_graph` stands), `helper`, `manager`, `wrapper`. `fractal` is
permitted: it denotes self-similarity, exactly.

## 2. The naming table
Confidence is evidence-backed: **high** = demonstrated by a prototype or by production code read
in this pass; **medium** = argued, not yet built; **low** = the shape is clear, the name is not.

| current | mathematical object | role | proposed | conf | action |
|---|---|---|---|---|---|
| `member_set` | finite set | representation contract | `graph_set_representation` | high | migrate |
| — (its `token`) | set identity | graph | `type(graph)` | high | migrate |
| `counted_set` | set with members 1..n | representation | `graph_set_counted_representation` | high | migrate |
| `subset_set` | a set, and S ⊆ A | — | **delete**: set graph + subset relation | high | delete |
| `subset_set % roll` | explicit member list | representation | `graph_set_listed_representation` | medium | migrate |
| `is_subobject_of` | transitive inclusion | algorithm | `subset_closure` (name low) | medium | migrate |
| `relation` | R ⊆ A₁×…×A_k | representation contract | `graph_relation_representation` | high | migrate |
| — (its `token`) | relation identity | graph | `type(graph)` | high | migrate |
| `slot` | one signature position | — | **delete**: a sequence of set graphs is homogeneous | high | delete |
| `arity`, `domain(k)` | signature length, k-th domain | view | `graph_relation_view` | high | migrate |
| `stored_relation` | tuple table | representation | `graph_relation_table_representation` | high | migrate |
| `binary_relation` | — | — | **delete**: arity is a signature fact | high | delete |
| `source`, `target` | domain 1, domain 2 | view | `graph_binary_relation_view` | high | migrate |
| `image`, `preimage` | fibres | algorithm | in the binary view, over a representation | high | migrate |
| `csr_relation` | compressed sparse rows | representation | `graph_relation_csr_representation` | high | migrate |
| `transposed_view` | slot permutation | view | `graph_relation_transpose_view` | high | migrate |
| `materialized()` | — | — | **delete**: a role marker | high | delete |
| `inclusion_of` | I_S ⊆ S × A | algorithm → representation | derivable; build only on demand | medium | defer |
| `relational_binding` | element ↦ legacy object | map | `graph_object_map` (name low) | medium | defer |
| `graph_profile` | the T/H schema reading | view | `graph_incidence_view` | medium | defer |
| `directed_incidence_view` | V/E graph over T,H | view | `graph_ordinary_view` | low | defer |
| `directed_adjacency_view` | A ⊆ V×V reading | view | `graph_adjacency_view` | medium | defer |
| `partitioner` | constructs ρ + parts | algorithm | `graph_partition_algorithm` | medium | defer |
| `assembler` | traverses ρ backwards | algorithm | `graph_assembly_algorithm` | medium | defer |
| `vglobal`/`eglobal` | local ↦ global | map | `graph_locality_map` | medium | defer |
| `vowner`/`eowner` | local ↦ part | map | `graph_ownership_map` | medium | defer |

`relational_binding` is deliberately **defer**: it is scaffolding for objects that are not yet
graphs, and it should be renamed when it is clear whether it survives the set and relation
migrations at all.

## 3. Transformations, classified by law
Only the classes the task defines, and only where evidence exists. **"proven" means a test in
this repository asserts it.**

| transformation | law | class | proven? |
|---|---|---|---|
| relation transpose, extension | T(T(R)) has R's tuples | **A involution** | **yes** — `graph-binary`, and `test/fractal-map/relation.f90` §4 |
| relation transpose, identity | T(T(R)) *is* R | A involution | **no today** — `transpose_of` mints a fresh token; the role-permutation prototype reaches it |
| dual | dual(dual(G)) = G | **A involution** | **not implemented at all** — AGENTS.md prose only |
| partition / assembly | `assemble ∘ partition = I` | **C section/retraction** | **only at nparts = 1** |
| refine / coarsen | `coarsen ∘ refine = I` | **D one-sided** | not proven either way |
| primal / adjoint | ⟨Ax,y⟩ = ⟨x,A*y⟩ | **E adjoint** | **yes** — derivative-action tower Gate B, duality 18 = 18 |
| inclusion I_S | S ↪ A, total functional injective | not a pair | by construction |

Four things this table refuses to say:

- **Transpose is not adjointness.** Transpose permutes slots of a signature; adjointness is an
  identity between inner products. The derivative tower proves the second; `graph-binary` proves
  the first. They are different laws that happen to reverse arrows.
- **Dual is not transpose.** `dual(dual(G)) = G` is an involution — the same *law class* as
  transpose — and that is all they share. Nothing implements dual; it stays a stated law for a
  future view, over a schema that defines which incidence is exchanged.
- **Assembly is not the inverse of partition.** See §4.
- **Refine/coarsen is unclassified.** Class D is the *hypothesis*; the test that would settle it
  does not exist.

## 4. Partition and assembly: what is actually proven
`test/graph-contract/test.f90` proves `assemble(partition(G)) == G` and
`assemble(partition(G,D)) == (G,D)` at **nparts = 1**, and the test itself calls that "the
weakest case". For nparts > 1 the partitioner hands back **one** part, so a single-part round
trip cannot be the identity — what is proven there instead is *ownership*: every whole-graph
cell is owned exactly once across the parts, under three rules and four part counts, and each
part borrows the cells across its cut.

So by measured evidence the pair is **class C, a section/retraction**: `partition` is a section
of `assemble` on a single part, and the general statement — assembly from *all* parts
reconstitutes the whole — **has no test**. Promoting it to class B needs exactly that test:
partition into k parts, assemble from all k, compare against G. That test should exist before
any of this migrates.

## 5. What ρ is, and where it currently lives
The whole↔part correspondence exists today as four arrays and two scalars **inside**
`stored_graph`:

```fortran
integer, allocatable :: vglobal(:), eglobal(:)   ! local -> global
integer, allocatable :: vowner(:),  eowner(:)    ! local -> owning part
integer :: nparts ;  logical :: cut
```

That is ρ living inside the object it describes. Decomposed by role, it is **four things**:

| thing | role | today |
|---|---|---|
| local ↦ global correspondence | **map** (functional, injective) | `vglobal`, `eglobal` |
| ownership: local ↦ part | **map** | `vowner`, `eowner` |
| owned / borrowed membership | **view** over the ownership map | `owned_vertices`, `borrowed_vertices` |
| the cut itself: how many parts | **view** state | `nparts`, `cut` |

`owned_vertices` and `borrowed_vertices` return `subset_set`, computed on demand from `vowner` —
an attribute map presented as a subset. Under the set map that becomes a view over a map, and no
subset object is built at all.

And the partition/assembly *algorithms* are a fifth and sixth thing, distinct from ρ and from
each other. §20's demand is met by this table: **partition relation, partition algorithm,
assembly algorithm and communication representation are four objects, not one type.** The
communication representation — who exchanges what with whom — is not present in this repository
at all today, and must not be conflated with the ownership map when it arrives.

## 6. Order of the remaining migration
Chosen by evidence available and blast radius, not by tower level.

1. **`subset_set` → set graph + subset relation.** One construction gate, five overrides, and
   every consumer already asks `is_subobject_of` rather than reading `host`. Blocked on nothing
   except the relation map's signature decision.
2. **`slot` → sequence of set graphs.** Mechanical once sets are graphs; deletes a wrapper type
   whose whole reason was Fortran's one-dynamic-type-per-array rule.
3. **`materialized()` → deleted**, when representation and view become separate modules.
4. **`binary_relation` → a view.** Larger, because `image_view`/`preimage_view` are the hot path
   and the deferred-primitive split must survive intact.
5. **Partition/assembly.** Last, and only after the missing assemble-from-all-parts test exists.

`interface_graph.f90`'s abstract `graph` is not in this list: it is an interface to the legacy
stack, and it migrates when that stack does.
