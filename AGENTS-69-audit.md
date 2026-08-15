# §69 audit of `src/`

Applied 2026-08-14 to production sources only, branch
`tower-graph-as-sets-relations`.  Towers under `test/` were read for evidence
but are not audited here.

Verdicts use the vocabulary of `AGENTS-69-candidate.md`.  One finding is
checked by a runnable probe; the rest are checked by reading, and each names
the file and line that decides it.

```text
V1   69.6/69.3   a query mints an identity                     CHECKED   severe
V2   69.4        domain() carries two meanings                 reading   severe
V3   69.3        three operations take their domain from context reading major
V4   69.3/69.1   resolved authority re-encoded as a flag       reading   major
V5   69.5        fabricated host sets                   reading   major
V6   69.1        cells/faces are aliases, not domains          reading   minor
V7   69.1        three lexical conventions for domain labels   reading   minor
V8   69.1        §0.1's canonical example names nothing real   reading   minor
V9   69.2        graph arguments declare no role, 143 sites    reading   major
V10  69.7        a completeness claim outlived its subject     reading   major
```

---

## V1 — a query mints an identity.  69.6, severe, CHECKED

`class_graph_reduction.f90:556`

```fortran
subroutine reduction_domain(this, input_graph, domain)
    ...
    allocate(domain, source=index_set('answer', 1))
```

`index_set(...)` is `create_index`, which calls `declare`, which mints a
fresh token (`graph_set.f90:299-307`).  Declaration is an act of
ownership.  This is a query, and it performs one on every call.

The consequence is not stylistic.  A pure query returns a different
mathematical object each time it is asked:

```text
$ ./probe
 set: two identical declarations equals     = F
 reduction: domain() equals itself across calls = F
```

```fortran
! probe.f90 - links against lib/libufvm.a
type(reduction)   :: r
type(stored_graph):: g
class(set), allocatable :: d1, d2

g = stored_graph(2, tails=[1], heads=[2])
call r % domain(g, d1)
call r % domain(g, d2)
print *, d1 % equals(d2)          ! F
```

Downstream, every consumer that routes by identity is defeated.  The
partitioner classifies a field by embedding
(`class_graph_partitioner.f90:481`), so a field on a reduction's domain
reaches

```text
error stop 'partition: this field does not live on this graph''s domains'
```

and two reduction answers can never be compared for domain equality — nor can
one reduction's answer with itself.

**The correct pattern is already in production, two files away.**
`class_graph_functional.f90:157` signs once, at construction, and stores it:

```fortran
this % home = index_set('functional', 1)     ! constructor
allocate(domain, source=this % home)           ! query, line 201
```

Repair: give `reduction` a `home` component signed in its constructor.

Law added by this finding: **a query may not sign.**  Section 69.6 states that
exactly one object may construct a structure; it did not say that construction
must happen at construction time.  It should.

---

## V2 — `domain()` carries two meanings.  69.4, severe

The two-implementation test over `operation_domain_interface`
(`graph_grammar.f90:533`):

```text
walk         % domain : all vertices of the graph handed in       where I READ
balance      % domain : all vertices of the graph handed in       where I READ
differential % domain : all edges or all vertices, by landing     where I LAND
reduction    % domain : a fresh one-entry set                 where I LAND
functional   % domain : this % home                               where I LAND
step         % domain : this % action % domain(...)               delegated
```

The sentences differ in far more than an axis word.  Two questions are being
answered by one verb: the domain an operation *reads over* and the domain its
answer *lands on*.  For `walk` and `balance` these coincide, which is why the
conflation survived.

Verdict `two meanings`.  Either split the binding — the honest reading, since
an operation genuinely has both a domain and a codomain — or name the one it
owes and make every concretion answer that one.  Note that `apply` already
returns `output`, so the codomain reading is the one with a consumer.

---

## V3 — three operations take their domain from the execution context.  69.3, major

```text
class_graph_walk.f90:143                 call input_graph % all_vertices(domain)
class_graph_balance.f90:132              call input_graph % all_vertices(domain)
class_graph_differential_operator.f90:432  input_graph % all_edges / all_vertices
```

A structural fact is read from the graph the action executes over.  This is
the pattern `graph_minimization` was repaired to refuse — "the unknown domain
U: where the answer lives, explicit at attach, identity preserved — never
inferred from the host" — and `class_graph_step.f90:199` names it outright:
"that domain, which asking the graph never could."

The defect is not that the answer is wrong today.  It is that context is
serving as authority, so the answer changes when the host changes, silently.

Verdict `inferred`.  Admissible repairs per 69.3: a declared seat with no
fallback, or a loud refusal.  If an operation's domain genuinely is the whole
host, the caller states that equality at its own call site.

---

## V4 — a resolved authority is re-encoded as a flag.  69.3 and 69.1, major

`class_graph_partitioner.f90:481` and `class_graph_assembler.f90:219` resolve
the question correctly, with the mechanism the architecture prescribes:

```fortran
if (dom % is_subobject_of(global_graph % vertex_set())) then
   call carry_field(global_data, dom, global_graph % vertex_set(), &
        &           part_graph, .true., part_data)
else if (dom % is_subobject_of(global_graph % edge_set())) then
   call carry_field(..., .false., part_data)
```

and then discard the answer, re-encoding it as an unlabelled boolean literal
that travels through four routines — `carry_field`, `gather_field`,
`global_of`, `owner_of` — each of which re-derives the set from the flag:

```fortran
logical, intent(in) :: on_vertices
if (on_vertices) then
   part_set = part_graph % vertex_set()
else
   part_set = part_graph % edge_set()
```

`dom` is passed to every one of those routines.  The declared supplier is in
the room and a boolean is consulted instead.

Two of the document's own laws name this:

- §6 — the embedding query, "never a side flag, never concrete-type
  inspection — is how routing decisions are made";
- §39 — "the new design should go beyond a two-valued `side`. ...  Do not
  replace `vertex/edge` with merely `A/B` and declare the problem solved."

`on_vertices` is `side` spelled as a logical.  A call site reading
`carry_field(global_data, dom, set, part_graph, .true., part_data)` cannot
be understood without opening the callee.

Verdict `inferred`, plus `alias` under 69.1 for the positional word.

---

## V5 — fabricated host sets.  69.5, major

```fortran
class_polynomial_form.f90:41
  this % subset = subset('polynomial-basis', &
       &                         index_set('polynomial-table', 4), [(m, m=1,4)])

class_harmonic_form.f90:50
  this % subset = subset('harmonic-basis',   &
       &                         index_set('harmonic-table', 3), [(m, m=1,3)])
```

Each host is constructed inline, sized to exactly equal the subset that
embeds in it, never stored, never named elsewhere, never asked a question.
It exists because `subset` is the only constructor admitting a chosen
roster, and it requires an ambient to embed in.

The fabrication test: the object exists so another type will compile.  It
fails all five clauses of §50.

The root cause is a missing concretion, and `graph_set.f90:145` already
predicted it: *"A domain that must list its members arrives as a second
concretion the day something needs it."*  That day was `polynomial-table`.

Verdict `fabricated` → `taxonomy review`, per 0.3.  The repair is an
architectural decision, not a patch: admit a set that holds an explicit
roster, or give these hosts a mathematical reason to exist.

---

## V6 — `cells` and `faces` are aliases.  69.1, minor

`class_graph_mesh.f90:131`

```fortran
cells = this % vertex_set()
faces = this % edge_set()
```

Two words, one token.  §44 reserves `cells` and `faces` for `X₃` and `X₂` of a
mesh's own sets; here they are local names for a graph's vertices and
edges.  Fields are then declared on them (`cell_volume`, `face_area`), so the
physics vocabulary reaches users while the identity is the graph's.

Verdict `alias`.  Harmless today, and it is the first thing that breaks when
a mesh grows a third set.

---

## V7 — three lexical conventions for domain labels.  69.1, minor

```text
snake    interior_vertices  owned_edges  overlap_vertices
kebab    polynomial-basis   harmonic-table
bare     answer  numbers  functional  basis  sources  sinks
```

§0.3: "Pluralization, prefixes, and suffixes must not create a second
ontology."  Three spellings for one class of name is a second ontology in the
weakest place — the one Section 0 does not currently reach, because it governs
identifiers and these are strings.

Separately, `numbers` (`graph_minimization.f90:198`) and `answer`
(`class_graph_reduction.f90:556`) are vague words where role words exist:
these are the unknown coordinates and a functional's value.

Verdict `rename`, reported and not patched.

---

## V8 — §0.1's canonical example names identifiers that do not exist.  69.1, minor

Section 0.1 presents, as *the* canonical statement for "the operation boundary
discussed by this architecture":

```text
gstructure : graph       = structural object over a set domain
ginput     : graph_field = value object over a set domain
goutput    : graph_field = value object over a set domain
```

`gstructure`, `ginput`, `goutput` occur zero times in `src/`.  The actual
boundary is `graph_grammar.f90:540`:

```fortran
subroutine operation_apply_interface(this, input_graph, input_data, output)
```

§0.6 clause 9 requires architectural documentation to use the canonical names.
Either §0.1 prescribes a rename nobody performed, or its example is fictional.
Section 0 must not be the first violation of Section 0.

Verdict `rename` — of the document or of the code, architect's call.

---

## V9 — graph arguments declare no role.  69.2, major

143 `class(graph), intent(...)` dummy arguments in `src/`.  Under the 69.2
normal form `<name> : graph = <role> for <what>`:

```text
coupling      authority for the dependent axis            declared
on            execution context for the action            declared in comment only
input_graph   ??? for ???                                 position, not role
g             ??? for ???                                 8 sites, a bare letter
```

`input` is a canonical role word for a *value* crossing an interface; for a
structural argument it states only that it arrived.  Whether a given
`input_graph` is operand, context, or authority is exactly what the reader
must recover, and it currently varies by callee — `walk` reads structure from
it (operand), `step` passes it through untouched (context).

`g` (`class_graph_differential_operator.f90`, 8 sites) states nothing at all;
§0.2 requires a variable name to carry the role of the object in the present
expression.

Verdict: `declared` at 1 site, unlabelled at the rest.  This is the largest
finding by count and the cheapest to fix, since the roles are already known
where the code is correct.

---

## V10 — a completeness claim outlived its subject.  69.7, major

AGENTS.md §37 states, of the support migration:

> ALL STAGES COMPLETE as of 2026-08-12 — `class_graph_support` is deleted,
> field domains are set, **the six side() routing sites are gone**

The `side()` *procedure* is gone.  The routing is not:

```text
graph_calculus.f90:99                       GRAPH_SIDE_VERTEX = 1
graph_calculus.f90:100                      GRAPH_SIDE_EDGE   = 2
class_graph_differential_operator.f90:221   integer :: landing = GRAPH_SIDE_VERTEX
class_graph_differential_operator.f90:432   if (this % landing == GRAPH_SIDE_EDGE)
class_graph_differential_operator.f90:813   if (this % landing == GRAPH_SIDE_EDGE)
class_graph.f90:57, class_graph_reduction.f90:86, graph_forms.f90:44   importers
```

The *code* is candid about this — `graph_calculus.f90:93` says plainly that
these are "the legacy differential operator's output-LANDING choices ... not
field-domain identity."  The architecture document is the party making the
unqualified claim.

Under 69.7 this is `asserted`: no artifact fails if a side constant spreads
to a seventh site.  Together with V4 — where the two-valued routing survived
the migration as a logical — the claim is checked by proxy at best.

Repair is a checker, not a rewrite: a gate that fixes the permitted importers
of `GRAPH_SIDE_*` and fails on a new one, selftested against a planted import.

---

## What passed

An audit that lists only failures is not an audit.

```text
69.1  transport re-declaration      restated-by-law   partitioner:505, tied to round trip
69.2  minimizer % coupling          declared          the model citizen; no fallback
69.3  minimizer % unknown_domain    declared          explicit at attach
69.4  step vs stencil dependencies  one meaning       differ only in the axis word
69.6  functional % home             single owner      signed at construction, stored
69.7  the numberless law            checked           reads sources, selftested
```

Five of the six passes are recent — the laws are being obeyed where they were
most recently applied, and violated where the code predates them.  That is
the ordinary shape of an architecture in migration, and it says the direction
is right.

## Ledger

```text
Names       : cells/faces alias · three label conventions · gstructure absent
Roles       : declared at 1 of 144 graph arguments
Authority   : inferred at walk, balance, differential; flag-encoded in transport
Axes        : domain() carries two meanings
Abstraction : polynomial-table, harmonic-table fabricated
Ownership   : reduction % domain signs on every call  (CHECKED)
Verification: §37's side() claim asserted, not checked
```
