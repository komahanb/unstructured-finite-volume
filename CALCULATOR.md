# AGENTS.md — Calculator Tower TDD Anchor

## Mission

Build the calculator-tower demonstration described in `test/calculator-tower/README.md`.

The entire demonstration uses one persistent mathematical statement:

\[
\boxed{(2+3)\times4=20}
\]

The purpose is larger than demonstrating arithmetic.

The calculator tower is the **vertical acceptance test for the architecture**:

```text
members
  ↓
relations
  ↓
relation algebra
  ↓
relational graph
  ↓
graph calculus
  ↓
field calculus
  ↓
discretization
  ↓
minimization
  ↓
constitution
  ↓
statement
  ↓
20
```

Every level must add exactly one new mathematical capability while preserving every truth established below it.

The calculator tower is therefore:

1. a tutorial;
2. a demonstration;
3. a regression suite;
4. a dependency-stratification test;
5. a test-driven-development anchor for the framework.


---

# Migration interlock — calculator TDD is mandatory during the refactor

`AGENTS.md` and this calculator specification are complementary authorities:

```text
AGENTS.md       architectural contract and migration direction
CALCULATOR.md   vertical mathematical acceptance oracle
```

Neither silently overrides the other.

If a calculator truth cannot be expressed using the architecture prescribed by
`AGENTS.md`, stop and report the conflict. Do not weaken the calculator test,
and do not add a calculator-specific production escape hatch.

The refactor is now developed on two independent axes:

```text
horizontal axis: implementation / migration phases
vertical axis:   calculator tower levels 0 → 9
```

A migration phase is not complete merely because the repository's historical
tests remain green. It must also satisfy every calculator level whose
mathematical capability it claims to implement or modify.

## Interlock law

Before production work at tower level \(n\):

1. all calculator levels \(0,\ldots,n-1\) must already be GREEN;
2. write the smallest Level-\(n\) calculator test that is RED for the missing
   capability;
3. implement the smallest general framework capability that makes it GREEN;
4. rerun all calculator levels \(0,\ldots,n\);
5. rerun all existing repository suites;
6. run the representative performance gate if a hot path changed;
7. report for architectural review before advancing the vertical frontier.

Existing production code above the current calculator frontier may remain in
place during migration. It does **not** count as proof that the skipped
calculator level exists.

If a lower calculator level exposes a missing abstraction, that lower-level
gap takes priority over further higher-level migration.

\[
\boxed{
\text{the calculator frontier controls architectural advancement}
}
\]

## Current migration consequence

The repository already contains substantial implementation corresponding to
carriers, relations, relational graphs, ordinary graph profiles, and support
subobjects.

Therefore the calculator tower must now be backfilled from Level 0 upward,
using the existing implementation wherever it genuinely satisfies the
calculator truth.

Expected behavior:

```text
Level 0  carriers             write test → should already be GREEN
Level 1  finite-arity relation write test → should already be GREEN
Level 2  relation algebra      write test → may expose the next real gap
Level 3  relational graph      only certify after Level 2 is GREEN
Level 4  graph calculus        only certify after Level 3 is GREEN
Level 5  field calculus        write RED test before invasive field migration
```

In particular, do **not** continue an invasive Level-5 field/support migration
past its RED test while Level 2 or Level 4 remains architecturally absent.

The calculator is now a real caller. A relation-algebra primitive previously
deferred until "a real caller requires it" is earned when the calculator's
Level-2 dependency derivation requires it.

## Phase/level mapping rule

Migration phase numbers and calculator level numbers are not required to be
identical.

For every migration phase, explicitly state:

```text
migration phase:
tower levels touched:
highest calculator level currently GREEN:
calculator tests added/changed:
production abstractions added/changed:
```

The mathematical tower order remains strict even when implementation migration
phases are arranged differently.

## No retroactive testing

Do not implement a migration phase and then invent a calculator test that
merely describes what the implementation happened to do.

The order is:

```text
architect-owned calculator truth
        ↓
RED test
        ↓
general production capability
        ↓
GREEN
        ↓
refactor
```

For capabilities already present before this TDD anchor was introduced, write
the test first and treat an unexpected failure as architectural evidence.

## Persistent-object law across both axes

The horizontal refactor and vertical calculator tower must operate on the same
mathematical object.

Do not create separate calculator representations for different phases.

The persistent structure remains:

\[
X=\{a,b,c,d,e\},\qquad
O=\{+,\times\},\qquad
P=\{in_1,in_2,out\},
\]

with

\[
R_{\mathrm{flow}}\subseteq O\times X\times P.
\]

Numerical field values first appear at Level 5.
Arithmetic meaning first appears at Level 8.

A migration implementation may change storage, views, indexing, or ownership,
but it must not change those calculator truths.

## Review gate for every migration checkpoint

Every agent report must now contain both views:

```text
MIGRATION
├── phase completed
├── production changes
├── existing repository suites
└── performance gate, if touched

CALCULATOR TOWER
├── highest GREEN level
├── new RED → GREEN truth
├── all lower levels still GREEN
├── forbidden higher-level imports: none
└── architectural observation / conflict
```

A phase affecting a calculator-relevant capability is not approved for the
next phase until this vertical report is present.


# 1. Authority and responsibility

## 1.1 The README is the narrative specification

`test/calculator-tower/README.md` is the canonical human-facing explanation.

Do not replace the calculator problem with another example.

Do not introduce a different calculator expression at different levels.

Do not simplify the tower by skipping levels.

Do not make numerical values structural objects merely to make a test easier.

The persistent expression is:

\[
(2+3)\times4.
\]

The persistent symbolic domains are conceptually:

\[
X=\{a,b,c,d,e\},
\]

\[
O=\{+,\times\},
\]

\[
P=\{in_1,in_2,out\}.
\]

The persistent numerical field is eventually:

\[
q(a)=2,\quad q(b)=3,\quad q(c)=5,\quad q(d)=4,\quad q(e)=20.
\]

---

## 1.2 The acceptance truths are architect-owned

The truths in this file are **not implementation suggestions**.

They are acceptance criteria.

If the implementation disagrees with a test truth:

```text
implementation changes
```

not:

```text
truth changes
```

Do not weaken, delete, bypass, mock away, or reinterpret a failing acceptance test merely to make the suite green.

If a truth appears inconsistent with the current framework, stop that implementation path and identify the architectural conflict.

A test may be changed only when the mathematical contract itself has explicitly changed.

---

## 1.3 The coding agent owns implementation

The agent should:

- create the calculator-tower test suites;
- write the smallest framework changes required for each test;
- preserve all existing suites;
- benchmark hot paths when the calculator work touches them;
- keep the README synchronized with actual demonstrated behavior;
- report architectural conflicts instead of silently routing around them.

The architect owns the oracle.

The agent owns making the implementation satisfy the oracle.

---

# 2. TDD discipline

For every level use:

```text
RED
  write the smallest test expressing the new level truth

GREEN
  implement only enough framework capability to satisfy it

REFACTOR
  remove duplication and improve representation without changing truth
```

Do not implement the whole calculator tower first and write tests afterward.

Do not pre-build abstractions because they "will probably be useful later."

Each level earns its implementation from a concrete calculator requirement.

---

# 3. The vertical dependency law

The calculator tower is also a dependency test.

A test at level \(n\) may import only modules from levels

\[
0,\ldots,n.
\]

It must not import a higher level.

The strongest acceptance rule is:

\[
\boxed{\text{level }n\text{ must not need level }n+1\text{ to state its truth}}
\]

If a lower-level calculator test requires a higher-level module, the tower layering is wrong.

Do not solve this by moving the test upward.

Repair the dependency boundary.

---

# 4. Persistent calculator fixture

Use symbolic **slots**, not numerical values, as structure.

Conceptually:

```text
X : value slots
    a  b  c  d  e

O : operations
    +  ×

P : ports
    in₁  in₂  out
```

The structural flow is:

```text
a ──in₁──▶ +
b ──in₂──▶ +
+ ──out──▶ c

c ──in₁──▶ ×
d ──in₂──▶ ×
× ──out──▶ e
```

The ternary relation is:

\[
R_{\mathrm{flow}}\subseteq O\times X\times P
\]

with exactly the six tuples

\[
\begin{aligned}
(+,a,in_1),\qquad &(+,b,in_2),\qquad (+,c,out),\\
(\times,c,in_1),\qquad &(\times,d,in_2),\qquad (\times,e,out).
\end{aligned}
\]

Do not encode these six facts as an ordinary vertex-edge graph at Level 1.

The ternary relation is the point of the test.

---

# 5. Recommended test tree

Create:

```text
test/calculator-tower/
├── README.md
├── run.sh
├── common/
│   └── calculator_assert.f90
├── level-0-carrier/
│   └── test.f90
├── level-1-relation/
│   └── test.f90
├── level-2-relation-algebra/
│   └── test.f90
├── level-3-graph/
│   └── test.f90
├── level-4-graph-calculus/
│   └── test.f90
├── level-5-field-calculus/
│   └── test.f90
├── level-6-discretization/
│   └── test.f90
├── level-7-minimization/
│   └── test.f90
├── level-8-constitution/
│   └── test.f90
└── level-9-statement/
    └── test.f90
```

`common/` may contain assertion helpers and dependency-free symbolic constants.

It must not become a back door through which higher-level framework modules enter lower-level tests.

Prefer tiny independent fixtures over a large magical shared fixture.

Some repetition in tests is acceptable when it makes the level boundary obvious.

---

# 6. Test output

The top-level runner should produce a compact vertical summary:

```text
calculator tower
├── level 0  carriers .............. PASS
├── level 1  relation .............. PASS
├── level 2  relation algebra ...... PASS
├── level 3  relational graph ...... PASS
├── level 4  graph calculus ........ PASS
├── level 5  field calculus ........ PASS
├── level 6  discretization ........ PASS
├── level 7  minimization .......... PASS
├── level 8  constitution .......... PASS
├── level 9  statement ............. PASS
└── mathematical result ............ 20
```

The physical Casio fx-115ES is a **manual external oracle**, not a CI dependency.

Do not attempt to automate the physical calculator.

The CI/software truth is:

\[
\boxed{20}
\]

and the README instructs the human to verify the same result independently on the Casio.

---

# 7. Level 0 — carriers

## New mathematical commitment

Level 0 answers:

> What members may exist?

Create independent carriers corresponding to:

\[
X,\quad O,\quad P.
\]

No relation exists yet.

No graph exists yet.

No arithmetic law exists yet.

No field values exist yet.

---

## Required truths

For every member \(x\) of every carrier:

\[
member(local\_index(x))=x.
\]

For every valid local index \(i\):

\[
local\_index(member(i))=i.
\]

Verify:

```text
X.same_as(O) = false
X.same_as(P) = false
O.same_as(P) = false
```

Verify exact cardinalities:

```text
|X| = 5
|O| = 2
|P| = 3
```

Verify an outsider is rejected by membership:

```text
X.has(outsider) = false
```

If a non-counted/listed carrier is used, verify carrier enumeration remains injective.

---

## Required negative truth

The test must not require:

```text
vertex
edge
graph
relation
field
residual
solver
```

to express Level 0.

---

## Acceptance meaning

Passing Level 0 proves:

\[
\boxed{\text{independent typed member domains exist}}
\]

without graph interpretation.

---

# 8. Level 1 — generic finite-arity relation

## New mathematical commitment

Level 1 answers:

> How may members of carriers be related?

Build:

\[
R_{\mathrm{flow}}\subseteq O\times X\times P.
\]

Its arity is exactly 3.

---

## Required truths

The relation contains exactly:

```text
(+, a, in₁)
(+, b, in₂)
(+, c, out)
(×, c, in₁)
(×, d, in₂)
(×, e, out)
```

Representative assertions:

```text
R.has(+, a, in₁)  = true
R.has(+, c, out)  = true
R.has(×, c, in₁)  = true
R.has(×, a, in₁)  = false
```

Signature truths:

```text
arity(R) = 3
domain(R,1) same_as O
domain(R,2) same_as X
domain(R,3) same_as P
```

Set semantics:

```text
adding an existing tuple again does not increase |R|
```

Foreign-domain members must be refused.

---

## Required negative truth

Do not reduce the relation to:

```text
ordinary graph
binary incidence
vertex-edge pair
```

merely to satisfy the test.

The ternary relationship is deliberate.

---

## Acceptance meaning

Passing Level 1 proves:

\[
\boxed{\text{the framework natively expresses a meaningful }k=3\text{ relation}}
\]

---

# 9. Level 2 — relation algebra

## New mathematical commitment

Level 2 answers:

> How can existing relations generate new relations?

Derive the operation-dependency relation

\[
D\subseteq O\times O.
\]

An operation \(o_1\) depends into \(o_2\) when an output value slot of \(o_1\) is an input value slot of \(o_2\).

For the calculator:

\[
\boxed{D=\{(+,\times)\}}.
\]

---

## Required truths

Exactly:

```text
D.has(+, ×) = true
D.has(×, +) = false
D.has(+, +) = false
D.has(×, ×) = false
|D| = 1
```

The result must be invariant to enumeration order of the six input tuples.

Test this explicitly by building the same \(R_{\mathrm{flow}}\) with a different tuple ordering and proving the same derived extension.

---

## Algebraic construction

Use the smallest algebraic operations actually required.

Conceptually the construction may be read as:

```text
restrict flow to output tuples
restrict flow to input tuples
join on X
project onto O × O
```

Do not add a general algebra API merely because the notation exists.

Each algebraic primitive must earn its place through this or another real caller.

---

## Required negative truth

Do not store \(D\) manually as a hard-coded calculator special case.

The test is specifically for **derivation**.

---

## Acceptance meaning

Passing Level 2 proves:

\[
\boxed{\text{relations can generate relations}}
\]

---

# 10. Level 3 — relational graph

## New mathematical commitment

Level 3 answers:

> How do carriers and relations coexist as one structure?

Construct:

\[
G=(\mathcal S,\mathcal R)
\]

with at least

\[
\mathcal S=\{X,O,P\}
\]

and

\[
\mathcal R\ni R_{\mathrm{flow}}.
\]

If Level 2 has materialized \(D\), it may also be admitted when appropriate.

---

## Required truths

Verify:

```text
G owns X
G owns O
G owns P
G owns R_flow
```

Every relation domain must resolve to a carrier owned by \(G\).

Two distinct relations with the same signature must be allowed to coexist.

Relations are addressed by identity, not signature alone.

A foreign carrier in a relation signature must be refused.

A borrowing/non-materialized relation view must not be accepted as graph-owned stable storage unless the graph ownership contract has explicitly evolved to make that safe.

---

## Most important negative truth

The Level-3 calculator test must compile without importing an ordinary-graph profile.

There must be no requirement for:

```text
vertex
edge
tail
head
```

in `relational_graph`.

---

## Acceptance meaning

Passing Level 3 proves:

\[
\boxed{\text{Graph is a collection of member sets and relations}}
\]

rather than:

\[
\boxed{\text{Graph means }(V,E)}
\]

---

# 11. Level 4 — graph calculus

## New mathematical commitment

Level 4 answers:

> What graph-theoretic interpretations and algorithms may be applied?

Interpret

\[
D=\{(+,\times)\}
\]

as a directed dependency graph over operation carrier \(O\).

---

## Required truths

Exactly:

```text
sources(D) = {+}
sinks(D) = {×}
reachable(+, ×) = true
reachable(×, +) = false
```

A valid topological walk must be:

```text
[+, ×]
```

No other operation order is permitted.

---

## Ordering law

If graph-calculus output has a canonical order, canonicalize by carrier enumeration/local index, not by raw integer member value.

The meaning of a member is not its numeric storage representation.

---

## Required negative truth

Do not teach the generic graph container itself about topological sorting.

Traversal/calculus acts on graph structure.

It is not graph storage.

---

## Acceptance meaning

Passing Level 4 proves:

\[
\boxed{\text{graph interpretation produces graph algorithms}}
\]

without contaminating the lower relational ontology.

---

# 12. Level 5 — field calculus and supports

## New mathematical commitment

Level 5 answers:

> What values live on a domain or subdomain?

The calculator introduces the first numerical values.

Define:

\[
q:X\to\mathbb R.
\]

Known input support:

\[
K=\{a,b,d\}\hookrightarrow X.
\]

Unknown/computed support:

\[
U=\{c,e\}\hookrightarrow X.
\]

Known values:

\[
q(a)=2,\qquad q(b)=3,\qquad q(d)=4.
\]

---

## Support law

A support is an indexable subobject:

\[
S\hookrightarrow A.
\]

It should behave as a `member_set` refinement/subdomain with an ambient carrier.

Required law:

\[
s\in S\implies s\in A.
\]

A field always lives on an indexable member domain.

Do not create a `field` API that must branch between unrelated "carrier" and "predicate" domain kinds.

---

## Required truths

For \(K\):

```text
|K| = 3
K.has(a) = true
K.has(b) = true
K.has(d) = true
K.has(c) = false
K.has(e) = false
```

For \(U\):

```text
|U| = 2
U.has(c) = true
U.has(e) = true
```

Round-trip indexing must hold on supports:

\[
member(local\_index(s))=s.
\]

Field truths:

```text
q(a) = 2
q(b) = 3
q(d) = 4
```

The test should demonstrate that field storage follows domain enumeration, not assumptions about contiguous/raw member identifiers.

---

## Required negative truth

Support must not pretend to be an edgeless graph.

Do not preserve obsolete support APIs merely because old characterization tests once required them.

Rewrite the marked transitional support tests to the new subobject law.

---

## Acceptance meaning

Passing Level 5 proves:

\[
\boxed{\text{values are functions over domains; supports are true subdomains}}
\]

---

# 13. Level 6 — discretization

## New mathematical commitment

Level 6 answers:

> How does relational structure become a discrete equation dependency pattern?

Introduce residual locations:

\[
r_c,\qquad r_e.
\]

The calculator's structural residual dependencies are:

\[
r_c=r_c(q_a,q_b,q_c),
\]

\[
r_e=r_e(q_c,q_d,q_e).
\]

---

## Required truths

Exactly:

```text
support(r_c) = {a,b,c}
support(r_e) = {c,d,e}
```

No other value slot may appear in either dependency set.

Equivalent structural Jacobian pattern:

```text
        a  b  c  d  e
r_c     ×  ×  ×
r_e           ×  ×  ×
```

---

## Important separation

At Level 6 it is enough to know **which quantities participate**.

Do not require arithmetic `+` or `×` semantics to state the sparsity/dependency truth.

The calculator structure must determine the dependency layout before the constitution determines coefficient/action meaning.

---

## Tangent/adjoint readiness

Where the framework exposes derivative structure, the same dependency representation should be usable by forward and reverse traversals.

Do not add a second independent dependency description just for adjoints.

---

## Acceptance meaning

Passing Level 6 proves:

\[
\boxed{\text{relations + domains compile into discrete operator structure}}
\]

---

# 14. Level 7 — minimization

## New mathematical commitment

Level 7 answers:

> Can a residual equation be solved?

Use the simplest residual compatible with the expected calculator solution:

\[
r_c=q(c)-5,
\]

\[
r_e=q(e)-4q(c).
\]

At this level these equations are supplied.

Their arithmetic origin is not yet the test.

---

## Required truths

Start from deliberately incorrect unknown values.

For example:

```text
q(c) = 0
q(e) = 0
```

Drive:

\[
r(q)\to0.
\]

Require:

\[
\|r\|\le\epsilon.
\]

Require:

\[
q(c)=5,
\]

\[
q(e)=20,
\]

within a tolerance appropriate to the real kind used by the framework.

---

## Required negative truth

Do not special-case the solver for the calculator.

The test must use the normal minimization/residual machinery.

The solver must not know that `+` means addition or `×` means multiplication.

---

## Acceptance meaning

Passing Level 7 proves:

\[
\boxed{\text{the solver solves equations independently of their physical meaning}}
\]

---

# 15. Level 8 — constitution

## New mathematical commitment

Level 8 answers:

> What laws do the operation symbols obey?

Only here assign arithmetic meaning:

\[
+(x,y)=x+y,
\]

\[
\times(x,y)=xy.
\]

From the Level-1 structure and Level-5 field values, generate:

\[
r_c
=
q(c)-\left(q(a)+q(b)\right),
\]

\[
r_e
=
q(e)-q(c)q(d).
\]

Substitution of known values yields:

\[
r_c=q(c)-5,
\]

\[
r_e=q(e)-4q(c).
\]

These must match the residual system used in Level 7.

---

## Required truths

At

\[
q(a)=2,\ q(b)=3,\ q(c)=5,\ q(d)=4,\ q(e)=20,
\]

verify:

\[
r_c=0,
\]

\[
r_e=0.
\]

Also verify the constitution-generated dependency support remains exactly the Level-6 support.

The constitution may add numerical meaning.

It must not alter topology.

---

## Required negative truth

Do not move arithmetic meaning into:

```text
member_set
relation
relational_graph
graph calculus
field
```

The symbol `+` is only a member of \(O\) below constitution.

---

## Acceptance meaning

Passing Level 8 proves:

\[
\boxed{\text{constitutive/model laws bind meaning to pre-existing structure}}
\]

---

# 16. Level 9 — statement

## New mathematical commitment

Level 9 answers:

> What problem is being asked?

Construct the complete calculator problem:

\[
\boxed{\text{Evaluate }(2+3)\times4.}
\]

The statement supplies or selects:

```text
structure      calculator relational graph
inputs         a=2, b=3, d=4
constitution   addition and multiplication
requested      output e
```

---

## Required truth

The framework must produce:

\[
\boxed{20}.
\]

The software result must not be injected as a hard-coded expected result into the evaluation path.

It may of course appear in the test assertion.

---

## Manual independent oracle

The README should instruct the human to type on the Casio fx-115ES:

```text
( 2 + 3 ) × 4 =
```

and observe:

```text
20
```

The physical calculator is independent confirmation that the entire vertical software composition preserves the intended mathematical statement.

---

## Acceptance meaning

Passing Level 9 proves:

\[
\boxed{\text{the full tower composes to the same mathematical truth}}
\]

---

# 17. Cross-level invariants

These are tested wherever practical.

## 17.1 Persistence

The calculator does not become a different object at each level.

The meanings accumulate.

```text
Level 0: a exists
Level 1: a participates in a relation
Level 3: that relation belongs to a graph
Level 5: a carries q(a)=2
Level 6: a participates in r_c
Level 8: its role is interpreted by addition
Level 9: it is the literal 2 in the requested expression
```

---

## 17.2 Lower truths remain true

A higher level must not invalidate a lower-level law.

Examples:

- constitution cannot change carrier identity;
- minimization cannot change graph topology;
- field values cannot alter relation membership;
- graph calculus cannot mutate relation extension;
- the statement cannot redefine support membership.

---

## 17.3 Structure before meaning

The framework must be able to describe:

```text
+ consumes a and b and produces c
```

before it knows:

```text
+ means numerical addition
```

Likewise:

```text
× consumes c and d and produces e
```

before multiplication is attached.

This is one of the most important calculator-tower proofs.

---

## 17.4 No duplicate sources of truth

Do not simultaneously store:

```text
flow relation
and a separately maintained dependency graph
and separately maintained residual dependencies
```

when one can be generated from another under the level contracts.

Derived structures may be materialized for performance only when their derivation and consistency laws are explicit.

---

## 17.5 Set semantics

Carriers and relations obey mathematical set semantics.

Tuple/member insertion order must not alter extensional truth.

Where APIs expose canonical order, use domain enumeration/local indices.

Do not infer semantic order from raw integer values.

---

# 18. Failure-driven architecture

A calculator test failure is architectural information.

When a level cannot be expressed cleanly, ask:

```text
What mathematical statement is missing from the level below?
```

Do not first ask:

```text
What special API can make this calculator test pass?
```

Examples:

- if Level 2 cannot derive dependency, relation algebra is missing a real primitive;
- if Level 5 cannot index a support, the support is not yet a proper domain;
- if Level 6 must know addition to determine sparsity, discretization and constitution are entangled;
- if Level 4 must import ordinary vertex/edge grammar, the relational abstraction leaked;
- if Level 9 has to manually construct all lower-level details, the higher-level composition contract is incomplete.

Use the calculator to expose these seams.

---

# 19. Do not overfit to the calculator

The calculator is deliberately tiny.

Its purpose is to expose the tower, not define the tower.

Do not introduce framework types such as:

```text
calculator_graph
addition_relation
multiplication_node
calculator_field
```

unless such a type is independently justified by a mathematical contract beyond this demo.

The calculator should instantiate general abstractions.

It should not become one.

---

# 20. No fake success

Forbidden:

- hard-coding `20` in production code;
- bypassing a level because the expression is easy to evaluate directly;
- constructing dependency \(+\to\times\) manually in the relation-algebra test;
- putting arithmetic semantics into the Level-1 relation;
- letting the Level-7 solver call the Level-8 calculator evaluator;
- changing test expectations to current implementation behavior without architectural review;
- marking a level PASS when the test is skipped;
- compiling all tests against all modules and thereby hiding forbidden dependencies.

---

# 21. Build strategy

## Existing lower levels

For levels already implemented by the current branch:

1. write the calculator test;
2. expect it to pass with little or no production change;
3. treat any failure as a useful mismatch between claimed capability and actual generality;
4. make the smallest correction;
5. run the full repository suite.

## Emerging higher levels

For levels not yet implemented:

1. write the first compiling RED test as soon as the lower API permits the truth to be stated;
2. add the minimal next-level contract;
3. make it GREEN;
4. refactor;
5. move upward one level.

Do not create dummy higher-level types solely so all ten test directories compile immediately.

A missing level should remain visibly unfinished until its real contract exists.

---

# 22. Test independence and shared fixtures

Prefer this balance:

```text
shared:
    assertion helpers
    dependency-free symbolic labels/constants

local to each level:
    framework object construction
    level-specific fixtures
    level-specific expected truths
```

Reason:

A giant common fixture can accidentally construct Level-8 objects and hand them to Level-2 tests, making the layering test meaningless.

The small duplication in rebuilding a five-slot calculator is worth the clarity.

---

# 23. Refusal tests

Add refusal tests when a level has a meaningful invalid construction.

At minimum consider:

```text
Level 0
    duplicate/invalid carrier declaration as defined by carrier contract

Level 1
    wrong arity
    foreign member
    undeclared domain

Level 3
    foreign relation domain
    duplicate owned identities
    non-materialized borrowing relation ownership

Level 5
    support member outside ambient carrier

Level 8
    constitution applied to incompatible operation/domain schema

Level 9
    statement missing a required input or requested output
```

Do not manufacture refusal tests for conditions that the mathematical level does not promise.

---

# 24. Performance

The calculator itself is too small to benchmark meaningfully.

Do not use it to claim hot-path performance.

When calculator-driven implementation changes touch:

```text
carrier lookup
relation fibres
graph traversal
field indexing
discrete operator application
```

run the existing representative benchmarks as regression gates.

The calculator validates semantics.

Large graph/mesh benchmarks validate performance.

Keep those roles separate.

---

# 25. Documentation discipline

Every implemented calculator level should keep its README subsection truthful.

If the implementation uses a materially different but mathematically equivalent construction, update the detailed subsection while preserving:

- the persistent expression;
- the level capability;
- the architect-owned truth;
- the dependency boundary.

Do not let README diagrams become aspirational after a level is marked PASS.

---

# 26. Commit discipline

Prefer one coherent commit per completed RED→GREEN→REFACTOR level or tightly related architectural correction.

Each commit report should include:

```text
LEVEL N
├── RED truth introduced
├── production capability added or reused
├── GREEN result
├── existing suites
└── architectural observation, if any
```

Do not bundle unrelated refactors into calculator-tower work.

---

# 27. Review gates

Before beginning level \(n+1\), report level \(n\) for architectural review.

The report should answer:

1. What new truth was tested?
2. Did the test require any higher-level import?
3. What production API, if any, was added?
4. Did an existing abstraction fail to generalize?
5. Did any existing behavior change?
6. Are all previous calculator levels still green?
7. Are all repository suites still green?
8. Was any performance-sensitive path touched?

Do not proceed automatically through all ten levels in one unattended redesign.

The calculator tower is intended to keep the architecture anchored through incremental review.

---

# 28. Suggested implementation order

The **calculator acceptance frontier** advances only in tower order:

```text
0  carrier
1  relation
2  relation algebra
3  relational graph
4  graph calculus
5  field calculus
6  discretization
7  minimization
8  constitution
9  statement
```

This order is itself part of the experiment.

The repository's migration phases may be sequenced differently because they
must preserve and progressively replace existing production code. That does
not permit the calculator acceptance frontier to skip a mathematical level.

Thus:

```text
production may already contain higher-level machinery
        ≠
higher calculator level is certified
```

If Level \(n\) is missing or RED, do not certify Level \(n+1\), even if related
production code already exists.

If a lower calculator level reveals that an abstraction previously scheduled
for a later migration phase is actually required now, the calculator has
provided a real caller. Implement the smallest general capability at the
correct mathematical level, then resume the migration.

Do not silently reorder the mathematical tower.

---

# 29. Definition of done

The calculator tower is complete only when:

```text
[ ] all ten level tests exist
[ ] every level test imports only its level and below
[ ] every architect-owned truth passes
[ ] all refusal tests pass
[ ] all pre-existing repository suites pass
[ ] no calculator-specific production abstractions were introduced without justification
[ ] no duplicate source of structural truth was introduced
[ ] README diagrams match implemented behavior
[ ] top-level run.sh reports all levels
[ ] final software result is 20
[ ] README documents Casio fx-115ES manual oracle = 20
```

The deepest acceptance criterion is:

\[
\boxed{
\text{the same calculator object gains one new mathematical capability at each level}
}
\]

while:

\[
\boxed{
\text{every lower-level truth remains unchanged}
}
\]

That is what the calculator tower is testing.c
