# 69. The laws that make code read as notation

Section 0 made naming a correctness law.  Naming alone is not sufficient.  A
reader with no repository history must be able to recover the mathematics from
the source alone, and a name can be canonical while the surrounding structure
still hides which object decided what.

This section states the remaining six obligations in the same form as Section
0: a law, the normal form or mechanical test that decides it, and the verdict
vocabulary that must appear in the handoff.

Each subsection answers one question:

```text
69.1  names        does one word denote one definition?
69.2  roles        does each argument declare what it is for?
69.3  authority    does every structural fact have a declared supplier?
69.4  axes         does one verb keep one meaning across concretions?
69.5  abstraction  does this object exist for its own mathematics?
69.6  ownership    does exactly one object sign this structure?
69.7  verification does a checked artifact fail when the claim stops holding?
```

A change is complete when all seven produce a verdict.  Silence is not a pass
in any of them.

---

## 69.1 Names — the audit extends to declared domains

Section 0 governs identifiers.  It does not yet govern the string a domain is
declared with, and that string is a name under exactly the same law:

```text
index_set('vertices', nv)
```

`'vertices'` is a public mathematical word.  Apply Section 0 to it.

Two failures are possible and they are not symmetric.

**Two labels, one identity.**  A violation.  If `cells` and `faces` are only
other words for `vertices` and `edges` — the same tokens, obtained by aliasing
a graph's own sets — then one distinction has two words, and the inverse
map from code to mathematics is destroyed at the domain level.  Either the
domains are genuinely distinct and must be declared, or the alias must go.

**One label, two identities.**  Admissible, but only when stated.  Transport
re-declares a domain under its source's name deliberately:

```text
sp = subset(dom % name(), part_set, kept(1:n))
```

Identity is not preserved across a partition boundary; extension and values
are.  That is a law, not an accident, so it is written where it happens and
tied to the round-trip acceptance law.  An unstated repetition of a label is a
naming defect.

Verdict: `canonical` / `alias` / `restated-by-law`.

---

## 69.2 Roles — a graph argument declares what it is for

A `class(graph)` dummy argument is not self-describing.  The same declared type
serves several mathematically distinct purposes, and they are a closed list:

```text
operand      the computation reads its structure arithmetically
context      the action executes over it and asks nothing structural of it
authority    it supplies a structural fact the routine is not allowed to infer
owner        it constructs and signs the structure
view         it reinterprets another graph's structure and owns none
```

Normal form, written before the argument is accepted:

```text
<name> : graph = <role> for <what>
```

For the minimization boundary:

```text
on       : graph = execution context for the action
coupling : graph = authority for the dependent axis
```

The test: if a procedure takes two `class(graph)` arguments and you cannot
write both lines without a slash or a near-synonym, the two are one role
duplicated or two roles conflated.  If a single graph argument needs two lines,
it is carrying two roles and must be split.

The law that makes the axes orthogonal rather than merely listed:

\[
\boxed{
\text{a role does not become another role by proximity}
}
\]

In particular, **context does not become authority**.  The graph an action
happens to execute over is not evidence about anything structural.  A routine
that reads structure from its execution context has silently collapsed two
axes that Section 12 of the fractal document requires to stay independent.

Verdict, per graph argument: one role word, or `conflated`.

---

## 69.3 Authority — every structural fact has a declared supplier

Enumerate the structural questions a routine asks.  Not the values it computes
— the facts it takes as given:

```text
which unknowns are coupled?
which domain does the answer live on?
which instants does the residual read?
which members does this part own?
```

For each one, name the seat that answers it.  If no seat exists, the routine is
inferring, and inference is the defect.  Three resolutions are admissible:

1. **give it a seat, with no fallback** — the fact arrives explicitly or not
   at all;
2. **refuse loudly** — a routine that needed the fact and was handed nothing
   says so, and stops;
3. **make the caller state the equality** — where two structures really are
   the same object, the caller asserts it at its own call site, where it is a
   claim that can be wrong rather than an assumption that cannot be seen.

Defaulting is not among them.  A seat that quietly falls back to a nearby
object is worse than no seat, because it converts a missing fact into a
plausible wrong one.

The pattern to refuse, written generally:

```text
if (.not. allocated(this % authority)) this % authority = this % context
```

That line is the entire defect.  It reads as a convenience and it is a
statement that context and authority are the same axis.

Verdict, per structural fact: `declared` / `refused` / `claimed-by-caller` /
`inferred`.  Only the first three pass.

---

## 69.4 Axes — one verb keeps one meaning across concretions

A deferred binding whose meaning depends on which axis it is asked about must
resolve that axis from the concrete type, never from the surrounding
application.

**The two-implementation test.**  Write the one-sentence answer for two
different concretions of the same binding:

```text
step_operator    % dependencies : the stencil on the INDEPENDENT axis
stencil_operator % dependencies : the stencil on the DEPENDENT   axis
```

If the two sentences differ in anything but the axis word, the contract has two
meanings and at least one concretion is answering the wrong question.  If they
differ only in the axis word, the contract is one contract.

This test catches a specific and durable error: answering a *true* relation
that is not the *owed* relation.  A step's succession

```text
1 -> 2 -> 3
```

and a step's stencil

```text
1 -> 3,  2 -> 3,  3 -> 3
```

are both real structure over the same instants.  Only the second says which
instants the residual reads, and only the second carries the self-arrow that
makes the newest instant an unknown rather than data.

\[
\boxed{
\text{a stencil is not a chronology}
}
\]

Second law of this subsection: **the axis is not the application.**  An
independent axis need not be time; a continuation coordinate or a parameter
takes the same seat.  Do not name the verb, the type, or the argument after
whichever application first needed it.

Verdict: `one meaning` / `axis inferred from application` / `two meanings`.

---

## 69.5 Abstraction — the object exists for its own mathematics

Section 50 gives five clauses by which a type earns existence.  It is missing
the test that catches the commonest modern failure, which is not a type
inflated from an English noun but an object fabricated to satisfy another
object's constructor.

**The fabrication test.**  For every declared object ask:

> Does this exist because the mathematics contains it, or because some other
> type would not compile without an argument of its shape?

An object that appears only as an argument, is never queried on its own terms,
and carries no law of its own, is a fixture wearing a production type.

The worked example is a subobject constructor that admits no arbitrary roster.
When the only way to declare a domain with chosen members is

\[
S\hookrightarrow A
\]

then a domain that has no natural host must invent one, and the invented host
exists solely to be embedded in.  Its size is chosen to equal its subset's
size; nothing ever asks it a question.  That object fails all five clauses of
Section 50 and the fabrication test names why.

Two admissible repairs, both architectural decisions under 0.3: admit a second
set concretion that holds an explicit roster, or accept the host and give
it a mathematical reason to exist.  Silence is the third option and it is not
admissible.

Verdict: `earned` / `fabricated` / `taxonomy review`.

---

## 69.6 Ownership — exactly one object signs each structure

Identity is minted once, at declaration, and a second signing is refused
loudly.  That law is enforced for sets, relations, and graphs.  Extend the
same discipline upward to every structure a routine can construct.

**The signing test.**  For one structure, count the objects that can bring it
into existence.  If more than one can, ownership is split, and the two will
diverge under partition, transport, or time.

Distinguish two relationships that look alike at a call site:

```text
owner       constructs the structure and signs it
passenger   holds a reference and may read it, and never signs
```

A clock that advances over a graph is a passenger.  A solver attached to a
mesh is a passenger.  Neither is entitled to declare a domain on the other's
behalf, and neither becomes the owner by being the object the user called.

Construction complexity belongs to builders.  A builder may mutate; the
structure it produces may not.  If a structural object exposes a mutator after
construction, ownership has leaked into its lifetime.

Verdict, per structure: `single owner` / `passenger` / `split`.

---

## 69.7 Verification — a claim is checked, not asserted

A claim written in prose is not a claim.  A claim is an artifact that fails
when the claim stops being true.

**The procedure.**  For every architectural claim a tower or a review makes:

1. name the artifact that fails when the claim stops holding;
2. check the claim itself, not a proxy for it;
3. plant a violation and confirm the artifact fails;
4. record the claim's history, including readings that were withdrawn.

Step 2 is where most verification is lost.  Refusing an import does not
establish that a level uses no numbers — a source can declare its own

```text
real(dp) :: w = 2.0_dp
```

and help itself.  The claim is established only by reading the sources for
what the claim actually forbids.  A test of a neighbouring fact is a proxy,
and a proxy that has never been distinguished from the claim is not evidence.

Step 3 is where the rest is lost.  A checker that has never failed is not
known to be able to fail.  Selftest it against a planted violation, and against
the near-misses that must still pass — an integer set size is not a
coefficient, and a comment containing the word "real" is prose.

Step 4 is the evidence trail.  When a reading turns out to be wrong, withdraw
it in the record rather than quietly correcting it.  A tower's observations are
worth what their failures are worth.

The same rule governs delegated work: a worker's reported `PASS` is a claim,
and it requires independent reproduction before it becomes a verdict.

Verdict, per claim: `checked` / `checked by proxy` / `asserted`.

---

## 69.8 Acceptance

Every implementation or review handoff carries one line per axis:

```text
Names       : canonical | <label -> verdict>
Roles       : <graph argument -> role> ...
Authority   : <structural fact -> declared | refused | claimed-by-caller>
Axes        : one meaning | <binding -> defect>
Abstraction : earned | <object -> fabricated>
Ownership   : single owner | <structure -> split>
Verification: checked | <claim -> proxy | asserted>
```

An axis that produced no finding says so.  An axis that produced a finding
outside the current change reports it and does not expand the patch.

The obligation these seven share is one sentence:

\[
\boxed{
\text{the source states the mathematics; nothing else is consulted to recover it}
}
\]
