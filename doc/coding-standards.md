# Mathematical naming

A name identifies the mathematical object, operation, domain or numerical
quantity represented by the code. Use the same name for the same concept
throughout its declaration and consumers. A name must remain meaningful
without a story about what an object knows, carries, answers or remembers.

## Objects and operations

- Name a set by its role: `retained_columns`, `active_members`,
  `fixed_indicator`, `halo_vertices`, `source_set` or `sink_set`.
- Name a numerical quantity by its definition: `initial_residual_norm`,
  `maximum_error`, `num_retained`, `iteration_limit` or `pivot_rows`.
- Name a procedure by its operation: `extend_directions`,
  `initialize_residual_history`, `record_residual_norm`, `solve_linear`,
  `linearize` or `restrict`.
- Name a logical value by the proposition it represents: `satisfied`,
  `tangent_defined`, `factors_current` or `converged`.
- Conventional mathematical symbols such as `q`, `r`, `lambda`, `i` and `j`
  are appropriate when their equations, domains or index sets are clear.
- Keep established mathematical and numerical vocabulary: Newton steps,
  Krylov bases, graph paths, trees, transposes and the `ceiling` intrinsic.
  A computer-specific object uses its exact technical role, such as
  `source_path`, `staged_tree`, `compiler` or `storage_size`.

Colloquial identifiers such as `kept`, `got`, `verdict`, `rooms`, `seat`,
`worst` and compounds containing them are refused. A retained column is
not a stored factorisation: replacing every occurrence with one generic
synonym would obscure two different objects. Read each declaration and
its uses before selecting the replacement.

Renames preserve equations, numerical constants and control flow. Update
all consumers, including keyword arguments and overriding Fortran
procedures, in the same change. Do not retain aliases for superseded names.

## Library structure

Library modules use namespace order with one of these prefixes:
`graph`, `token`, `relation`, `field`, `operation`, `transform`, `map`,
`view` or `util`. Types use mathematical adjective-noun order without a
module namespace prefix. Use `num_` for cardinalities and distinguish a
member's identity from its local index.

Every library source has an entry in `src/OBJECTS`. That list places each
module after its dependencies. Library readers use nouns rather than a
`get_` prefix. Library imports do not rename symbols, and a module does
not re-export an imported name unless it defines an extending generic.
The historical `dp => real64` kind import remains the checker's explicit
exception. Use British spelling, including `centre`, `colouring`,
`neighbour` and `fibre`.

The packed application follows the type, reader and spelling rules; its
multi-module source does not use the library's file-prefix rule.

## Enforcement and verification

Run `./check_naming.sh`. Its canonical vocabulary list is checked by
`tools/check_vocabulary.py` against every non-ignored Fortran source,
Python identifiers and shell declarations, including tests and build
scripts. Identifier components separated by underscores are checked;
Fortran case and continuations cannot conceal a refused word. Standard
language syntax and imported APIs retain the names required by their
language or dependency.

Production Fortran comments and diagnostic strings also follow the
closed vocabulary. Test descriptions and deliberate invalid-source
strings are not declarations. This permits refusal tests to express the
invalid input they are required to reject.

`python3 test/naming/test.py` verifies the checker, including a staged
invalid file concealed by a corrected working copy. `githooks/pre-commit`
checks an export of the staged tree on every commit. Install it with
`git config core.hooksPath githooks`; the current repository uses that
hook directory. `./verify.sh` runs naming checks and their regressions
before building and running the numerical suites.

The checker enforces the explicit lexical rules. Review must still verify
that each name denotes the object's actual mathematical role; a finite
word list cannot establish that correspondence by itself.
