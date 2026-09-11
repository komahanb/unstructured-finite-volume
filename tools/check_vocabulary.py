#!/usr/bin/env python3
"""Check mathematical vocabulary in source identifiers and production prose."""

import ast
import os
from pathlib import Path
import re
import subprocess
import sys


IDENTIFIER = re.compile(r"[A-Za-z][A-Za-z_0-9]*")
FORTRAN_SUFFIXES = {".f", ".for", ".ftn", ".f90", ".f95", ".f03", ".f08"}


def source_paths(root):
    """Use the same files as Git; staged exports contain only committed paths."""
    repository = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "--show-toplevel"],
        capture_output=True, text=True, check=False,
    )
    if repository.returncode == 0 and Path(repository.stdout.strip()).resolve() == root:
        output = subprocess.check_output(
            ["git", "-C", str(root), "ls-files", "-z", "--cached", "--others", "--exclude-standard"]
        )
        paths = {root / os.fsdecode(path) for path in output.split(b"\0") if path}
    else:
        paths = {path for path in root.rglob("*") if path.is_file() and ".git" not in path.parts}
    return sorted(path for path in paths if path.is_file() and (
        path.suffix.lower() in FORTRAN_SUFFIXES | {".sh", ".py"}
        or path.relative_to(root).as_posix() == "githooks/pre-commit"
    ))


def identifier_components(identifier):
    # Fortran case never changes identity. Camel-case components supplement,
    # rather than replace, the components of the case-folded identifier.
    components = set(identifier.lower().split("_"))
    separated = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", identifier)
    separated = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", separated)
    components.update(separated.lower().split("_"))
    return components


def fortran_code(source, fixed_form=False):
    """Mask literals and comments without changing source offsets."""
    characters = list(source)
    quote = None
    offset = 0
    for line in source.splitlines(keepends=True):
        if ((not fixed_form and line.lstrip().startswith("!"))
                or (fixed_form and line and line[0] in "cC*!")):
            characters[offset:offset + len(line.rstrip("\r\n"))] = " " * len(line.rstrip("\r\n"))
            offset += len(line)
            continue
        column = 6 if fixed_form else 0
        while column < len(line):
            char = line[column]
            if quote:
                if char not in "\r\n":
                    characters[offset + column] = " "
                if char == quote:
                    if column + 1 < len(line) and line[column + 1] == quote:
                        characters[offset + column + 1] = " "
                        column += 1
                    else:
                        quote = None
            elif char in "'\"":
                quote = char
                characters[offset + column] = " "
            elif char == "!":
                end = len(line.rstrip("\r\n"))
                characters[offset + column:offset + end] = " " * (end - column)
                break
            column += 1
        offset += len(line)
    return "".join(characters)


def is_mathematical_intrinsic(code, start, end):
    """CEILING is an exact mathematical intrinsic, not a user-defined label."""
    if code[start:end].lower() != "ceiling":
        return False
    prefix = code[:start].rsplit(";", 1)[-1]
    if re.fullmatch(r"\s*intrinsic\s*(?:::)?\s*(?:[a-z]\w*\s*,\s*)*", prefix, re.I):
        return True
    suffix = code[end:].lstrip()
    if (not suffix.startswith("(") or prefix.rstrip().endswith("%")
            or re.search(r"\b(function|subroutine|interface|call)\s*$", prefix, re.I)):
        return False
    # Reject an entity named CEILING, while permitting the intrinsic in an
    # entity's bounds or initialization expression.
    if "::" in prefix:
        entity_prefix = prefix.split("::", 1)[1]
    else:
        declaration = re.match(
            r"\s*(?:integer|real|complex|logical|character|double\s+precision|dimension)"
            r"(?:\s*\([^)]*\)|\s*\*\s*\d+)?\s+", prefix, re.I,
        )
        entity_prefix = prefix[declaration.end():] if declaration else None
    if entity_prefix is not None:
        depth = 0
        entity_start = 0
        for index, char in enumerate(entity_prefix):
            depth += (char == "(") - (char == ")")
            if char == "," and depth == 0:
                entity_start = index + 1
        if not entity_prefix[entity_start:].strip():
            return False
    # An array assignment or statement-function definition is not an
    # invocation of the mathematical intrinsic.
    depth = 0
    for index, char in enumerate(suffix):
        depth += (char == "(") - (char == ")")
        if depth == 0:
            remainder = suffix[index + 1:].lstrip()
            return not remainder.startswith("=") or remainder.startswith("==")
    return True


def fortran_statements(code, fixed_form=False):
    """Join continuations and retain the original offsets of statement text."""
    statement = ""
    offsets = []
    first_line = 1
    continued = False
    offset = 0
    for line_number, physical_line in enumerate(code.splitlines(keepends=True), 1):
        line = physical_line.rstrip("\r\n")
        indices = list(range(offset, offset + len(line)))
        offset += len(physical_line)
        if not line.strip():
            continue
        if fixed_form:
            continuation = len(line) > 5 and line[5] not in " 0"
            if statement and not continuation:
                yield statement, first_line, offsets
                statement, offsets = "", []
            line, indices = line[6:72], indices[6:72]
            if not statement:
                first_line = line_number
            statement += line.rstrip()
            offsets.extend(indices[:len(line.rstrip())])
            continue
        if not statement:
            first_line = line_number
        if continued:
            indentation = len(line) - len(line.lstrip())
            line = line.lstrip()
            indices = indices[indentation:]
            if line.startswith("&"):
                line = line[1:]
                indices = indices[1:]
        continued = line.rstrip().endswith("&")
        if continued:
            line = line.rstrip()[:-1]
            indices = indices[:len(line)]
        statement += line
        offsets.extend(indices)
        if continued:
            continue
        yield statement, first_line, offsets
        statement, offsets = "", []
    if statement:
        yield statement, first_line, offsets


def fortran_identifiers(code, fixed_form=False):
    for statement, first_line, _ in fortran_statements(code, fixed_form):
        for match in IDENTIFIER.finditer(statement):
            yield match.group(), first_line, statement, match.start(), match.end()


def shell_tokens(source):
    """Tokenize shell commands, keeping quotations and here-documents as data."""
    position = 0
    line = 1
    documents = []
    document_operator = None
    operators = (";;&", "<<-", "<<<", "&&", "||", ";;", ";&", "<<", ">>",
                 "((", "))", "<&", ">&", "&>", ">|")
    while position < len(source):
        char = source[position]
        if char in " \t\r":
            position += 1
            continue
        if source.startswith("\\\n", position):
            line += 1
            position += 2
            continue
        if char == "#":
            end = source.find("\n", position)
            position = len(source) if end < 0 else end
            continue
        if char == "\n":
            yield "operator", "\n", "\n", line
            position += 1
            line += 1
            for delimiter, strip_tabs in documents:
                while position < len(source):
                    end = source.find("\n", position)
                    end = len(source) if end < 0 else end
                    content = source[position:end]
                    position = min(end + 1, len(source))
                    line += 1
                    if (content.lstrip("\t") if strip_tabs else content) == delimiter:
                        break
            documents.clear()
            continue
        if char in ";&|(){}<>":
            operator = next((value for value in operators if source.startswith(value, position)), char)
            yield "operator", operator, operator, line
            position += len(operator)
            if operator in {"<<", "<<-"}:
                document_operator = operator
            continue
        start, first_line = position, line
        value = ""
        quote = None
        while position < len(source):
            char = source[position]
            if quote:
                if char == quote:
                    quote = None
                elif char == "\\" and quote == '"' and position + 1 < len(source):
                    position += 1
                    char = source[position]
                    if char != "\n":
                        value += char
                else:
                    value += char
                line += char == "\n"
                position += 1
            elif char in "'\"":
                quote = char
                position += 1
            elif char == "\\" and position + 1 < len(source):
                position += 1
                char = source[position]
                line += char == "\n"
                if char != "\n":
                    value += char
                position += 1
            elif char in " \t\r\n;&|(){}<>":
                break
            else:
                value += char
                position += 1
        raw = source[start:position]
        if document_operator:
            documents.append((value, document_operator == "<<-"))
            document_operator = None
            yield "document", value, raw, first_line
        else:
            yield "word", value, raw, first_line


def shell_identifiers(source):
    """Read bindings in shell command context; ordinary arguments are data."""
    tokens = list(shell_tokens(source))
    name = re.compile(r"[A-Za-z_][A-Za-z_0-9]*\Z")
    assignment = re.compile(r"([A-Za-z_][A-Za-z_0-9]*)(?:\[[^\]]*\])?\+?=")
    command_start = True
    declaration = False
    named_binding = False
    array_depth = 0
    arithmetic_depth = 0
    previous_assignment = False
    for index, (kind, value, raw, line) in enumerate(tokens):
        if kind == "document":
            continue
        if kind == "operator":
            if value == "((":
                arithmetic_depth += 1
            elif arithmetic_depth and value == "))":
                arithmetic_depth -= 1
            elif value == "(" and (previous_assignment or array_depth):
                array_depth += 1
            elif array_depth and value == ")":
                array_depth -= 1
            elif value in {";", "&&", "||", "|", "&", "\n", "(", ")", "{", "}", ";;", ";&", ";;&"}:
                command_start, declaration, named_binding = True, False, False
            previous_assignment = False
            continue
        if array_depth:
            continue
        if arithmetic_depth:
            for identifier in re.findall(r"[A-Za-z_][A-Za-z_0-9]*", value):
                yield identifier, line
            continue
        if named_binding:
            if name.fullmatch(value):
                yield value, line
            named_binding, command_start = False, False
            continue
        matched = assignment.match(value if declaration else raw)
        if declaration:
            if matched or name.fullmatch(value):
                yield matched.group(1) if matched else value, line
            previous_assignment = bool(matched and value.endswith("="))
            continue
        if not command_start:
            continue
        if matched:
            yield matched.group(1), line
            previous_assignment = value.endswith("=")
            continue
        previous_assignment = False
        if value in {"if", "then", "elif", "else", "while", "until", "do", "!", "time"}:
            continue
        if value in {"local", "declare", "typeset", "export", "readonly"}:
            declaration = True
            continue
        if value in {"for", "select", "function"}:
            named_binding = True
            continue
        following = [token[1] for token in tokens[index + 1:index + 3]]
        if name.fullmatch(raw) and following == ["(", ")"]:
            yield value, line
        command_start = False


def python_identifiers(source, path):
    tree = ast.parse(source, filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            yield node.name, node.lineno
        elif isinstance(node, ast.arg):
            yield node.arg, node.lineno
        elif isinstance(node, ast.Name):
            yield node.id, node.lineno
        elif isinstance(node, ast.Attribute) and isinstance(node.ctx, ast.Store):
            yield node.attr, node.lineno
        elif isinstance(node, ast.alias):
            yield node.asname or node.name.split(".")[0], node.lineno
        elif isinstance(node, ast.ExceptHandler) and node.name:
            yield node.name, node.lineno
        elif isinstance(node, (ast.MatchAs, ast.MatchStar)) and node.name:
            yield node.name, node.lineno
        elif isinstance(node, ast.MatchMapping) and node.rest:
            yield node.rest, node.lineno
        elif isinstance(node, (ast.Global, ast.Nonlocal)):
            for name in node.names:
                yield name, node.lineno


def check(root, banned):
    failures = set()
    paths = source_paths(root)
    vocabulary = re.compile(r"\b(?:" + "|".join(map(re.escape, sorted(banned))) + r")\b", re.I)
    for path in paths:
        source = path.read_text()
        relative_path = path.relative_to(root)
        if path.suffix.lower() in FORTRAN_SUFFIXES:
            fixed_form = path.suffix.lower() in {".f", ".for", ".ftn"}
            code = fortran_code(source, fixed_form)
            mathematical_offsets = set()
            for statement, line, offsets in fortran_statements(code, fixed_form):
                for match in IDENTIFIER.finditer(statement):
                    identifier = match.group()
                    if is_mathematical_intrinsic(statement, match.start(), match.end()):
                        mathematical_offsets.add(offsets[match.start()])
                        continue
                    if banned.intersection(identifier_components(identifier)):
                        failures.add((str(relative_path), line, "identifier", identifier))
            # Preserve the existing stricter prose rule for production Fortran.
            if relative_path.parts[0] in {"src", "application"}:
                for match in vocabulary.finditer(source):
                    if (match.group().lower() == "ceiling"
                            and match.start() in mathematical_offsets):
                        continue
                    line = source.count("\n", 0, match.start()) + 1
                    failures.add((str(relative_path), line, "vocabulary", match.group()))
        else:
            identifiers = python_identifiers(source, path) if path.suffix.lower() == ".py" else shell_identifiers(source)
            for identifier, line in identifiers:
                if banned.intersection(identifier_components(identifier)):
                    failures.add((str(relative_path), line, "identifier", identifier))
    for path, line, category, identifier in sorted(failures):
        print(f"{path}:{line}: colloquial {category}: {identifier}")
    if failures:
        print(f" FAIL : mathematical naming rejected {len(failures)} occurrences")
        return 1
    print(f" PASS : mathematical identifiers in all {len(paths)} source files and production vocabulary")
    return 0


if __name__ == "__main__":
    sys.exit(check(Path(sys.argv[1]).resolve(), set(os.environ["NAMING_BANNED"].split("|"))))
