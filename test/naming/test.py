#!/usr/bin/env python3
"""Regression tests for mathematical vocabulary and staged-tree enforcement."""

import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCANNER = ROOT / "tools/check_vocabulary.py"
BANNED = re.search(r"^banned='([^']+)'", (ROOT / "check_naming.sh").read_text(), re.M).group(1)


class NamingChecks(unittest.TestCase):
    def check_source(self, path, source, accepted):
        with tempfile.TemporaryDirectory(prefix="ufvm-naming-") as directory:
            root = Path(directory)
            target = root / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(source)
            completed = subprocess.run(
                ["python3", str(SCANNER), str(root)], capture_output=True, text=True,
                env=dict(os.environ, NAMING_BANNED=BANNED), check=False,
            )
            self.assertEqual(completed.returncode == 0, accepted, completed.stdout + completed.stderr)
            if not accepted:
                self.assertIn("colloquial", completed.stdout)

    def test_fortran_identifiers(self):
        for identifier in ["kept", "kept_inner", "inner_kept", "num_kept_columns", "KEPT_INNER", "KePt_inner"]:
            with self.subTest(identifier=identifier):
                self.check_source("test/example.f90", f"integer :: {identifier}\n", False)
        self.check_source("src/util_example.f90", "real :: initial_residual_norm\ninteger :: retained_columns(2)\n", True)
        self.check_source("test/example.f90", "type :: held_context\nend type\n", False)
        self.check_source("test/example.f90", "subroutine check_budget(n)\nend subroutine\n", False)

    def test_continuations_and_literals(self):
        for source in ["integer :: ke&\n&pt_inner\n", "integer :: ke&\n! comment\n\n&pt_inner\n"]:
            self.check_source("test/example.f90", source, False)
        self.check_source("test/example.f", "      integer ke\n     &pt_inner\n", False)
        self.check_source("test/example.f90", "! kept_inner is invalid input\nprint *, 'it''s kept_inner!'\n", True)
        self.check_source("test/example.f90", 'print *, "text!"; integer :: kept_inner\n', False)
        self.check_source("src/util_example.f90", "! a kept value\ninteger :: retained\n", False)
        self.check_source("application/example.f90", "print *, 'the answer'\n", False)

    def test_mathematical_intrinsic(self):
        for source in ["n = ceiling(1.2)\n", "integer, parameter :: n = ceiling(1.2)\n", "intrinsic :: ceiling\n"]:
            self.check_source("src/util_example.f90", source, True)
        self.check_source("src/util_example.f90", "integer :: ceiling\n", False)
        self.check_source("test/example.f90", "integer :: ceiling_gb\n", False)

    def test_shell_declarations(self):
        for source in ["kept_inner=1\n", "f() { local kept_inner; }\n", "local x kept_inner=1\n",
                       "declare -i kept_inner=1\n", "true; kept_inner=1\n", "if true; then kept_inner=1; fi\n",
                       "check_budget() { :; }\n"]:
            with self.subTest(source=source):
                self.check_source("test/example.sh", source, False)
        self.check_source("test/example.sh", "cat <<'DATA'\nkept_inner=1\nDATA\n", True)
        self.check_source("test/example.sh", "text='\nkept_inner=1\n'\nprintf '%s' \"$text\"\n", True)

    def test_python_declarations(self):
        for source in ["kept_inner = 1\n", "def check_budget(x): pass\n", "class HeldContext: pass\n",
                       "try: pass\nexcept ValueError as kept_inner: pass\n",
                       "class A:\n def f(self): self.kept_inner = 1\n"]:
            with self.subTest(source=source):
                self.check_source("test/example.py", source, False)
        self.check_source("test/example.PY", "kept_inner = 1\n", False)
        self.check_source("test/example.py", "import ast\nfor node in ast.walk(ast.parse('x=1')): pass\n", True)
        self.check_source("test/example.py", "source = 'kept_inner = 1'\n", True)

    def test_staged_tree_is_authoritative(self):
        with tempfile.TemporaryDirectory(prefix="ufvm-naming-index-") as directory:
            root = Path(directory)
            for relative in ["check_naming.sh", "tools/check_vocabulary.py", "githooks/pre-commit"]:
                target = root / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(ROOT / relative, target)
            (root / "src").mkdir()
            source = root / "src/util_example.f90"
            source.write_text("module util_example\ninteger :: retained_columns\nend module\n")
            # OBJECTS entries follow the library's one-module-per-line format.
            (root / "src/OBJECTS").write_text("objects= \\\nutil_example.o\n")
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            subprocess.run(["git", "-C", str(root), "add", "."], check=True)
            source.write_text("module util_example\ninteger :: kept_columns\nend module\n")
            accepted = subprocess.run(["bash", "githooks/pre-commit"], cwd=root, capture_output=True, text=True, check=False)
            self.assertEqual(accepted.returncode, 0, accepted.stdout + accepted.stderr)
            subprocess.run(["git", "-C", str(root), "add", "src/util_example.f90"], check=True)
            source.write_text("module util_example\ninteger :: retained_columns\nend module\n")
            refused = subprocess.run(["bash", "githooks/pre-commit"], cwd=root, capture_output=True, text=True, check=False)
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("kept_columns", refused.stdout)


if __name__ == "__main__":
    unittest.main()
