#!/usr/bin/env python3

from pathlib import Path
import unittest

import pipeline


def task(*, naming_required=False):
    return {
        "protocol_version": pipeline.PROTOCOL_VERSION,
        "task_id": "protocol-unit-test",
        "objective": "Exercise report validation.",
        "base_ref": "HEAD",
        "policy_files": ["AGENTS.md"],
        "sealed_facts": ["No repository access is required."],
        "scope": {
            "naming_audit_required": naming_required,
            "allowed_paths": ["src/allowed/**"],
            "forbidden_paths": ["src/allowed/forbidden.f90"],
        },
        "requirements": [{"id": "R1", "statement": "One obligation."}],
        "required_checks": [
            {"id": "T1", "command": "true", "expected": "Exit status 0."}
        ],
        "stop_conditions": ["A required check cannot run."],
        "deliverables": ["A structured report."],
        "execution": {
            "mode": "write",
            "model": "sonnet",
            "effort": "low",
            "timeout_seconds": 60,
            "max_budget_usd": None,
            "allow_network": False,
        },
    }


def manifest():
    return {
        "task_id": "protocol-unit-test",
        "base_head": "a" * 40,
    }


def report():
    return {
        "protocol_version": pipeline.PROTOCOL_VERSION,
        "task_id": "protocol-unit-test",
        "status": "candidate",
        "summary": "Candidate report.",
        "base_head": "a" * 40,
        "changed_files": [
            {"path": "src/allowed/example.f90", "purpose": "Test candidate."}
        ],
        "requirements": [
            {"id": "R1", "status": "satisfied", "evidence": ["Observed."]}
        ],
        "checks": [
            {
                "id": "T1",
                "command": "true",
                "status": "pass",
                "exit_code": 0,
                "evidence": ["Exit status 0."],
            }
        ],
        "naming_audit": {
            "status": "canonical",
            "surface": ["src/allowed/example.f90"],
            "findings": [],
        },
        "assumptions": [],
        "risks": [],
        "questions": [],
        "proposed_next_steps": [],
    }


class ScopeTests(unittest.TestCase):
    def test_directory_and_glob_patterns(self):
        self.assertTrue(pipeline.pattern_matches("src/a.f90", "src"))
        self.assertTrue(pipeline.pattern_matches("src/a.f90", "src/**"))
        self.assertFalse(pipeline.pattern_matches("test/a.f90", "src/**"))

    def test_forbidden_path_overrides_allowed_path(self):
        violations = pipeline.paths_in_scope(
            ["src/allowed/forbidden.f90"], task()["scope"]
        )
        self.assertEqual(
            ["forbidden path changed: src/allowed/forbidden.f90"], violations
        )

    def test_dirty_overlap_ignores_unrelated_changes(self):
        overlap = pipeline.overlapping_allowed_paths(
            ["AGENTS.md", "src/allowed/example.f90"], task()["scope"]
        )
        self.assertEqual(["src/allowed/example.f90"], overlap)


class ReportTests(unittest.TestCase):
    def test_complete_candidate_is_valid(self):
        violations = pipeline.validate_report(
            report(), task(naming_required=True), manifest(), ["src/allowed/example.f90"]
        )
        self.assertEqual([], violations)

    def test_candidate_cannot_hide_unpassed_check_or_question(self):
        value = report()
        value["checks"][0]["status"] = "not_run"
        value["checks"][0]["exit_code"] = None
        value["questions"] = ["May I substitute another check?"]
        violations = pipeline.validate_report(
            value, task(), manifest(), ["src/allowed/example.f90"]
        )
        self.assertIn("candidate has unpassed required check: T1", violations)
        self.assertIn("candidate contains unresolved decision questions", violations)

    def test_blocked_report_may_preserve_unrun_check_and_question(self):
        value = report()
        value["status"] = "blocked"
        value["checks"][0]["status"] = "not_run"
        value["checks"][0]["exit_code"] = None
        value["questions"] = ["The exact command is unavailable."]
        violations = pipeline.validate_report(
            value, task(), manifest(), ["src/allowed/example.f90"]
        )
        self.assertEqual([], violations)

    def test_required_naming_audit_rejects_not_applicable(self):
        value = report()
        value["naming_audit"] = {
            "status": "not_applicable",
            "surface": [],
            "findings": [],
        }
        violations = pipeline.validate_report(
            value, task(naming_required=True), manifest(), ["src/allowed/example.f90"]
        )
        self.assertIn(
            "changed candidate reports naming audit as not_applicable", violations
        )


class RequestTests(unittest.TestCase):
    def test_revision_exposes_worker_copy_not_canonical_review(self):
        value = {
            "task_id": "protocol-unit-test",
            "base_head": "a" * 40,
            "worktree": "/worker",
            "policy": [],
            "task_path": "/run/worker-input/task.json",
            "last_review": "/run/reviews/review-001.json",
            "worker_review": "/run/worker-input/reviews/review-001.json",
        }
        request = pipeline.worker_request(
            Path("/run"), value, {"candidate_digest": "digest"}
        )
        self.assertIn(value["worker_review"], request)
        self.assertNotIn(value["last_review"], request)


if __name__ == "__main__":
    unittest.main()
