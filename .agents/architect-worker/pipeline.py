#!/usr/bin/env python3
"""Deterministic transport between a Codex architect and a Claude worker."""

from __future__ import annotations

import argparse
import datetime as dt
import fnmatch
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import uuid


PROTOCOL_VERSION = "architect-worker-v1"
TASK_ID_RE = re.compile(r"^[a-z0-9][a-z0-9-]{2,63}$")
HERE = Path(__file__).resolve().parent
WORKER_PROMPT = HERE / "WORKER.md"
REPORT_SCHEMA = HERE / "worker-report.schema.json"
TASK_TEMPLATE = HERE / "task.template.json"
REVIEW_TEMPLATE = HERE / "review.template.json"


class ProtocolError(RuntimeError):
    pass


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def run(
    command: list[str],
    *,
    cwd: Path,
    check: bool = True,
    text: bool = True,
    stdout=None,
    stderr=None,
    timeout: int | None = None,
):
    result = subprocess.run(
        command,
        cwd=cwd,
        check=False,
        text=text,
        stdout=stdout if stdout is not None else subprocess.PIPE,
        stderr=stderr if stderr is not None else subprocess.PIPE,
        timeout=timeout,
    )
    if check and result.returncode != 0:
        detail = result.stderr.strip() if text and result.stderr else ""
        raise ProtocolError(f"command failed ({result.returncode}): {' '.join(command)}\n{detail}")
    return result


def git(repo: Path, *args: str, check: bool = True) -> str:
    return run(["git", *args], cwd=repo, check=check).stdout.strip()


def repository_root() -> Path:
    return Path(git(Path.cwd(), "rev-parse", "--show-toplevel")).resolve()


def common_git_dir(repo: Path) -> Path:
    value = git(repo, "rev-parse", "--path-format=absolute", "--git-common-dir")
    return Path(value).resolve()


def state_paths(repo: Path, task_id: str | None = None) -> tuple[Path, Path, Path | None]:
    state_root = common_git_dir(repo) / "architect-worker"
    workers_root = repo.parent / f".{repo.name}-claude-workers"
    run_dir = state_root / "runs" / task_id if task_id else None
    return state_root, workers_root, run_dir


def load_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ProtocolError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ProtocolError(f"expected a JSON object in {path}")
    return value


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def write_new_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def append_event(run_dir: Path, state: str, **detail) -> None:
    event = {"at": now(), "state": state, **detail}
    with (run_dir / "events.jsonl").open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(event, sort_keys=True) + "\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_relative_path(repo: Path, value: str) -> Path:
    if not value or Path(value).is_absolute():
        raise ProtocolError(f"path must be nonempty and repository-relative: {value!r}")
    candidate = (repo / value).resolve()
    try:
        candidate.relative_to(repo)
    except ValueError as exc:
        raise ProtocolError(f"path escapes repository: {value}") from exc
    return candidate


def validate_task(task: dict) -> None:
    required = {
        "protocol_version",
        "task_id",
        "objective",
        "base_ref",
        "policy_files",
        "sealed_facts",
        "scope",
        "requirements",
        "required_checks",
        "stop_conditions",
        "deliverables",
        "execution",
    }
    missing = sorted(required - task.keys())
    if missing:
        raise ProtocolError(f"task packet missing keys: {', '.join(missing)}")
    if task["protocol_version"] != PROTOCOL_VERSION:
        raise ProtocolError("task packet protocol_version does not match launcher")
    if not TASK_ID_RE.fullmatch(task["task_id"]):
        raise ProtocolError("task_id must match ^[a-z0-9][a-z0-9-]{2,63}$")
    if not isinstance(task["objective"], str) or not task["objective"].strip():
        raise ProtocolError("objective must be a nonempty string")
    if not isinstance(task["scope"], dict):
        raise ProtocolError("scope must be an object")
    for key in ("allowed_paths", "forbidden_paths"):
        values = task["scope"].get(key)
        if not isinstance(values, list) or not all(isinstance(item, str) for item in values):
            raise ProtocolError(f"scope.{key} must be an array of strings")
    naming_required = task["scope"].get("naming_audit_required", False)
    if not isinstance(naming_required, bool):
        raise ProtocolError("scope.naming_audit_required must be true or false")
    execution = task["execution"]
    if not isinstance(execution, dict) or execution.get("mode") not in {"read-only", "write"}:
        raise ProtocolError("execution.mode must be read-only or write")
    if execution.get("effort", "xhigh") not in {"low", "medium", "high", "xhigh", "max"}:
        raise ProtocolError("execution.effort is invalid")
    timeout = execution.get("timeout_seconds", 3600)
    if not isinstance(timeout, int) or timeout < 30:
        raise ProtocolError("execution.timeout_seconds must be an integer >= 30")
    budget = execution.get("max_budget_usd")
    if budget is not None and (not isinstance(budget, (int, float)) or budget <= 0):
        raise ProtocolError("execution.max_budget_usd must be null or positive")
    for collection in ("policy_files", "sealed_facts", "stop_conditions", "deliverables"):
        if not isinstance(task[collection], list) or not all(
            isinstance(item, str) and item.strip() for item in task[collection]
        ):
            raise ProtocolError(f"{collection} must be an array of nonempty strings")
    validate_id_records(task["requirements"], "requirements", ("statement",))
    validate_id_records(task["required_checks"], "required_checks", ("command", "expected"))


def validate_id_records(records, label: str, required_fields: tuple[str, ...]) -> None:
    if not isinstance(records, list):
        raise ProtocolError(f"{label} must be an array")
    seen: set[str] = set()
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get("id"), str):
            raise ProtocolError(f"each {label} entry needs a string id")
        if record["id"] in seen:
            raise ProtocolError(f"duplicate {label} id: {record['id']}")
        seen.add(record["id"])
        for field in required_fields:
            if not isinstance(record.get(field), str) or not record[field].strip():
                raise ProtocolError(f"{label} {record['id']} needs nonempty {field}")


def pattern_matches(path: str, pattern: str) -> bool:
    normalized = pattern.rstrip("/")
    if not any(mark in normalized for mark in "*?["):
        return path == normalized or path.startswith(normalized + "/")
    return fnmatch.fnmatchcase(path, normalized)


def paths_in_scope(paths: list[str], scope: dict) -> list[str]:
    allowed = scope["allowed_paths"]
    forbidden = scope["forbidden_paths"]
    violations = []
    for path in paths:
        if any(pattern_matches(path, item) for item in forbidden):
            violations.append(f"forbidden path changed: {path}")
        elif not any(pattern_matches(path, item) for item in allowed):
            violations.append(f"out-of-scope path changed: {path}")
    return violations


def overlapping_allowed_paths(paths: list[str], scope: dict) -> list[str]:
    return sorted(
        path
        for path in paths
        if any(pattern_matches(path, pattern) for pattern in scope["allowed_paths"])
    )


def changed_paths(repo: Path, base_head: str) -> list[str]:
    tracked = run(
        ["git", "diff", "--name-only", "-z", base_head, "--"],
        cwd=repo,
        text=False,
    ).stdout.split(b"\0")
    untracked = run(
        ["git", "ls-files", "--others", "--exclude-standard", "-z"],
        cwd=repo,
        text=False,
    ).stdout.split(b"\0")
    values = {
        os.fsdecode(item)
        for item in tracked + untracked
        if item
    }
    return sorted(values)


def current_checkout_changes(repo: Path) -> list[str]:
    head = git(repo, "rev-parse", "HEAD")
    return changed_paths(repo, head)


def candidate_index(worktree: Path, paths: list[str]) -> tuple[list[dict], str]:
    records = []
    digest = hashlib.sha256()
    for relative in paths:
        path = worktree / relative
        if path.is_symlink():
            record = {"path": relative, "kind": "symlink", "target": os.readlink(path)}
        elif path.is_file():
            record = {
                "path": relative,
                "kind": "file",
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        elif path.exists():
            record = {"path": relative, "kind": "other"}
        else:
            record = {"path": relative, "kind": "deleted"}
        encoded = json.dumps(record, sort_keys=True).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
        records.append(record)
    return records, digest.hexdigest()


def snapshot_candidate(turn_dir: Path, worktree: Path, base_head: str, paths: list[str]) -> str:
    patch = run(
        ["git", "diff", "--binary", "--no-ext-diff", base_head, "--"],
        cwd=worktree,
        text=False,
    ).stdout
    (turn_dir / "diff.patch").write_bytes(patch)
    records, digest = candidate_index(worktree, paths)
    write_new_json(turn_dir / "candidate.json", {"digest": digest, "files": records})
    with tarfile.open(turn_dir / "candidate.tar.gz", "w:gz") as archive:
        for relative in paths:
            path = worktree / relative
            if path.exists() or path.is_symlink():
                archive.add(path, arcname=relative, recursive=False)
    return digest


def load_manifest(repo: Path, task_id: str) -> tuple[Path, dict]:
    if not TASK_ID_RE.fullmatch(task_id):
        raise ProtocolError("invalid task_id")
    _, _, run_dir = state_paths(repo, task_id)
    assert run_dir is not None
    manifest_path = run_dir / "manifest.json"
    if not manifest_path.is_file():
        raise ProtocolError(f"unknown task: {task_id}")
    return run_dir, load_json(manifest_path)


def save_manifest(run_dir: Path, manifest: dict) -> None:
    manifest["updated_at"] = now()
    write_json(run_dir / "manifest.json", manifest)


def manifest_task_path(run_dir: Path, manifest: dict) -> Path:
    return Path(manifest.get("task_path", run_dir / "task.json"))


def immutable_input_violations(run_dir: Path, manifest: dict) -> list[str]:
    violations = []
    task_path = manifest_task_path(run_dir, manifest)
    if not task_path.is_file() or sha256_file(task_path) != manifest["task_sha256"]:
        violations.append("immutable task packet changed during worker execution")
    for policy in manifest["policy"]:
        snapshot = Path(policy["snapshot"])
        if not snapshot.is_file() or sha256_file(snapshot) != policy["sha256"]:
            violations.append(f"immutable policy snapshot changed: {policy['path']}")
    review_path = manifest.get("last_review")
    review_sha256 = manifest.get("last_review_sha256")
    if review_path and review_sha256:
        review = Path(review_path)
        if not review.is_file() or sha256_file(review) != review_sha256:
            violations.append("immutable architect review changed during worker execution")
    worker_review_path = manifest.get("worker_review")
    worker_review_sha256 = manifest.get("worker_review_sha256")
    if worker_review_path and worker_review_sha256:
        worker_review = Path(worker_review_path)
        if not worker_review.is_file() or sha256_file(worker_review) != worker_review_sha256:
            violations.append("immutable worker-visible review changed during worker execution")
    return violations


def command_doctor(_: argparse.Namespace) -> None:
    repo = repository_root()
    schema = load_json(REPORT_SCHEMA)
    if schema.get("title") != "Architect-worker Claude report":
        raise ProtocolError("worker report schema is not the expected file")
    validate_task(load_json(TASK_TEMPLATE))
    review_template = load_json(REVIEW_TEMPLATE)
    if review_template.get("protocol_version") != PROTOCOL_VERSION:
        raise ProtocolError("review template protocol_version does not match launcher")
    claude_path = shutil.which("claude")
    if not claude_path:
        raise ProtocolError("claude executable is not available")
    version = run([claude_path, "--version"], cwd=repo).stdout.strip()
    result = {
        "protocol_version": PROTOCOL_VERSION,
        "repository": str(repo),
        "head": git(repo, "rev-parse", "HEAD"),
        "claude": claude_path,
        "claude_version": version,
        "worker_contract": str(WORKER_PROMPT),
        "report_schema": str(REPORT_SCHEMA),
        "task_template": str(TASK_TEMPLATE),
        "review_template": str(REVIEW_TEMPLATE),
        "status": "ready",
    }
    print(json.dumps(result, indent=2, sort_keys=True))


def command_prepare(args: argparse.Namespace) -> None:
    repo = repository_root()
    task_path = Path(args.task).resolve()
    task = load_json(task_path)
    validate_task(task)
    task_id = task["task_id"]
    state_root, workers_root, run_dir = state_paths(repo, task_id)
    assert run_dir is not None
    worktree = workers_root / task_id
    if run_dir.exists() or worktree.exists():
        raise ProtocolError(f"task already exists: {task_id}")

    overlap = overlapping_allowed_paths(current_checkout_changes(repo), task["scope"])
    if overlap:
        raise ProtocolError(
            "main checkout has changes overlapping delegated scope:\n" + "\n".join(overlap)
        )

    base_head = git(repo, "rev-parse", f"{task['base_ref']}^{{commit}}")
    policy_sources = []
    for relative in task["policy_files"]:
        source = validate_relative_path(repo, relative)
        if not source.is_file():
            raise ProtocolError(f"policy file does not exist: {relative}")
        policy_sources.append((relative, source))

    workers_root.mkdir(parents=True, exist_ok=True)
    git(repo, "worktree", "add", "--detach", str(worktree), base_head)

    run_dir.mkdir(parents=True)
    worker_input_dir = run_dir / "worker-input"
    worker_input_dir.mkdir()
    task_snapshot = worker_input_dir / "task.json"
    shutil.copyfile(task_path, task_snapshot)
    policy_records = []
    for relative, source in policy_sources:
        target = worker_input_dir / "policy" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        policy_records.append(
            {"path": relative, "snapshot": str(target), "sha256": sha256_file(target)}
        )

    manifest = {
        "protocol_version": PROTOCOL_VERSION,
        "task_id": task_id,
        "state": "DELEGATED",
        "repository": str(repo),
        "worktree": str(worktree),
        "base_head": base_head,
        "session_id": str(uuid.uuid4()),
        "task_path": str(task_snapshot),
        "task_sha256": sha256_file(task_snapshot),
        "worker_input_dir": str(worker_input_dir),
        "policy": policy_records,
        "turns": 0,
        "reviews": 0,
        "candidate_digest": None,
        "last_review": None,
        "archived": False,
        "created_at": now(),
        "updated_at": now(),
    }
    manifest["session_history"] = [manifest["session_id"]]
    write_new_json(run_dir / "manifest.json", manifest)
    append_event(run_dir, "DELEGATED", base_head=base_head, worktree=str(worktree))
    print(json.dumps({"task_id": task_id, "state": "DELEGATED", "run_dir": str(run_dir), "worktree": str(worktree)}, indent=2))


def worker_request(run_dir: Path, manifest: dict, revision: dict | None) -> str:
    policies = "\n".join(f"- {item['snapshot']}" for item in manifest["policy"])
    request = f"""Architect-worker delegation

Protocol: {PROTOCOL_VERSION}
Task ID: {manifest['task_id']}
Base HEAD: {manifest['base_head']}
Detached worktree: {manifest['worktree']}

Read these immutable policy snapshots first:
{policies}

Then read the immutable task packet:
- {manifest_task_path(run_dir, manifest)}

Execute only that packet in the current worktree. Return the structured worker
report. The launcher, not you, will decide whether the actual diff matches the
report and scope.
"""
    if revision is not None:
        request += f"""

This is a revision turn. The original task is unchanged. Read the architect's
latest review and address only its bounded findings:
- {manifest.get('worker_review', manifest['last_review'])}

Candidate digest reviewed by the architect: {revision['candidate_digest']}
"""
    return request


def claude_command(
    task: dict,
    manifest: dict,
    run_dir: Path,
    request: str,
    *,
    resume: bool,
) -> list[str]:
    execution = task["execution"]
    mode = execution["mode"]
    tools = "Read,Glob,Grep" if mode == "read-only" else "Read,Glob,Grep,Edit,Write,Bash"
    permission_mode = "dontAsk" if mode == "read-only" else "acceptEdits"
    schema = REPORT_SCHEMA.read_text(encoding="utf-8")
    worker_contract = WORKER_PROMPT.read_text(encoding="utf-8")
    command = [
        "claude",
        "-p",
        "--model",
        execution.get("model", "opus"),
        "--effort",
        execution.get("effort", "xhigh"),
        "--safe-mode",
        "--setting-sources",
        "project",
        "--disable-slash-commands",
        "--no-chrome",
        "--permission-mode",
        permission_mode,
        "--tools",
        tools,
        "--add-dir",
        str(Path(manifest.get("worker_input_dir", run_dir)).resolve()),
        "--append-system-prompt",
        worker_contract,
        "--output-format",
        "json",
        "--json-schema",
        schema,
        "--name",
        f"worker-{manifest['task_id']}",
    ]
    if not execution.get("allow_network", False):
        command.extend(["--disallowedTools", "WebSearch,WebFetch"])
    budget = execution.get("max_budget_usd")
    if budget is not None:
        command.extend(["--max-budget-usd", str(budget)])
    if resume:
        command.extend(["--resume", manifest["session_id"]])
    else:
        command.extend(["--session-id", manifest["session_id"]])
    command.append(request)
    return command


def validate_report(report: dict, task: dict, manifest: dict, actual_paths: list[str]) -> list[str]:
    violations = []
    if report.get("protocol_version") != PROTOCOL_VERSION:
        violations.append("worker report protocol_version mismatch")
    if report.get("task_id") != manifest["task_id"]:
        violations.append("worker report task_id mismatch")
    if report.get("base_head") != manifest["base_head"]:
        violations.append("worker report base_head mismatch")

    reported_paths = sorted(
        item.get("path", "")
        for item in report.get("changed_files", [])
        if isinstance(item, dict)
    )
    if reported_paths != actual_paths:
        violations.append(
            f"reported changed paths {reported_paths} do not equal actual paths {actual_paths}"
        )

    required_ids = {item["id"] for item in task["requirements"]}
    reported_requirements = {
        item.get("id"): item
        for item in report.get("requirements", [])
        if isinstance(item, dict)
    }
    missing_requirements = sorted(required_ids - reported_requirements.keys())
    if missing_requirements:
        violations.append(f"requirements absent from report: {missing_requirements}")

    required_checks = {item["id"]: item["command"] for item in task["required_checks"]}
    reported_checks = {
        item.get("id"): item
        for item in report.get("checks", [])
        if isinstance(item, dict)
    }
    for check_id, command in required_checks.items():
        if reported_checks.get(check_id, {}).get("command") != command:
            violations.append(f"required check missing or command changed: {check_id}")

    if report.get("status") == "candidate":
        for requirement_id in required_ids:
            if reported_requirements.get(requirement_id, {}).get("status") != "satisfied":
                violations.append(f"candidate has unsatisfied requirement: {requirement_id}")
        for check_id in required_checks:
            check = reported_checks.get(check_id, {})
            if check.get("status") != "pass" or check.get("exit_code") != 0:
                violations.append(f"candidate has unpassed required check: {check_id}")
        if report.get("questions"):
            violations.append("candidate contains unresolved decision questions")

    if task["execution"]["mode"] == "read-only" and actual_paths:
        violations.append("read-only task changed files")
    naming = report.get("naming_audit", {})
    if (
        actual_paths
        and task["scope"].get("naming_audit_required", False)
        and naming.get("status") == "not_applicable"
    ):
        violations.append("changed candidate reports naming audit as not_applicable")
    violations.extend(paths_in_scope(actual_paths, task["scope"]))
    return violations


def execute_turn(
    repo: Path,
    run_dir: Path,
    manifest: dict,
    *,
    continuation: str,
) -> None:
    task = load_json(manifest_task_path(run_dir, manifest))
    worktree = Path(manifest["worktree"])
    if not worktree.is_dir():
        raise ProtocolError("worker worktree is missing")
    if continuation not in {"dispatch", "revise", "retry"}:
        raise ProtocolError(f"invalid continuation mode: {continuation}")
    expected_states = {"DELEGATED"} if continuation == "dispatch" else {"REVISION_REQUESTED"}
    if manifest["state"] not in expected_states:
        raise ProtocolError(f"cannot execute from state {manifest['state']}")

    is_revision = continuation in {"revise", "retry"}
    resume_session = continuation == "revise"
    revision = load_json(Path(manifest["last_review"])) if is_revision else None
    turn_number = manifest["turns"] + 1
    turn_dir = run_dir / f"turn-{turn_number:03d}"
    turn_dir.mkdir(parents=True)
    request = worker_request(run_dir, manifest, revision)
    (turn_dir / "request.md").write_text(request, encoding="utf-8")

    manifest["state"] = "RUNNING"
    save_manifest(run_dir, manifest)
    append_event(run_dir, "RUNNING", turn=turn_number, continuation=continuation)

    command = claude_command(task, manifest, run_dir, request, resume=resume_session)
    stdout_path = turn_dir / "envelope.json"
    stderr_path = turn_dir / "claude.stderr"
    timeout_seconds = task["execution"].get("timeout_seconds", 3600)
    return_code = None
    failure = None
    try:
        with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
            result = run(
                command,
                cwd=worktree,
                check=False,
                text=False,
                stdout=stdout,
                stderr=stderr,
                timeout=timeout_seconds,
            )
            return_code = result.returncode
    except subprocess.TimeoutExpired:
        failure = f"Claude exceeded timeout of {timeout_seconds} seconds"

    manifest["turns"] = turn_number
    if failure or return_code != 0:
        manifest["state"] = "FAILED"
        save_manifest(run_dir, manifest)
        append_event(run_dir, "FAILED", turn=turn_number, reason=failure or f"exit {return_code}")
        raise ProtocolError(failure or f"Claude exited with status {return_code}")

    envelope = load_json(stdout_path)
    report = envelope.get("structured_output")
    if not isinstance(report, dict):
        manifest["state"] = "FAILED"
        save_manifest(run_dir, manifest)
        append_event(run_dir, "FAILED", turn=turn_number, reason="missing structured_output")
        raise ProtocolError("Claude envelope has no structured_output object")
    write_new_json(turn_dir / "report.json", report)

    actual_head = git(worktree, "rev-parse", "HEAD")
    actual_paths = changed_paths(worktree, manifest["base_head"])
    digest = snapshot_candidate(turn_dir, worktree, manifest["base_head"], actual_paths)
    violations = validate_report(report, task, manifest, actual_paths)
    violations.extend(immutable_input_violations(run_dir, manifest))
    if actual_head != manifest["base_head"]:
        violations.append(f"worker moved HEAD from {manifest['base_head']} to {actual_head}")

    if violations:
        state = "PROTOCOL_VIOLATION"
        write_new_json(turn_dir / "protocol-violations.json", {"violations": violations})
    else:
        state = {
            "candidate": "CANDIDATE",
            "blocked": "BLOCKED",
            "failed": "FAILED",
        }.get(report.get("status"), "FAILED")

    manifest["state"] = state
    manifest["candidate_digest"] = digest
    manifest["last_report"] = str(turn_dir / "report.json")
    save_manifest(run_dir, manifest)
    append_event(
        run_dir,
        state,
        turn=turn_number,
        candidate_digest=digest,
        changed_paths=actual_paths,
        violations=violations,
    )
    print(
        json.dumps(
            {
                "task_id": manifest["task_id"],
                "state": state,
                "turn": turn_number,
                "candidate_digest": digest,
                "changed_paths": actual_paths,
                "report": str(turn_dir / "report.json"),
                "violations": violations,
            },
            indent=2,
            sort_keys=True,
        )
    )


def command_dispatch(_: argparse.Namespace) -> None:
    repo = repository_root()
    task_id = _.task_id
    run_dir, manifest = load_manifest(repo, task_id)
    execute_turn(repo, run_dir, manifest, continuation="dispatch")


def validate_review(review: dict, task: dict, manifest: dict) -> str:
    required = {
        "protocol_version",
        "task_id",
        "candidate_digest",
        "decision",
        "summary",
        "findings",
        "independent_checks",
        "integration",
        "next_steps",
    }
    missing = sorted(required - review.keys())
    if missing:
        raise ProtocolError(f"review missing keys: {missing}")
    if review["protocol_version"] != PROTOCOL_VERSION:
        raise ProtocolError("review protocol_version mismatch")
    if review["task_id"] != manifest["task_id"]:
        raise ProtocolError("review task_id mismatch")
    if review["candidate_digest"] != manifest["candidate_digest"]:
        raise ProtocolError("review candidate_digest is stale or incorrect")
    decision = review["decision"]
    transitions = {
        "revision_requested": ({"CANDIDATE", "BLOCKED", "FAILED", "PROTOCOL_VIOLATION", "VERIFIED", "INTEGRATED"}, "REVISION_REQUESTED"),
        "verified": ({"CANDIDATE"}, "VERIFIED"),
        "integrated": ({"VERIFIED"}, "INTEGRATED"),
        "complete": ({"VERIFIED"} if task["execution"]["mode"] == "read-only" else {"INTEGRATED"}, "COMPLETE"),
        "rejected": ({"CANDIDATE", "BLOCKED", "FAILED", "PROTOCOL_VIOLATION", "VERIFIED"}, "REJECTED"),
        "escalated": ({"CANDIDATE", "BLOCKED", "FAILED", "PROTOCOL_VIOLATION", "VERIFIED", "INTEGRATED"}, "ESCALATED"),
    }
    if decision not in transitions:
        raise ProtocolError(f"unknown review decision: {decision}")
    allowed, state = transitions[decision]
    if manifest["state"] not in allowed:
        raise ProtocolError(f"decision {decision} is invalid from state {manifest['state']}")

    check_by_id = {
        item.get("id"): item for item in review["independent_checks"] if isinstance(item, dict)
    }
    if decision in {"verified", "complete"}:
        if review["findings"]:
            raise ProtocolError(f"decision {decision} cannot contain unresolved findings")
        for required_check in task["required_checks"]:
            check = check_by_id.get(required_check["id"])
            if not check or check.get("command") != required_check["command"]:
                raise ProtocolError(f"architect review lacks required check {required_check['id']}")
            if check.get("status") != "pass" or check.get("exit_code") != 0:
                raise ProtocolError(f"architect check has not passed: {required_check['id']}")
    if decision == "integrated" and review.get("integration", {}).get("status") != "integrated":
        raise ProtocolError("integrated decision requires integration.status=integrated")
    return state


def command_review(args: argparse.Namespace) -> None:
    repo = repository_root()
    run_dir, manifest = load_manifest(repo, args.task_id)
    task = load_json(manifest_task_path(run_dir, manifest))
    review = load_json(Path(args.review).resolve())
    new_state = validate_review(review, task, manifest)
    review_number = manifest["reviews"] + 1
    target = run_dir / "reviews" / f"review-{review_number:03d}.json"
    write_new_json(target, review)
    worker_input_dir = Path(manifest.get("worker_input_dir", run_dir))
    if worker_input_dir == run_dir:
        worker_review = target
    else:
        worker_review = worker_input_dir / "reviews" / target.name
        worker_review.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(target, worker_review)
    manifest["reviews"] = review_number
    manifest["last_review"] = str(target)
    manifest["last_review_sha256"] = sha256_file(target)
    manifest["worker_review"] = str(worker_review)
    manifest["worker_review_sha256"] = sha256_file(worker_review)
    manifest["state"] = new_state
    save_manifest(run_dir, manifest)
    append_event(run_dir, new_state, review=review_number, decision=review["decision"])
    print(json.dumps({"task_id": args.task_id, "state": new_state, "review": str(target)}, indent=2))


def command_revise(args: argparse.Namespace) -> None:
    repo = repository_root()
    run_dir, manifest = load_manifest(repo, args.task_id)
    execute_turn(repo, run_dir, manifest, continuation="revise")


def command_retry(args: argparse.Namespace) -> None:
    repo = repository_root()
    run_dir, manifest = load_manifest(repo, args.task_id)
    if manifest["state"] != "REVISION_REQUESTED":
        raise ProtocolError(f"cannot retry from state {manifest['state']}")
    previous_session = manifest["session_id"]
    manifest["session_id"] = str(uuid.uuid4())
    manifest.setdefault("session_history", [previous_session]).append(manifest["session_id"])
    save_manifest(run_dir, manifest)
    append_event(
        run_dir,
        "REVISION_REQUESTED",
        event="SESSION_REPLACED",
        previous_session=previous_session,
        new_session=manifest["session_id"],
    )
    execute_turn(repo, run_dir, manifest, continuation="retry")


def command_status(args: argparse.Namespace) -> None:
    repo = repository_root()
    run_dir, manifest = load_manifest(repo, args.task_id)
    print(json.dumps({**manifest, "run_dir": str(run_dir)}, indent=2, sort_keys=True))


def command_archive(args: argparse.Namespace) -> None:
    repo = repository_root()
    run_dir, manifest = load_manifest(repo, args.task_id)
    if manifest["state"] not in {"COMPLETE", "REJECTED"}:
        raise ProtocolError("archive requires COMPLETE or REJECTED")
    if manifest.get("archived"):
        raise ProtocolError("task worktree is already archived")
    worktree = Path(manifest["worktree"]).resolve()
    _, workers_root, _ = state_paths(repo, args.task_id)
    try:
        worktree.relative_to(workers_root.resolve())
    except ValueError as exc:
        raise ProtocolError("refusing to remove worktree outside managed worker root") from exc
    if changed_paths(worktree, manifest["base_head"]):
        raise ProtocolError("worker worktree is not clean; preserve or integrate its candidate first")
    git(repo, "worktree", "remove", str(worktree))
    manifest["archived"] = True
    save_manifest(run_dir, manifest)
    append_event(run_dir, manifest["state"], archived=True, worktree=str(worktree))
    print(json.dumps({"task_id": args.task_id, "state": manifest["state"], "archived": True}, indent=2))


def command_discard(args: argparse.Namespace) -> None:
    repo = repository_root()
    run_dir, manifest = load_manifest(repo, args.task_id)
    if manifest["state"] != "REJECTED":
        raise ProtocolError("discard requires REJECTED")
    if manifest.get("archived"):
        raise ProtocolError("task worktree is already archived")
    worktree = Path(manifest["worktree"]).resolve()
    _, workers_root, _ = state_paths(repo, args.task_id)
    try:
        worktree.relative_to(workers_root.resolve())
    except ValueError as exc:
        raise ProtocolError("refusing to discard worktree outside managed worker root") from exc
    turn_dir = run_dir / f"turn-{manifest['turns']:03d}"
    candidate = load_json(turn_dir / "candidate.json")
    if candidate.get("digest") != manifest.get("candidate_digest"):
        raise ProtocolError("latest candidate snapshot does not match manifest digest")
    if not (turn_dir / "candidate.tar.gz").is_file():
        raise ProtocolError("latest candidate archive is missing")
    git(repo, "worktree", "remove", "--force", str(worktree))
    manifest["archived"] = True
    manifest["candidate_discarded"] = True
    save_manifest(run_dir, manifest)
    append_event(
        run_dir,
        "REJECTED",
        archived=True,
        candidate_discarded=True,
        worktree=str(worktree),
    )
    print(
        json.dumps(
            {
                "task_id": args.task_id,
                "state": "REJECTED",
                "archived": True,
                "candidate_discarded": True,
            },
            indent=2,
        )
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    doctor = commands.add_parser("doctor", help="verify local transport dependencies")
    doctor.set_defaults(function=command_doctor)
    prepare = commands.add_parser("prepare", help="freeze a task packet and create its worktree")
    prepare.add_argument("task")
    prepare.set_defaults(function=command_prepare)
    dispatch = commands.add_parser("dispatch", help="start the Claude worker")
    dispatch.add_argument("task_id")
    dispatch.set_defaults(function=command_dispatch)
    status = commands.add_parser("status", help="show the current protocol state")
    status.add_argument("task_id")
    status.set_defaults(function=command_status)
    review = commands.add_parser("review", help="record an architect decision")
    review.add_argument("task_id")
    review.add_argument("review")
    review.set_defaults(function=command_review)
    revise = commands.add_parser("revise", help="resume Claude with the latest revision request")
    revise.add_argument("task_id")
    revise.set_defaults(function=command_revise)
    retry = commands.add_parser("retry", help="start a fresh Claude session after transport failure")
    retry.add_argument("task_id")
    retry.set_defaults(function=command_retry)
    archive = commands.add_parser("archive", help="remove a clean terminal worker worktree")
    archive.add_argument("task_id")
    archive.set_defaults(function=command_archive)
    discard = commands.add_parser("discard", help="remove a snapshotted rejected candidate worktree")
    discard.add_argument("task_id")
    discard.set_defaults(function=command_discard)
    return result


def main() -> int:
    try:
        args = parser().parse_args()
        args.function(args)
        return 0
    except (ProtocolError, OSError, subprocess.TimeoutExpired) as exc:
        print(f"protocol error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
