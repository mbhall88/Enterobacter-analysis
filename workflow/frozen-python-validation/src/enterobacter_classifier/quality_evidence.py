"""Attach pinned local CheckM2 evidence, without interpreting a quality gate.

This development adapter checks artifact consistency, not execution authenticity
or biological purity. Native estimates remain separate from structural metrics.
It supports the recorded CheckM2 1.1.0/database3 producer contract only.
"""
from __future__ import annotations

import copy
import gzip
import hashlib
import json
import math
import re
from pathlib import Path

from .evidence import explicit_percent, read_tsv, sha256
from .quality import SCHEMA_VERSION as STRUCTURAL_SCHEMA

SCHEMA = "1.0.0-dev-measured-quality-evidence"
DATABASE_SHA256 = "1b86ef3eac0813c1853f53182c17657045e3763d66f384ec95747261a63ae46f"
NATIVE_STRUCTURE = {"Genome_Size": "total_bp", "Total_Contigs": "record_count",
                    "Max_Contig_Length": "max_record_bp", "Contig_N50": "record_n50_bp"}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _pinned(path: str, digest: str) -> Path:
    _require(isinstance(path, str) and bool(path), "quality artifact path required")
    _require(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
             "quality artifact SHA256 required")
    file = Path(path)
    _require(file.is_absolute() and file.is_file(), "quality artifact must be an existing absolute path")
    _require(sha256(file) == digest, "quality artifact SHA256 mismatch")
    return file


def attach_checkm2(structural: dict, evidence_path: Path, *, expected_sha256: str) -> dict:
    """Return a copy with verified estimates; no assignment or scope inference.

    ``evidence_path`` is one JSON object from the recorded producer contract,
    with explicit manifest paths/columns. JSONL extraction must be recorded by
    its caller. Query bytes, native row, manifest membership and four structural
    metrics must agree. The execution record is a producer attestation whose
    hash is checked; this function cannot independently prove execution history.
    """
    _require(isinstance(structural, dict) and structural.get("schema_version") == STRUCTURAL_SCHEMA
             and structural.get("status") == "structurally_valid", "structural assembly assessment required")
    _require(structural.get("biological_quality", {}).get("status") == "unknown",
             "refusing to replace existing biological-quality evidence")
    evidence_file = _pinned(str(evidence_path.absolute()), expected_sha256)
    e = json.loads(evidence_file.read_text())
    _require(isinstance(e, dict) and e.get("schema_version") == SCHEMA, "unsupported measured-quality schema")
    _require((e.get("tool"), e.get("tool_version"), e.get("database_version"), e.get("database_sha256")) ==
             ("CheckM2", "1.1.0", "3", DATABASE_SHA256), "unsupported CheckM2 tool/database provenance")
    _require(e.get("evidence_origin") == "local_calculation" and e.get("status") == "estimated",
             "local measured-quality evidence required")
    _require(e.get("purity_claim") is False and e.get("residual_sequence_explained") is False
             and e.get("biological_quality_gate") == "not_defined"
             and e.get("taxonomic_assignment") == "not_evaluated", "quality estimates cannot assert a gate, purity, scope or assignment")
    query = _pinned(structural.get("query_path"), structural.get("query_sha256"))
    _require(e.get("query_sha256") == structural["query_sha256"], "quality evidence belongs to another query")
    files = {}
    for name in ("native_output", "input_manifest", "installed_files_inventory", "execution_record"):
        files[name] = _pinned(e.get(name + "_path"), e.get(name + "_sha256"))
    fields = e.get("input_manifest_fields")
    _require(isinstance(fields, dict) and set(fields) == {"identifier", "path", "sha256"}
             and all(isinstance(v, str) and bool(v) for v in fields.values())
             and len(set(fields.values())) == 3, "explicit distinct manifest columns required")
    manifest = read_tsv(files["input_manifest"], set(fields.values()))
    identifiers = [r[fields["identifier"]] for r in manifest]
    _require(bool(identifiers) and all(identifiers) and len(set(identifiers)) == len(identifiers),
             "quality manifest IDs must be nonempty and unique")
    key = e.get("input_id")
    _require(isinstance(key, str) and key in identifiers and e.get("native_row_key") == key,
             "quality row/manifest identity mismatch")
    member = next(r for r in manifest if r[fields["identifier"]] == key)
    _require(member[fields["sha256"]] == structural["query_sha256"]
             and Path(member[fields["path"]]).is_absolute()
             and Path(member[fields["path"]]).resolve() == query.resolve(), "manifest query path/hash mismatch")
    normalization_fields = {"normalization_record_path", "normalization_record_sha256",
                            "execution_input_manifest_path", "execution_input_manifest_sha256"}
    normalized_pin = None
    if normalization_fields & e.keys():
        _require(normalization_fields <= e.keys(), "incomplete normalization provenance")
        for name in ("normalization_record", "execution_input_manifest"):
            files[name] = _pinned(e[name + "_path"], e[name + "_sha256"])
        normalization = json.loads(files["normalization_record"].read_text())
        _require(isinstance(normalization, dict)
                 and normalization.get("schema_version") == "1.0.0-dev-checkm2-normalization"
                 and normalization.get("mode") == "ascii_universal_newlines_ensure_terminal_lf"
                 and normalization.get("mode_version") == "1", "unsupported normalization contract")
        original = _pinned(normalization.get("original_path"), normalization.get("original_sha256"))
        normalized = _pinned(normalization.get("normalized_path"), normalization.get("normalized_sha256"))
        normalized_pin = (normalized, normalization["normalized_sha256"])
        _require(original.resolve() == query.resolve() and normalization["original_sha256"] == e["query_sha256"],
                 "normalization original/query mismatch")
        for field in ("record_count", "total_bp"):
            _require(type(normalization.get(field)) is int and normalization[field] == structural["structural_metrics"][field],
                     "normalization structural metric mismatch")
        digest = hashlib.sha256()
        opener = gzip.open if query.suffix == ".gz" else open
        with opener(query, "rt", encoding="ascii", newline=None) as handle:
            for line in handle:
                digest.update(line.encode("ascii"))
                if not line.endswith("\n"):
                    digest.update(b"\n")
        _require(digest.hexdigest() == normalization["normalized_sha256"],
                 "normalized execution copy changes content beyond declared line endings")
        execution_inputs = read_tsv(files["execution_input_manifest"], {"input_id", "path", "sha256"})
        _require(len(execution_inputs) == 1 and execution_inputs[0]["input_id"] == key
                 and Path(execution_inputs[0]["path"]).is_absolute()
                 and Path(execution_inputs[0]["path"]).resolve() == normalized.resolve()
                 and execution_inputs[0]["sha256"] == normalization["normalized_sha256"],
                 "execution manifest/normalized query mismatch")
    native = read_tsv(files["native_output"], {"Name", "Completeness", "Contamination",
                      "Completeness_Model_Used", *NATIVE_STRUCTURE})
    _require(len(native) == len(manifest) and {r["Name"] for r in native} == set(identifiers),
             "native rows must match the complete frozen manifest exactly once")
    row = next(r for r in native if r["Name"] == key)
    _require(row == e.get("native_row"), "native quality row differs from supplied evidence")
    complete = explicit_percent(row["Completeness"], "CheckM2 completeness")
    try:
        contamination = float(row["Contamination"])
    except (ValueError, TypeError) as exc:
        raise ValueError("invalid native contamination estimate") from exc
    _require(math.isfinite(contamination) and contamination >= 0, "invalid native contamination estimate")
    # Native CheckM2 contamination is not capped at100. Values are estimates,
    # not the measured fraction of foreign bases in an assembly.
    _require(type(e.get("completeness_percent")) in (float, int)
             and e["completeness_percent"] == complete
             and type(e.get("contamination_percent")) in (float, int)
             and e["contamination_percent"] == contamination, "native/supplied quality estimates disagree")
    _require(bool(row["Completeness_Model_Used"]) and row["Completeness_Model_Used"] == e.get("completeness_model"),
             "completeness model provenance mismatch")
    for native_key, structural_key in NATIVE_STRUCTURE.items():
        _require(re.fullmatch(r"[0-9]+", row[native_key]) is not None
                 and int(row[native_key]) == structural["structural_metrics"][structural_key],
                 f"native/structural sequence mismatch: {native_key}")
    inventory = json.loads(files["installed_files_inventory"].read_text())
    _require(isinstance(inventory, list) and bool(inventory)
             and all(isinstance(r, dict) and r.get("path") and r.get("sha256") for r in inventory),
             "installed source/model inventory is missing")
    _require(len({r["path"] for r in inventory}) == len(inventory), "duplicate installed inventory paths")
    for record in inventory:
        _pinned(record["path"], record["sha256"])
    job = json.loads(files["execution_record"].read_text())
    _require(isinstance(job, dict) and job.get("status") == "passed" and job.get("rows") == len(manifest)
             and bool(job.get("finished_utc")), "quality execution is not recorded as complete")
    report_hash = job.get("report_sha256", job.get("native_report_sha256"))
    _require(report_hash == e["native_output_sha256"], "quality execution/output mismatch")
    if "report_sha256" in job and "native_report_sha256" in job:
        _require(job["report_sha256"] == job["native_report_sha256"], "conflicting quality execution report hashes")
    for field in ("input_manifest_sha256", "database_sha256", "normalization_record_sha256", "execution_input_manifest_sha256"):
        if field in job:
            _require(job[field] == e.get(field), f"quality execution/{field} mismatch")
    commands = job.get("commands")
    if commands is None:
        commands = [{"argv": job.get("argv"), "exit_status": job.get("exit_status")}]
    _require(isinstance(commands, list) and bool(commands) and all(isinstance(c, dict)
             and type(c.get("exit_status")) is int and c["exit_status"] == 0
             and isinstance(c.get("argv"), list) and bool(c["argv"])
             and all(isinstance(v, str) for v in c["argv"]) for c in commands),
             "successful quality command records required")
    _require(any("predict" in c["argv"] for c in commands), "quality prediction command missing")
    for name, path in files.items():
        _require(sha256(path) == e[name + "_sha256"], "quality artifact changed during import")
    if normalized_pin is not None:
        _require(sha256(normalized_pin[0]) == normalized_pin[1], "normalized query changed during import")
    _require(sha256(query) == e["query_sha256"] and sha256(evidence_file) == expected_sha256,
             "query or quality envelope changed during import")
    result = copy.deepcopy(structural)
    result["biological_quality"] = {
        "status": "estimated", "completeness_percent": complete,
        "contamination_percent": contamination, "quality_gate": "not_defined",
        "purity_established": False, "residual_sequence_explained": False,
        "evidence": [{"path": str(evidence_file), "sha256": expected_sha256,
                      "tool": e["tool"], "tool_version": e["tool_version"],
                      "database_version": e["database_version"], "database_sha256": e["database_sha256"],
                      "completeness_model": e["completeness_model"], "input_id": key,
                      "verification": "artifact_and_row_consistency_checked; producer_attested_execution; not_a_purity_test"}],
    }
    return result
