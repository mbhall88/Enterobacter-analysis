"""Apply the approved unresolved-scope reporting invariant to saved evidence.

This is a reporting contract, not a scope detector or species caller. A supplied
assessment's artifacts are hash-checked; its biological interpretation remains
the responsibility of the producing method, which still needs validation.
"""
from __future__ import annotations

import argparse
import copy
import json
import re
from pathlib import Path

from .evidence import sha256

POLICY_ID = "unexplained-residual-unassigned-2026-10-06"
APPROVAL = "https://github.com/mbhall88/Enterobacter/issues/15#issuecomment-6004724669"
TAXON_FIELDS = ("species_id", "species_name", "subspecies_id", "subspecies_name", "lineage_id", "genomovar_id")


def query_digest(evidence: dict) -> str:
    """Recognise the two explicitly recorded development provenance layouts."""
    provenance = evidence.get("provenance", {})
    candidates = [v for v in [provenance.get("query_sha256"), provenance.get("sha256", {}).get("query")] if v]
    if not candidates or len(set(candidates)) != 1 or not re.fullmatch(r"[0-9a-f]{64}", candidates[0]):
        raise ValueError("evidence needs one unambiguous original-query SHA256")
    return candidates[0]


def apply_scope_assessment(evidence: dict, assessment: dict, *, evidence_sha256: str) -> dict:
    """Withhold taxonomy when scope checks explicitly end with unexplained DNA.

    The function never assigns a species or relaxes an ANI/AF threshold. A
    positive assessment cannot promote an unassigned candidate. Checks must be
    traceable; the CLI additionally verifies their local artifact digests.
    """
    if not isinstance(evidence, dict) or not isinstance(assessment, dict):
        raise ValueError("evidence and assessment must be JSON objects")
    if evidence.get("record_type") != "development_evidence_only":
        raise ValueError("unsupported evidence record type")
    if "scope_assessment" in evidence:
        raise ValueError("apply a fresh assessment to original evidence; do not overwrite assessment history")
    if assessment.get("schema_version") != "1.0.0-scope-assessment":
        raise ValueError("unsupported scope-assessment schema")
    if assessment.get("query_sha256") != query_digest(evidence):
        raise ValueError("scope assessment belongs to a different query")
    if not re.fullmatch(r"[0-9a-f]{64}", evidence_sha256) or assessment.get("evidence_sha256") != evidence_sha256:
        raise ValueError("scope assessment belongs to a different evidence record")
    outcome = assessment.get("outcome")
    if outcome not in {"not_assessed", "explained_after_checks", "unexplained_after_checks"}:
        raise ValueError("unsupported assessment outcome")
    assignment = evidence.get("assignment")
    if not isinstance(assignment, dict):
        raise ValueError("evidence needs an explicit assignment state")
    if outcome != "unexplained_after_checks" and any(assignment.get(k) is not None for k in TAXON_FIELDS):
        raise ValueError("evidence-only input cannot carry positive taxonomic assignments")
    checks = assessment.get("checks")
    if not isinstance(checks, list) or any(not isinstance(c, dict) for c in checks):
        raise ValueError("scope assessment needs an explicit check list")
    if outcome == "not_assessed" and checks:
        raise ValueError("not_assessed cannot claim completed checks")
    if outcome != "not_assessed" and not checks:
        raise ValueError("post-check assessment requires completed check provenance")
    if not isinstance(assessment.get("reason"), str) or not assessment["reason"].strip():
        raise ValueError("scope assessment needs a human-readable reason")
    identifiers = set()
    for check in checks:
        if any(not isinstance(check.get(k), str) or not check[k].strip()
               for k in ["check_id", "method", "method_version", "artifact_path", "sha256"]):
            raise ValueError("incomplete check provenance")
        if check["check_id"] in identifiers or check.get("status") != "completed":
            raise ValueError("duplicate or incomplete scope check")
        if not re.fullmatch(r"[0-9a-f]{64}", check["sha256"]):
            raise ValueError("scope check needs a valid artifact SHA256")
        identifiers.add(check["check_id"])
    result = copy.deepcopy(evidence)
    result["scope_assessment"] = {k: copy.deepcopy(assessment[k]) for k in [
        "schema_version", "query_sha256", "evidence_sha256", "outcome", "reason", "checks"]}
    # Verification facts belong to this adapter, never to caller metadata.
    result["scope_assessment"]["artifact_hashes_verified"] = 0
    result["scope_assessment"]["assessment_file_sha256"] = None
    result["scope_assessment"]["assessment_source_path"] = None
    result["scope_assessment"]["verification"] = "producer_attested_interpretation; not_a_validated_scope_detector"
    result["reporting_policy"] = {"id": POLICY_ID, "approval": APPROVAL,
                                  "numerical_thresholds_approved": False}
    if outcome == "unexplained_after_checks":
        result["assignment"] = {k: None for k in TAXON_FIELDS}
        result["assignment"].update(status="unresolved", reason="unexplained_residual_sequence")
        result["assignment"]["explanation"] = (
            "Additional sequence remains unexplained after the recorded checks. "
            "Nearest-reference candidates are evidence only; no species is assigned.")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", required=True, type=Path)
    parser.add_argument("--assessment", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite output")
    try:
        evidence = json.loads(args.evidence.read_text())
        assessment = json.loads(args.assessment.read_text())
        result = apply_scope_assessment(evidence, assessment, evidence_sha256=sha256(args.evidence))
        for check in result["scope_assessment"]["checks"]:
            source = Path(check["artifact_path"])
            if not source.is_absolute():
                source = args.assessment.parent / source
            if sha256(source) != check["sha256"]:
                raise ValueError(f"scope check artifact changed: {check['check_id']}")
            check["artifact_path"] = str(source.resolve())
        result["scope_assessment"]["artifact_hashes_verified"] = len(assessment["checks"])
        result["scope_assessment"]["assessment_file_sha256"] = sha256(args.assessment)
        result["scope_assessment"]["assessment_source_path"] = str(args.assessment.resolve())
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
    except (OSError, ValueError, TypeError, KeyError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
