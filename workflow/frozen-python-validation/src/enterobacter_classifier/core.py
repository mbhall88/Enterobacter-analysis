"""Deterministic evidence validation and provisional ANI assignment.

This module intentionally makes no biological naming decision.  It records the
nearest supported references and emits an unresolved/provisional result until a
held-out evaluation supplies a validated policy.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable

from . import OUTPUT_SCHEMA_VERSION


@dataclass(frozen=True)
class Reference:
    reference_id: str
    species: str
    subspecies: str = ""
    operational_cluster_id: str = ""


def _number(value: Any) -> float | None:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def percent(value: Any) -> float | None:
    """Accept a fraction or percentage and reject values outside both domains."""
    value = _number(value)
    if value is None or value < 0:
        return None
    if value <= 1:
        return value * 100
    if value <= 100:
        return value
    return None


def classify_hits(
    hits: Iterable[dict[str, Any]], references: dict[str, Reference], config: dict[str, Any]
) -> dict[str, Any]:
    """Validate hits and produce a traceable, non-production assignment record."""
    threshold = _number(config.get("species_ani_threshold"))
    coverage_threshold = _number(config.get("minimum_alignment_fraction_percent"))
    if threshold is None or coverage_threshold is None:
        raise ValueError("config requires finite species and coverage thresholds")

    observed = list(hits)
    valid: list[dict[str, Any]] = []
    rejected: list[dict[str, str]] = []
    for hit in observed:
        ref_id = str(hit.get("reference_id", "")).strip()
        ani = _number(hit.get("ani"))
        af_ref = percent(hit.get("af_ref"))
        af_query = percent(hit.get("af_query"))
        if ref_id not in references:
            rejected.append({"reference_id": ref_id, "reason": "reference_absent_from_database"})
        elif ani is None or ani < 0 or ani > 100:
            rejected.append({"reference_id": ref_id, "reason": "invalid_ani"})
        elif af_ref is None or af_query is None:
            rejected.append({"reference_id": ref_id, "reason": "missing_or_invalid_reciprocal_coverage"})
        elif min(af_ref, af_query) < coverage_threshold:
            rejected.append({"reference_id": ref_id, "reason": "below_reciprocal_coverage_threshold"})
        else:
            valid.append({"reference_id": ref_id, "ani": ani, "af_ref": af_ref, "af_query": af_query})

    valid.sort(key=lambda row: (-row["ani"], row["reference_id"]))
    base = {
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "policy_status": config.get("policy_status", "unspecified"),
        "thresholds": {"species_ani": threshold, "minimum_reciprocal_coverage_percent": coverage_threshold},
        "n_hits_observed": len(observed), "n_hits_accepted": len(valid), "rejected_hits": rejected,
        "nearest_reference": None, "competing_reference": None,
        "species": "unresolved", "subspecies": "unresolved", "operational_cluster_id": "",
        "assignment_status": "unresolved", "abstention_reason": "",
    }
    if not observed:
        base["abstention_reason"] = "no_hits"
        return base
    if not valid:
        base["abstention_reason"] = "no_valid_reciprocal_evidence"
        return base

    def trace(row: dict[str, Any]) -> dict[str, Any]:
        ref = references[row["reference_id"]]
        return {**row, "species": ref.species, "subspecies": ref.subspecies,
                "operational_cluster_id": ref.operational_cluster_id}

    nearest = trace(valid[0])
    competitors = [trace(row) for row in valid[1:] if references[row["reference_id"]].species != nearest["species"]]
    base["nearest_reference"] = nearest
    base["competing_reference"] = competitors[0] if competitors else None
    tied_species = {trace(row)["species"] for row in valid if row["ani"] == nearest["ani"]}
    if nearest["ani"] < threshold:
        base["abstention_reason"] = "below_provisional_species_ani_threshold"
        return base
    if len(tied_species) > 1:
        base["abstention_reason"] = "cross_species_tie_at_nearest_ani"
        return base

    base.update({"species": nearest["species"], "subspecies": nearest["subspecies"],
                 "operational_cluster_id": nearest["operational_cluster_id"]})
    if config.get("allow_confident_assignments", False):
        base["assignment_status"] = "provisional_candidate"
    else:
        base["abstention_reason"] = "policy_not_validated"
    return base
