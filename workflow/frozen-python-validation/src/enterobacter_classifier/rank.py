"""Disabled experimental rank evidence from validated assembly evidence.

``evaluate_rank(evidence, rule)`` is pure: it neither reads a reference package
nor verifies external artifacts. The caller must first validate the package and
its evidence join (as the development workflow/reporting context does). This
module rechecks record invariants, compares stable concept IDs, and reports
candidate evidence separately from the unchanged, unassigned input state.
Every numerical guide is required explicitly; none is a production default.
"""
from __future__ import annotations

import copy
from decimal import Decimal, InvalidOperation
import hashlib
import json
import re

from .evidence import SCHEMA_VERSION as EVIDENCE_SCHEMA, REFERENCE_FIELDS, validate_references
from .scope import TAXON_FIELDS, query_digest

SCHEMA_VERSION = "1.0.0-dev-experimental-rank"
RULE_SCHEMA_VERSION = "1.0-development-rank-proposal"
NUMERIC_FIELDS = (
    "species_ani_percent", "species_both_af_percent", "subspecies_ani_percent",
    "subspecies_both_af_percent", "subspecies_margin_percentage_points",
)
RULE_FIELDS = {
    "schema_version", "rule_id", "production_enabled", "input_kind",
    "require_within_parent_competitor", "margin_witness", "competitors",
    "scope_QC_eligibility", *NUMERIC_FIELDS,
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _number(value: object, name: str, *, positive: bool = False) -> Decimal:
    _require(type(value) in (str, int, float), f"{name}: explicit finite percentage required")
    try:
        result = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError(f"{name}: explicit finite percentage required") from exc
    _require(result.is_finite() and 0 <= result <= 100, f"{name}: value outside [0, 100]")
    _require(not positive or result > 0, f"{name}: returned ANI requires positive aligned fraction")
    return result


def _text(value: object) -> bool:
    return isinstance(value, str) and bool(value) and value == value.strip()


def _canonical_digest(value: dict) -> str:
    try:
        encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    except (ValueError, TypeError) as exc:
        raise ValueError("evidence and rule must contain finite JSON values") from exc
    return hashlib.sha256(encoded).hexdigest()


def validate_rule(rule: dict) -> dict[str, Decimal]:
    """Validate an explicit, complete disabled rule; return its decimal values."""
    _require(isinstance(rule, dict) and set(rule) == RULE_FIELDS,
             "explicit complete rank-rule object required; no numeric defaults")
    _require(rule["schema_version"] == RULE_SCHEMA_VERSION and rule["production_enabled"] is False,
             "development-only rule; production must remain disabled")
    _require(rule["input_kind"] == "assembly", "read rank adapter is not implemented")
    _require(_text(rule["rule_id"]), "explicit rule identifier required")
    _require(type(rule["require_within_parent_competitor"]) is bool,
             "within-parent requirement must be an explicit boolean")
    _require(rule["margin_witness"] == "same_reference_meets_ANI_AF_and_margin" and
             rule["competitors"] == "all_raw_other_groups_including_foreign_parents",
             "same-witness and all-raw-competitor semantics are required")
    _require(rule["scope_QC_eligibility"] == "separate_not_established_here",
             "rank evidence cannot establish QC or scope eligibility")
    return {key: _number(rule[key], key) for key in NUMERIC_FIELDS}


def _parent_key(hit: dict) -> tuple[str, str]:
    if hit["species_id"]:
        return "species", hit["species_id"]
    if hit["lineage_id"]:
        return "lineage", hit["lineage_id"]
    # Unmapped contextual references are not silently merged into one taxon.
    return "untyped_reference", hit["reference_id"]


def _group_key(hit: dict) -> tuple[str, str]:
    if hit["subspecies_id"]:
        return "subspecies", hit["subspecies_id"]
    if hit["species_id"]:
        return "species_only", hit["species_id"]
    return _parent_key(hit)


def _group(hit: dict) -> dict:
    kind, identifier = _group_key(hit)
    return {"kind": kind, "concept_id": None if kind == "untyped_reference" else identifier,
            "species_id": hit["species_id"] or None, "subspecies_id": hit["subspecies_id"] or None,
            "name": hit["subspecies_name"] or hit["species_name"] or None}


def _hit(hit: dict | None) -> dict | None:
    if hit is None:
        return None
    return {"reference_id": hit["reference_id"], "species_id": hit["species_id"] or None,
            "subspecies_id": hit["subspecies_id"] or None, "lineage_id": hit["lineage_id"] or None,
            "reference_role": hit["reference_role"], "reference_name_status": hit["name_status"],
            "ani_percent": float(hit["ani"]), "aligned_fraction_query_percent": float(hit["aq"]),
            "aligned_fraction_reference_percent": float(hit["ar"])}


def _best(hits: list[dict]) -> dict | None:
    if not hits:
        return None
    top = hits[0]["ani"]
    return {"ani_percent": float(top), "tied_references": [_hit(h) for h in hits if h["ani"] == top]}


def _normalize(evidence: dict) -> list[dict]:
    _require(isinstance(evidence, dict), "evidence must be an object")
    _require(evidence.get("schema_version") == EVIDENCE_SCHEMA and
             evidence.get("record_type") == "development_evidence_only", "unsupported evidence schema/type")
    _require(evidence.get("input_kind") == "assembly", "read rank adapter is not implemented")
    _require(evidence.get("units") == {"ani": "percent", "aligned_fraction": "percent"},
             "explicit ANI and AF percent units required")
    _require(_text(evidence.get("query_search_path")), "explicit query search path required")
    provenance = evidence.get("provenance")
    _require(isinstance(provenance, dict), "evidence provenance required")
    _require(provenance.get("input_kind") == "assembly" and
             type(provenance.get("search_exit_code")) is int and provenance["search_exit_code"] == 0,
             "only successful assembly-search evidence is accepted")
    _require(_text(provenance.get("engine_version")) and _text(provenance.get("search_parameters_record")),
             "engine version and search-parameter provenance required")
    if "reference_component_id" in provenance:
        _require(provenance["reference_component_id"] == "primary_whole", "rank evidence requires the primary whole component")
    if "reference_package_sha256" in provenance:
        value = provenance["reference_package_sha256"]
        _require(isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None,
                 "invalid reference package digest")
    assignment = evidence.get("assignment")
    _require(isinstance(assignment, dict) and assignment.get("status") in ("not_evaluated", "unresolved", "unassigned"),
             "rank evidence requires an unassigned input state")
    _require(_text(assignment.get("reason")) and all(k in assignment and assignment[k] is None for k in TAXON_FIELDS),
             "unassigned input must have a reason and empty taxon fields")
    qc = evidence.get("qc")
    _require(isinstance(qc, dict) and _text(qc.get("status")), "separate supplied QC state required")
    biological = qc.get("biological_quality")
    if "biological_quality" in qc:
        _require(isinstance(biological, dict) and _text(biological.get("status")), "invalid biological-quality state")
    if "query_sha256" in qc:
        _require(qc["query_sha256"] == query_digest(evidence), "QC/evidence query digest mismatch")
    scope = evidence.get("scope_assessment")
    if "scope_assessment" in evidence:
        _require(isinstance(scope, dict) and scope.get("outcome") in
                 ("not_assessed", "explained_after_checks", "unexplained_after_checks"), "invalid supplied scope state")
        _require(scope.get("query_sha256") == query_digest(evidence), "scope/evidence query digest mismatch")
        if scope["outcome"] == "unexplained_after_checks":
            _require(assignment["status"] == "unresolved" and assignment["reason"] == "unexplained_residual_sequence",
                     "apply the approved unresolved-residual guard before rank evaluation")
    original = evidence.get("ranked_reference_evidence")
    _require(isinstance(original, list), "ranked evidence must be a list")
    normalized, refs, ranks = [], [], {}
    for h in original:
        _require(isinstance(h, dict), "reference evidence must be an object")
        metadata = {k: h.get(k) for k in REFERENCE_FIELDS if k != "search_path"}
        _require(all(isinstance(v, str) for v in metadata.values()), "explicit reference/concept metadata required")
        metadata["search_path"] = metadata["reference_id"]
        refs.append(metadata)
        for rank in ("species", "subspecies", "lineage"):
            cid = h[rank + "_id"]
            if cid:
                _require(cid not in ranks or ranks[cid] == rank, "one concept ID cannot represent multiple ranks")
                ranks[cid] = rank
        n = dict(h, ani=_number(h.get("ani_percent"), "ANI"),
                 aq=_number(h.get("aligned_fraction_query_percent"), "AF query", positive=True),
                 ar=_number(h.get("aligned_fraction_reference_percent"), "AF reference", positive=True))
        _require(_number(h.get("minimum_directional_af_percent"), "minimum AF") == min(n["aq"], n["ar"]),
                 "minimum aligned fraction mismatch")
        normalized.append(n)
    if refs:
        validate_references(refs)
    hits = sorted(normalized, key=lambda h: (-h["ani"], h["reference_id"]))
    _require([h["reference_id"] for h in original] == [h["reference_id"] for h in hits], "evidence is not canonically ranked")
    for key in ("reference_count", "raw_hit_count", "unique_hit_count", "not_returned_reference_count", "duplicate_rows_collapsed"):
        _require(type(evidence.get(key)) is int and evidence[key] >= 0, f"invalid evidence count: {key}")
    _require(evidence["reference_count"] > 0 and evidence["unique_hit_count"] == len(hits) and
             evidence["reference_count"] - len(hits) == evidence["not_returned_reference_count"] and
             evidence["raw_hit_count"] - len(hits) == evidence["duplicate_rows_collapsed"], "inconsistent evidence counts")
    nearest = [h["reference_id"] for h in hits if h["ani"] == hits[0]["ani"]] if hits else []
    _require(evidence.get("nearest_reference_ids") == nearest, "nearest reference IDs mismatch")
    _require(evidence.get("evidence_status") == ("estimates_returned" if hits else "no_estimate_returned"),
             "evidence status mismatch")
    return hits


def _passes(hit: dict, ani: Decimal, af: Decimal) -> bool:
    return hit["ani"] >= ani and min(hit["aq"], hit["ar"]) >= af


def _species(hits: list[dict], p: dict) -> tuple[dict, str | None]:
    out = {"candidate_id": None, "candidate_name": None, "reason": "no_estimate",
           "raw_nearest": _best(hits), "qualifying_witness": None, "best_other_parent": None,
           "competing_parents": []}
    if not hits:
        return out, None
    top = hits[0]
    other = [h for h in hits if _parent_key(h) != _parent_key(top)]
    out["best_other_parent"] = _best(other)
    parents = {}
    for h in other:
        parents.setdefault(_parent_key(h), []).append(h)
    out["competing_parents"] = [dict(kind=key[0], concept_id=key[1] if key[0] != "untyped_reference" else None,
                                      returned_reference_count=len(hs), best=_best(hs))
                                for key, hs in parents.items()]
    if len({_parent_key(h) for h in hits if h["ani"] == top["ani"]}) > 1:
        out["reason"] = "different_species_nearest_tie"
    elif not top["species_id"]:
        out["reason"] = "nearest_reference_has_no_species"
    elif not other:
        out["reason"] = "species_competitor_unmeasured"
    elif top["ani"] < p["species_ani_percent"]:
        out["reason"] = "nearest_below_species_ANI"
    else:
        eligible = [h for h in hits if h["species_id"] == top["species_id"] and
                    _passes(h, p["species_ani_percent"], p["species_both_af_percent"])]
        if not eligible:
            out["reason"] = "species_coverage_not_met"
        else:
            out.update(candidate_id=top["species_id"], candidate_name=top["species_name"],
                       reason="diagnostic_candidate_only", qualifying_witness=_hit(eligible[0]))
    return out, out["candidate_id"]


def _subspecies(hits: list[dict], parent: str | None, p: dict, rule: dict) -> dict:
    out = {"candidate_id": None, "candidate_name": None, "reason": "parent_species_unresolved",
           "raw_nearest_group": None, "raw_nearest_ani_percent": None,
           "ani_af_witness": None, "qualifying_witness": None,
           "best_within_parent_competitor": None, "best_global_group_competitor": None,
           "best_foreign_parent_competitor": None, "raw_global_gap_percentage_points": None,
           "raw_within_parent_gap_percentage_points": None, "qualifying_witness_gap_percentage_points": None,
           "competing_groups": []}
    if parent is None:
        return out
    ph = [h for h in hits if h["species_id"] == parent]
    top = ph[0]
    out.update(raw_nearest_group=_group(top), raw_nearest_ani_percent=float(top["ani"]))
    # Retain lower competitors even when the nearest group is ambiguous or
    # species-only. These are observations, not a chosen subgroup or margin.
    other = [h for h in hits if _group_key(h) != _group_key(top)]
    groups = {}
    for h in other:
        groups.setdefault(_group_key(h), []).append(h)
    out["competing_groups"] = [dict(group=_group(hs[0]), returned_reference_count=len(hs), best=_best(hs))
                               for hs in groups.values()]
    if len({_group_key(h) for h in ph if h["ani"] == top["ani"]}) > 1:
        out["reason"] = "within_parent_nearest_group_tie"
        return out
    if not top["subspecies_id"]:
        out["reason"] = "nearest_reference_has_no_subspecies"
        return out
    within = [h for h in ph if _group_key(h) != _group_key(top)]
    foreign = [h for h in hits if h["species_id"] != parent]
    out.update(best_within_parent_competitor=_best(within), best_global_group_competitor=_best(other),
               best_foreign_parent_competitor=_best(foreign))
    if within:
        out["raw_within_parent_gap_percentage_points"] = float(top["ani"] - within[0]["ani"])
    if other:
        out["raw_global_gap_percentage_points"] = float(top["ani"] - other[0]["ani"])
    if top["ani"] < p["subspecies_ani_percent"]:
        out["reason"] = "nearest_group_below_subspecies_ANI"
        return out
    eligible = [h for h in ph if _group_key(h) == _group_key(top) and
                _passes(h, p["subspecies_ani_percent"], p["subspecies_both_af_percent"])]
    if not eligible:
        out["reason"] = "subspecies_coverage_not_met"
        return out
    witness = eligible[0]
    out["ani_af_witness"] = _hit(witness)
    if other:
        out["qualifying_witness_gap_percentage_points"] = float(witness["ani"] - other[0]["ani"])
    if rule["require_within_parent_competitor"] and not within:
        out["reason"] = "within_parent_competitor_unmeasured"
    elif not other:
        out["reason"] = "global_group_competitor_unmeasured"
    elif top["ani"] - other[0]["ani"] < p["subspecies_margin_percentage_points"]:
        out["reason"] = "global_group_margin_below_guide"
    elif witness["ani"] - other[0]["ani"] < p["subspecies_margin_percentage_points"]:
        out["reason"] = "qualifying_reference_margin_below_guide"
    else:
        out.update(candidate_id=top["subspecies_id"], candidate_name=top["subspecies_name"],
                   reason="diagnostic_candidate_only", qualifying_witness=_hit(witness))
    return out


def evaluate_rank(evidence: dict, rule: dict) -> dict:
    """Return experimental candidate hierarchy, preserving unassigned state.

    Concept IDs define identity. Names are checked for consistent presentation,
    never used to group references. Species-only references remain no-subgroup
    competitors; foreign-parent and low-AF estimates remain competitors too.
    Exact ANI ties are retained. QC/scope status is reported as supplied, never
    promoted by a taxonomic match or independently validated by this function.
    """
    p = validate_rule(rule)
    hits = _normalize(evidence)
    evidence_digest, rule_digest = _canonical_digest(evidence), _canonical_digest(rule)
    species, parent = _species(hits, p)
    subspecies = _subspecies(hits, parent, p, rule)
    scope = evidence.get("scope_assessment", {})
    qc = evidence["qc"]
    return {
        "schema_version": SCHEMA_VERSION, "record_type": "experimental_rank_candidate_evidence_only",
        "rule_id": rule["rule_id"], "production_enabled": False, "rule": copy.deepcopy(rule),
        "canonical_rule_sha256": rule_digest, "canonical_input_evidence_sha256": evidence_digest,
        "hash_semantics": "canonical_JSON_content; not_source_file_hashes_or_external_artifact_verification",
        "input_kind": "assembly", "query_search_path": evidence["query_search_path"],
        "engine_version": evidence["provenance"]["engine_version"],
        "reference_package_sha256": evidence["provenance"].get("reference_package_sha256"),
        "reference_component_id": evidence["provenance"].get("reference_component_id"),
        "species": species, "subspecies": subspecies,
        "candidate_state": "subspecies_candidate_evidence" if subspecies["candidate_id"] else
                           "species_only_candidate_evidence" if parent else "parent_unresolved",
        "assignment": copy.deepcopy(evidence["assignment"]),
        "QC_scope_eligibility": {
            "status": "separate_not_established_here", "supplied_qc_status": qc["status"],
            "supplied_biological_quality_status": qc.get("biological_quality", {}).get("status", "not_supplied"),
            "supplied_scope_outcome": scope.get("outcome", "not_assessed"),
            "verification": "input_statuses_only; no_QC_scope_artifact_or_biological_validation_by_rank",
        },
        "limitations": [
            "Candidate evidence is not an assignment, accuracy estimate or adopted numerical policy.",
            "Missing-reference estimates are censored, not zero ANI; absent competitors are not infinite margins.",
            "The caller validates reference package/file provenance; this pure function validates record invariants only.",
            "QC/scope prerequisites remain separate; favorable rank evidence cannot explain a problematic residual.",
        ],
    }
