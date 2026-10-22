"""Disabled eligibility screen for the frozen candidate policy.

The screen reports whether saved evidence passes the candidate policy frozen on
6 October 2026 (DEC-005) for independent validation. It never assigns taxonomy,
enables production calling or establishes purity: the original assignment is
copied unchanged and hypothetical IDs are reported in separate fields.

It derives from the independently reviewed run-local prototype
(``results/tool-planning/2026-10-06/shadow-eligibility-astra``, SHA256
3efb34de0357c302982f9bdb196532340b0c97d5728802f7a570fb2a391e2d66). By the same
decision, query assembly quality is not part of the policy: CheckM2 is not run
on genomes being classified, and a supplied estimate is reported only. The
prototype's QC guides and QC reasons are therefore gone. A later decision
(6 October 2026) adds a report-only near-boundary note for a nearest single
species at 95% <= ANI < 96%; it is not a call and changes no state, so the call
policy is unchanged while the record and configuration schemas are version 2.1.
Envelope schemas are unchanged, so saved ANIm and warning envelopes remain valid. ``analyze`` adds a
``workflow_verification`` block where the file wrapper adds ``file_verification``.

``evaluate_eligibility`` is pure and caller-bound: a supplied digest is an
attestation, not proof that a file was read. ``EligibilityContext`` verifies
saved files, native QC artifacts and optional envelopes. The optional
reciprocal-ANIm and problem-warning envelopes bind an exact saved evidence-file
digest, so they are post-analysis inputs only; ``analyze`` runs the built-in
screen without them. No native program is executed here.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from decimal import Decimal, InvalidOperation
from pathlib import Path

from .evidence import sha256
from .rank import evaluate_rank, validate_rule
from .reporting import ReportingContext
from .scope import TAXON_FIELDS, query_digest

SCHEMA = "2.1-disabled-frozen-candidate-eligibility"
CONFIG_SCHEMA = "2.1-frozen-candidate-config"
SUPERSEDED_CONFIG_SCHEMAS = ("2.0-frozen-candidate-config",)
INVENTORY_SCHEMA = "1.0-shadow-primary-inventory"
DIAGNOSTICS_SCHEMA = "1.0-shadow-diagnostics"
WARNINGS_SCHEMA = "1.0-shadow-problem-warnings"
ELIGIBLE = "screen_eligible_under_explicit_test_rule"
CONFIG_FIELDS = ("schema_version", "policy_status", "decision_record", "production_enabled", "input_kind",
                 "rank_rule_canonical_sha256", "query_quality", "alert_ani_minimum", "alert_reference_af_minimum",
                 "near_boundary_decision_record", "near_boundary_ani_minimum", "near_boundary_both_af_minimum")
CONCEPT_FIELDS = ("species_id", "species_name", "subspecies_id", "subspecies_name", "lineage_id", "reference_role", "name_status")
ENVELOPE_BINDING = ("schema_version", "query_sha256", "evidence_sha256", "reference_package_sha256")
# The single frozen candidate point (Michael, 6 October 2026). These are not
# defaults: a configuration must state every value, and any other value is
# rejected rather than being evaluated under this record type.
POLICY_STATUS = "frozen_candidate_for_independent_validation"
DECISION_RECORD = "https://github.com/mbhall88/Enterobacter/issues/3#issuecomment-6014123942"
# Exact strings, as in the frozen file, so the frozen point has one canonical digest.
FROZEN_ALERT_POINT = {"alert_ani_minimum": "96", "alert_reference_af_minimum": "80"}
# Report-only band below the species call threshold (Michael, 6 October 2026).
NEAR_BOUNDARY_DECISION_RECORD = "https://github.com/mbhall88/Enterobacter/issues/3#issuecomment-6016074353"
FROZEN_NEAR_BOUNDARY_POINT = {"near_boundary_ani_minimum": "95", "near_boundary_both_af_minimum": "80"}
FROZEN_RANK_POINT = {"species_ani_percent": 96, "species_both_af_percent": 80, "subspecies_ani_percent": 98,
                     "subspecies_both_af_percent": 80, "subspecies_margin_percentage_points": 1}
LIMITATIONS = ["No biological assignment or production gate.",
               "Query assembly quality is not assessed: a supplied CheckM2 estimate is reported only, "
               "and moderate contamination or missing sequence is not detected.",
               "A near-boundary note names the closest species below the call threshold; it is not a call.",
               "No alert is not proof of purity or explanation of all unaligned sequence.",
               "Only rank-required missing contrasts block; other unreturned estimates remain unmeasured.",
               "Source names, MLST and genomovar cannot rescue eligibility."]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _digest(value: dict) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _is_hash(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch("[0-9a-f]{64}", value) is not None


def _number(value: object, name: str, *, maximum: int | None = 100) -> Decimal:
    _require(type(value) in (str, int, float), name + ": finite numeric value required")
    try:
        result = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError(name + ": invalid number") from exc
    _require(result.is_finite() and result >= 0 and (maximum is None or result <= maximum), name + ": invalid range")
    return result


def _exact_fields(obj: object, keys, label: str) -> None:
    _require(isinstance(obj, dict) and set(obj) == set(keys), label + ": exact fields required")


def validate_config(config: dict, rule: dict) -> None:
    """Accept only the disabled, frozen candidate configuration bound to this exact rule."""
    rule_values = validate_rule(rule)
    _require(not (isinstance(config, dict) and config.get("schema_version") in SUPERSEDED_CONFIG_SCHEMAS),
             "superseded configuration: schema " + CONFIG_SCHEMA + " adds the near-boundary report band")
    _exact_fields(config, CONFIG_FIELDS, "configuration")
    _require(config["schema_version"] == CONFIG_SCHEMA and config["production_enabled"] is False
             and config["input_kind"] == "assembly", "disabled assembly configuration required")
    _require(config["policy_status"] == POLICY_STATUS and config["decision_record"] == DECISION_RECORD,
             "frozen candidate policy record required")
    _require(config["query_quality"] == "reported_not_required", "query QC is reported only under this policy")
    _require(config["rank_rule_canonical_sha256"] == _digest(rule), "rank rule content mismatch")
    _require(all(config[k] == v for k, v in FROZEN_ALERT_POINT.items()),
             "this policy is frozen at alert96/refAF80 (exact string values)")
    _require(config["near_boundary_decision_record"] == NEAR_BOUNDARY_DECISION_RECORD
             and all(config[k] == v for k, v in FROZEN_NEAR_BOUNDARY_POINT.items()),
             "the near-boundary report band is frozen at ANI95/bothAF80 (exact string values)")
    _require(rule["production_enabled"] is False, "production rule rejected")
    _require(all(rule_values[k] == v for k, v in FROZEN_RANK_POINT.items()), "frozen rank point required")
    _require(rule["require_within_parent_competitor"] is True, "measured within-parent comparator required")


def _validate_core_objects(evidence: dict) -> None:
    """Reject malformed core containers before any adapter dereferences them."""
    _require(isinstance(evidence, dict), "evidence object required")
    for key in ("assignment", "provenance", "qc"):
        _require(isinstance(evidence.get(key), dict), key + ": object required")
    assignment = evidence["assignment"]
    _require(assignment.get("status") in ("not_evaluated", "unresolved", "unassigned")
             and isinstance(assignment.get("reason"), str) and bool(assignment["reason"].strip())
             and all(k in assignment and assignment[k] is None for k in TAXON_FIELDS),
             "unassigned input state with empty taxon fields required")
    provenance, qc = evidence["provenance"], evidence["qc"]
    for outer, key in ((provenance, "sha256"), (qc, "biological_quality"), (qc, "structural_metrics"),
                       (evidence, "scope_assessment"), (evidence, "reporting_policy")):
        if key in outer:
            _require(isinstance(outer[key], dict), key + ": object required")
    if "query_sha256" in provenance:
        _require(_is_hash(provenance["query_sha256"]), "invalid original-query digest")
    if "query" in provenance.get("sha256", {}):
        _require(_is_hash(provenance["sha256"]["query"]), "invalid original-query digest")
    ranked = evidence.get("ranked_reference_evidence")
    _require(isinstance(ranked, list) and all(isinstance(h, dict) for h in ranked), "ranked evidence objects required")
    biological = qc.get("biological_quality", {})
    if "evidence" in biological:
        _require(isinstance(biological["evidence"], list) and all(isinstance(x, dict) for x in biological["evidence"]),
                 "QC provenance must contain objects")
    query_digest(evidence)
    scope = _validate_scope(evidence)
    if assignment["reason"] == "unexplained_residual_sequence":
        _require(scope == "unexplained_after_checks", "unexplained assignment reason requires its bound scope veto")


def _validate_inventory(evidence: dict, inventory: dict) -> dict:
    _exact_fields(inventory, ("schema_version", "package_sha256", "component_id", "references", "representation"), "inventory")
    _require(inventory["schema_version"] == INVENTORY_SCHEMA and inventory["component_id"] == "primary_whole"
             and _is_hash(inventory["package_sha256"]), "invalid primary inventory")
    provenance = evidence["provenance"]
    _require(provenance.get("reference_package_sha256") == inventory["package_sha256"]
             and provenance.get("reference_component_id") == "primary_whole", "package/component mismatch")
    references = inventory["references"]
    _require(isinstance(references, list) and len(references) == evidence["reference_count"], "reference inventory count mismatch")
    _require(all(isinstance(r, dict) and all(isinstance(r.get(k), str) for k in ("reference_id", *CONCEPT_FIELDS))
                 for r in references), "inventory references need explicit string IDs and concept metadata")
    refs = {r["reference_id"]: r for r in references}
    _require(len(refs) == len(references), "duplicate inventory reference")
    for r in refs.values():
        _require(_is_hash(r.get("source_sha256")), "reference source digest required")
    for h in evidence["ranked_reference_evidence"]:
        _require(h["reference_id"] in refs, "reference outside primary inventory")
        for k in CONCEPT_FIELDS:
            _require(h[k] == refs[h["reference_id"]][k], "reference concept/name/role mismatch")
    _require(isinstance(inventory["representation"], dict), "representation required")
    for rank in ("species", "subspecies"):
        for cid in {r[rank + "_id"] for r in refs.values()} - {""}:
            x = inventory["representation"].get(cid)
            members = sorted(r["reference_id"] for r in refs.values() if r[rank + "_id"] == cid)
            _require(isinstance(x, dict) and x.get("rank") == rank and x.get("reference_assemblies") == len(members)
                     and isinstance(x.get("reference_ids"), list) and all(isinstance(m, str) for m in x["reference_ids"])
                     and sorted(x["reference_ids"]) == members, "representation count/member mismatch")
    return refs


def _validate_scope(evidence: dict) -> str:
    """Return the supplied scope outcome; an existing unexplained veto stays binding."""
    s = evidence.get("scope_assessment")
    if s is None:
        return "not_assessed"
    _require(s.get("outcome") in {"not_assessed", "explained_after_checks", "unexplained_after_checks"}, "invalid scope outcome")
    _require(s.get("query_sha256") == query_digest(evidence), "scope query mismatch")
    _require(s.get("schema_version") == "1.0.0-scope-assessment" and _is_hash(s.get("evidence_sha256"))
             and isinstance(s.get("reason"), str) and s["reason"].strip(), "scope provenance missing")
    checks = s.get("checks")
    _require(isinstance(checks, list) and bool(checks) == (s["outcome"] != "not_assessed"), "scope checks inconsistent")
    ids = set()
    for c in checks:
        _require(isinstance(c, dict) and c.get("status") == "completed"
                 and all(isinstance(c.get(k), str) and c[k].strip() for k in ("check_id", "method", "method_version", "artifact_path"))
                 and _is_hash(c.get("sha256")), "scope check provenance invalid")
        _require(c["check_id"] not in ids, "duplicate scope check")
        ids.add(c["check_id"])
    return s["outcome"]


def metric_triple(hit: dict) -> dict:
    return dict(ani_percent=hit["ani_percent"], af_query_percent=hit["aligned_fraction_query_percent"],
                af_reference_percent=hit["aligned_fraction_reference_percent"])


def _meets_species_gate(triple: dict, rule: dict) -> bool:
    return (_number(triple["ani_percent"], "ANI") >= _number(rule["species_ani_percent"], "speciesANI")
            and min(_number(triple["af_query_percent"], "queryAF"), _number(triple["af_reference_percent"], "referenceAF"))
            >= _number(rule["species_both_af_percent"], "speciesAF"))


def near_species_boundary(raw: dict | None, rule: dict, config: dict, species_names: dict) -> dict | None:
    """Report-only note: a nearest single species just below the species call ANI.

    Applies when the nearest reference(s) share one species, their ANI lies in
    [band minimum, species call ANI) and one of them has both aligned fractions
    at or above the band AF floor; the first such tied reference is the witness.
    It names that species but never makes or enables a call. ``species_names``
    maps reference IDs to species names; reports re-derive the note with it.
    """
    if raw is None:
        return None
    tied = raw["tied_references"]
    species_ids = {h["species_id"] for h in tied}
    ani = _number(raw["ani_percent"], "ANI")
    call_ani = _number(rule["species_ani_percent"], "speciesANI")
    band_ani = _number(config["near_boundary_ani_minimum"], "near_boundary_ani_minimum")
    band_af = _number(config["near_boundary_both_af_minimum"], "near_boundary_both_af_minimum")
    witnesses = [h for h in tied if min(_number(h["aligned_fraction_query_percent"], "queryAF"),
                                        _number(h["aligned_fraction_reference_percent"], "referenceAF")) >= band_af]
    if len(species_ids) != 1 or None in species_ids or not band_ani <= ani < call_ani or not witnesses:
        return None
    w = witnesses[0]
    name = species_names.get(w["reference_id"])
    _require(isinstance(name, str) and bool(name), "near-boundary witness species name required")
    return dict(state="closest_species_below_call_threshold", is_call=False, species_id=w["species_id"], species_name=name,
                reference_id=w["reference_id"], tied_reference_ids=[h["reference_id"] for h in tied],
                ani_percent=w["ani_percent"], af_query_percent=w["aligned_fraction_query_percent"],
                af_reference_percent=w["aligned_fraction_reference_percent"],
                band=dict(ani_minimum=config["near_boundary_ani_minimum"], species_call_ani=rule["species_ani_percent"],
                          both_af_minimum=config["near_boundary_both_af_minimum"]),
                interpretation="closest species, below the call threshold; not a call",
                message="closest species: %s at %s%% ANI, below the %s%% call threshold; not a call"
                        % (name, w["ani_percent"], rule["species_ani_percent"]))


def _bound_items(envelope: dict, items_key: str, schema: str, evidence: dict, evidence_sha256: str, label: str) -> list:
    """Return an envelope's items after checking it binds this query, evidence file and package."""
    _exact_fields(envelope, (*ENVELOPE_BINDING, items_key), label + " envelope")
    _require(envelope["schema_version"] == schema and envelope["query_sha256"] == query_digest(evidence)
             and envelope["evidence_sha256"] == evidence_sha256
             and envelope["reference_package_sha256"] == evidence["provenance"]["reference_package_sha256"],
             label + " binding mismatch")
    _require(isinstance(envelope[items_key], list), label + " list required")
    return envelope[items_key]


def _artifact_pins(artifacts: object, label: str) -> None:
    _require(isinstance(artifacts, list) and bool(artifacts), label + " source artifacts required")
    for a in artifacts:
        _exact_fields(a, ("path", "sha256"), label + " artifact")
        _require(isinstance(a["path"], str) and Path(a["path"]).is_absolute() and _is_hash(a["sha256"]),
                 "invalid " + label + " artifact pin")


def _diagnostic_results(envelope: dict | None, evidence: dict, rule: dict, refs: dict,
                        evidence_sha256: str, parent: str | None) -> list[dict]:
    """Same-reference skani versus reciprocal ANIm at the configured species gate.

    Any supplied, correctly bound gate disagreement blocks, including support
    for a competing parent. Neither estimator is treated as truth.
    """
    if envelope is None:
        return []
    items = _bound_items(envelope, "diagnostics", DIAGNOSTICS_SCHEMA, evidence, evidence_sha256, "diagnostic")
    hits = {h["reference_id"]: h for h in evidence["ranked_reference_evidence"]}
    results, ids = [], set()
    for d in items:
        _exact_fields(d, ("diagnostic_id", "method", "method_version", "reference_id", "reference_source_sha256",
                          "query_scope", "reference_scope", "skani_version", "skani", "directions", "artifacts"), "diagnostic")
        _require(isinstance(d["diagnostic_id"], str) and d["diagnostic_id"] and d["diagnostic_id"] not in ids,
                 "duplicate/invalid diagnostic ID")
        ids.add(d["diagnostic_id"])
        _require(d["method"] == "reciprocal_ANIm" and isinstance(d["method_version"], str) and d["method_version"].strip(),
                 "diagnostic method missing")
        _require(d["query_scope"] == "whole_assembly_as_supplied" and d["reference_scope"] == "whole_assembly_as_deposited"
                 and d["skani_version"] == evidence["provenance"]["engine_version"], "diagnostic scope/version mismatch")
        rid = d["reference_id"]
        _require(isinstance(rid, str) and rid in hits and d["reference_source_sha256"] == refs[rid]["source_sha256"],
                 "diagnostic reference/hash mismatch")
        h = hits[rid]
        _exact_fields(d["skani"], ("ani_percent", "af_query_percent", "af_reference_percent"), "skani triple")
        for k, v in metric_triple(h).items():
            _require(_number(d["skani"][k], k) == _number(v, k), "diagnostic skani triple mismatch")
        directions = d["directions"]
        _require(isinstance(directions, list) and len(directions) == 2 and all(isinstance(x, dict) for x in directions)
                 and all(isinstance(x.get("direction_id"), str) for x in directions)
                 and {x.get("direction_id") for x in directions} == {"reference_to_query", "query_to_reference"},
                 "two unique reciprocal directions required")
        skani_passes = _meets_species_gate(d["skani"], rule)
        observed = []
        for x in directions:
            _exact_fields(x, ("direction_id", "ani_percent", "af_query_percent", "af_reference_percent"), "method triple")
            observed.append(dict(**copy.deepcopy(x), meets_same_species_gate=_meets_species_gate(x, rule)))
        _artifact_pins(d["artifacts"], "diagnostic")
        disagreement = any(x["meets_same_species_gate"] != skani_passes for x in observed)
        results.append(dict(diagnostic_id=d["diagnostic_id"], reference_id=rid, method=d["method"],
                            method_version=d["method_version"], skani=copy.deepcopy(d["skani"]),
                            skani_meets_same_species_gate=skani_passes, directions=observed,
                            gate_disagreement=disagreement, relevant_to_nearest_parent=bool(parent) and h["species_id"] == parent,
                            blocks_shadow=disagreement, artifacts=copy.deepcopy(d["artifacts"]),
                            interpretation="method-dependent comparison; neither estimator is universal truth"))
    return results


def _problem_warnings(envelope: dict | None, evidence: dict, evidence_sha256: str) -> list[dict]:
    """Bound external concerns block; resolution is deliberately not implemented."""
    if envelope is None:
        return []
    items = _bound_items(envelope, "warnings", WARNINGS_SCHEMA, evidence, evidence_sha256, "external warning")
    out, ids = [], set()
    for w in items:
        _exact_fields(w, ("warning_id", "reason", "source", "source_version", "artifacts"), "external warning")
        _require(all(isinstance(w[k], str) and w[k].strip() for k in ("warning_id", "reason", "source", "source_version"))
                 and w["warning_id"] not in ids, "invalid/duplicate warning metadata")
        ids.add(w["warning_id"])
        _require(isinstance(w["artifacts"], list) and bool(w["artifacts"]), "external warning needs artifact pins")
        _artifact_pins(w["artifacts"], "warning")
        out.append(dict(**copy.deepcopy(w), interpretation="supplied problem warning; not measured detector success", blocks_shadow=True))
    return out


def _failed_search_record(evidence: dict) -> dict:
    return dict(schema_version=SCHEMA, production_enabled=False, evaluation_status="failed", eligibility_state="not_assessed",
                blocking_reasons=["native_search_failed"], would_conditionally_assign=None, near_species_boundary=None,
                assignment=copy.deepcopy(evidence.get("assignment")), purity_established=False, residual_sequence_explained=False)


def evaluate_eligibility(evidence: dict, rank_rule: dict, config: dict, *, reference_inventory: dict,
                         evidence_sha256: str, diagnostics: dict | None = None,
                         problem_warnings: dict | None = None) -> dict:
    """Return one disabled eligibility record; invalid contracts raise ValueError.

    States: ``failed`` (declared native search failure), ``incomplete`` (a
    rank-required parent contrast absent), or ``completed`` with
    ``screened_ineligible`` / ``screen_eligible_under_explicit_test_rule``.
    An established veto stays ineligible even when other evidence is missing.
    A biological QC estimate, if present, is checked for integrity and reported;
    it never changes eligibility.
    """
    _require(_is_hash(evidence_sha256), "source evidence SHA256 required")
    validate_config(config, rank_rule)
    _validate_core_objects(evidence)
    e = copy.deepcopy(evidence)
    if e.get("provenance", {}).get("search_exit_code") not in (None, 0):
        return _failed_search_record(e)
    rank = evaluate_rank(e, rank_rule)
    refs = _validate_inventory(e, reference_inventory)
    scope = _validate_scope(e)
    _require(e["provenance"].get("sequence_scope") == "whole_assembly_as_supplied"
             and e["provenance"].get("engine_version") == "0.3.2", "whole-input skani0.3.2 evidence required")
    alert_ani = _number(config["alert_ani_minimum"], "alert_ani_minimum")
    alert_reference_af = _number(config["alert_reference_af_minimum"], "alert_reference_af_minimum")

    biological = e["qc"].get("biological_quality", {})
    qc_status = biological.get("status", "unknown")
    _require(qc_status in {"estimated", "unknown", "missing", "not_assessed", "failed"}, "unrecognised biological QC status")
    reasons, incomplete, measurements = [], [], None
    if qc_status == "estimated":
        # Reported only; validated so a reported value belongs to this query.
        _require(e["qc"].get("query_sha256") == query_digest(e), "QC query hash required/mismatched")
        completeness = _number(biological.get("completeness_percent"), "completeness")
        # CheckM2 contamination can exceed 100%; it is a valid estimate, not a fraction of bases.
        contamination = _number(biological.get("contamination_percent"), "contamination", maximum=None)
        _require(biological.get("purity_established") is False and biological.get("residual_sequence_explained") is False,
                 "QC cannot claim purity or explain scope")
        measurements = dict(completeness_percent=float(completeness), contamination_percent=float(contamination),
                            estimate_can_exceed100=True)

    species = rank["species"]
    parent = species["candidate_id"]
    raw = species["raw_nearest"]
    nearest_parent = raw["tied_references"][0]["species_id"] if raw else None
    if not parent:
        # A missing parent contrast is incomplete evidence; other rank reasons are completed screens.
        if species["reason"] == "species_competitor_unmeasured":
            incomplete.append(species["reason"])
        else:
            reasons.append(species["reason"])

    # Fixed two-parent alert: single references of two distinct species, each
    # at the alert ANI and reference AF. Query AF has no floor here.
    triggers = [dict(reference_id=h["reference_id"], species_id=h["species_id"], **metric_triple(h))
                for h in e["ranked_reference_evidence"]
                if h["species_id"] and _number(h["ani_percent"], "ANI") >= alert_ani
                and _number(h["aligned_fraction_reference_percent"], "refAF") >= alert_reference_af]
    alert = len({h["species_id"] for h in triggers}) >= 2
    if alert:
        reasons.append("two_parent_support_alert")
    if scope == "unexplained_after_checks":
        reasons.append("unexplained_residual_sequence")
    method_results = _diagnostic_results(diagnostics, e, rank_rule, refs, evidence_sha256, parent or nearest_parent)
    method_block = any(d["blocks_shadow"] for d in method_results)
    if method_block:
        reasons.append("coverage_method_disagreement")
    warnings = _problem_warnings(problem_warnings, e, evidence_sha256)
    if warnings:
        reasons.append("external_problem_warning")

    if reasons:
        state = "screened_ineligible"
    elif incomplete:
        state = "not_assessed"
    else:
        state = ELIGIBLE
    status = "incomplete" if incomplete else "completed"
    screen = ("problem_unexplained_after_checks" if scope == "unexplained_after_checks" else
              "problem_flag_review_required" if alert or warnings or method_block else
              "incomplete" if incomplete else "screen_completed_no_problem_flag")
    near = near_species_boundary(raw, rank_rule, config, {h["reference_id"]: h["species_name"] for h in e["ranked_reference_evidence"]})
    _require(near is None or parent is None, "near-boundary note cannot accompany a species candidate")
    hypothetical = None
    if state == ELIGIBLE:
        hypothetical = dict(species_id=parent, species_name=species["candidate_name"],
                            subspecies_id=rank["subspecies"]["candidate_id"], subspecies_name=rank["subspecies"]["candidate_name"])
    represented = [reference_inventory["representation"][x] for x in (parent, rank["subspecies"]["candidate_id"]) if x]
    completed_checks = ["rank_and_inventory_invariants", "two_parent_support_hypothesis"]
    if diagnostics is not None:
        completed_checks.append("optional_bound_method_comparison")
    if problem_warnings is not None:
        completed_checks.append("supplied_problem_warning_binding")
    return dict(
        schema_version=SCHEMA, record_type="disabled_shadow_eligibility_not_assignment", production_enabled=False,
        evaluation_status=status, eligibility_state=state, blocking_reasons=list(dict.fromkeys(reasons + incomplete)),
        would_conditionally_assign=hypothetical, near_species_boundary=near,
        assignment=copy.deepcopy(e["assignment"]), rank_evidence=rank,
        QC=measurements, QC_status=qc_status, QC_role="reported_not_required", scope_screen_state=screen, supplied_scope_outcome=scope,
        purity_established=False, residual_sequence_explained=False,
        two_parent_alert=dict(triggered=alert, qualifying_same_reference_triples=triggers, query_AF_floor=None,
                              interpretation="fixed developmental alert; minor or unrepresented mixtures may escape"),
        method_diagnostics=method_results, external_problem_warnings=warnings,
        external_warning_status="supplied_bound_warnings" if problem_warnings is not None else "not_supplied_optional",
        warning_resolution="not_implemented; supplied warning cannot be cleared by old scope explanation",
        reference_representation=copy.deepcopy(represented), unmeasured_reference_count=e["not_returned_reference_count"],
        completed_checks=completed_checks,
        provenance=dict(evidence_file_sha256=evidence_sha256, canonical_evidence_sha256=_digest(evidence),
                        canonical_rule_sha256=_digest(rank_rule), canonical_config_sha256=_digest(config),
                        package_sha256=reference_inventory["package_sha256"], query_sha256=query_digest(e),
                        verification="pure_contract_checks; external_artifacts_and_origin_require_file_wrapper"),
        limitations=list(LIMITATIONS))


def primary_inventory(context: ReportingContext) -> dict:
    """Primary-component references and representation from a verified package."""
    refs = [dict(row, source_sha256=context.reference_source_sha256(row["reference_id"]))
            for row in context.evidence_references()]
    concepts = {r[k] for r in refs for k in ("species_id", "subspecies_id")} - {""}
    return dict(schema_version=INVENTORY_SCHEMA, package_sha256=context.package_sha256, component_id="primary_whole",
                references=refs, representation={cid: context.representation(cid) for cid in sorted(concepts)})


def _pinned_json(path: Path | str | None, expected: str | None) -> dict:
    _require(isinstance(path, (str, Path)) and bool(str(path)), "pinned file path required")
    p = Path(path).absolute()
    _require(_is_hash(expected) and p.is_file() and sha256(p) == expected, "file hash mismatch: " + str(p))
    return json.loads(p.read_text())


class EligibilityContext:
    """Verified primary package plus pinned rule/configuration for saved evidence.

    Reusable across many saved workflow outputs. It verifies query/package joins,
    any existing scope assessment, native QC artifacts and optional envelope
    artifacts by hash. Historical execution and derived ANIm metrics remain
    producer attestations; nothing is recalculated.
    """

    def __init__(self, package_path: Path, rank_rule_path: Path, config_path: Path, *,
                 expected_package_sha256: str, expected_rank_rule_sha256: str, expected_config_sha256: str):
        self.rule = _pinned_json(rank_rule_path, expected_rank_rule_sha256)
        self.config = _pinned_json(config_path, expected_config_sha256)
        validate_config(self.config, self.rule)
        self.reporting_context = ReportingContext(Path(package_path), expected_sha256=expected_package_sha256)
        self.inventory = primary_inventory(self.reporting_context)
        self._pins = [(Path(package_path), expected_package_sha256), (Path(rank_rule_path), expected_rank_rule_sha256),
                      (Path(config_path), expected_config_sha256)]

    def evaluate(self, evidence_path: Path, *, expected_evidence_sha256: str,
                 diagnostics_path: Path | None = None, expected_diagnostics_sha256: str | None = None,
                 problem_warnings_path: Path | None = None, expected_problem_warnings_sha256: str | None = None) -> dict:
        for path, digest in self._pins:
            _require(sha256(path) == digest, "context input changed")
        path = Path(evidence_path).absolute()
        e = _pinned_json(path, expected_evidence_sha256)
        _validate_core_objects(e)
        # Validates query bytes, complete primary joins, counts and any bound scope assessment.
        self.reporting_context.report(path)
        qc = e["qc"]
        biological = qc.get("biological_quality", {})
        if biological.get("status") == "estimated":
            saved = _pinned_json(e["provenance"].get("quality_evidence_path"), e["provenance"].get("quality_evidence_sha256"))
            _require(saved == qc, "embedded/saved QC disagreement")
            claims = biological.get("evidence")
            _require(isinstance(claims, list) and len(claims) == 1, "one pinned native QC envelope required")
            from .quality_evidence import attach_checkm2
            structural = copy.deepcopy(qc)
            structural["biological_quality"] = dict(status="unknown", completeness_percent=None, contamination_percent=None, evidence=[])
            rebuilt = attach_checkm2(structural, Path(claims[0]["path"]), expected_sha256=claims[0]["sha256"])
            _require(rebuilt == qc, "native QC replay disagrees")
        _require(bool(diagnostics_path) == bool(expected_diagnostics_sha256), "optional diagnostic path/hash must be paired")
        diagnostics = _pinned_json(diagnostics_path, expected_diagnostics_sha256) if diagnostics_path else None
        _require(bool(problem_warnings_path) == bool(expected_problem_warnings_sha256), "optional warning path/hash must be paired")
        warnings = _pinned_json(problem_warnings_path, expected_problem_warnings_sha256) if problem_warnings_path else None
        result = evaluate_eligibility(e, self.rule, self.config, reference_inventory=self.inventory,
                                      evidence_sha256=expected_evidence_sha256, diagnostics=diagnostics, problem_warnings=warnings)
        # The pure evaluation validated envelope schemas before artifact paths are read.
        for envelope, key in ((diagnostics, "diagnostics"), (warnings, "warnings")):
            if envelope is not None:
                for item in envelope[key]:
                    for a in item["artifacts"]:
                        _require(Path(a["path"]).is_file() and sha256(Path(a["path"])) == a["sha256"],
                                 "optional diagnostic/warning artifact hash mismatch")
        _require(sha256(path) == expected_evidence_sha256, "evidence changed during replay")
        result["file_verification"] = dict(
            evidence_sha256=expected_evidence_sha256, query_and_package_checked=True,
            QC_native_artifacts_checked=biological.get("status") == "estimated",
            diagnostic_artifacts_checked=diagnostics is not None, external_warning_artifacts_checked=warnings is not None,
            diagnostic_interpretation="producer_attested_derived_metrics; no_native_ANIm_recalculation",
            execution_history="hash_bound_producer_attestation_not_independent_proof")
        return result


def main(argv: list[str] | None = None) -> int:
    """Evaluate one saved evidence file; every failure record exits 2.

    Failure records are input/artifact contract failures; saved evidence from a
    failed native search is refused upstream as one. Completed and incomplete
    screens exit 0.
    """
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("evidence", "package", "rank-rule", "config", "output"):
        p.add_argument("--" + name, required=True, type=Path)
    for name in ("evidence", "package", "rank-rule", "config"):
        p.add_argument("--" + name + "-sha256", required=True)
    p.add_argument("--diagnostics", type=Path, help="Optional reciprocal-ANIm envelope bound to the saved evidence SHA256")
    p.add_argument("--diagnostics-sha256")
    p.add_argument("--problem-warnings", type=Path, help="Optional problem-warning envelope bound to the saved evidence SHA256")
    p.add_argument("--problem-warnings-sha256")
    a = p.parse_args(argv)
    if a.output.exists():
        p.error("refusing to overwrite output")
    try:
        context = EligibilityContext(a.package, a.rank_rule, a.config, expected_package_sha256=a.package_sha256,
                                     expected_rank_rule_sha256=a.rank_rule_sha256, expected_config_sha256=a.config_sha256)
        result = context.evaluate(a.evidence, expected_evidence_sha256=a.evidence_sha256,
                                  diagnostics_path=a.diagnostics, expected_diagnostics_sha256=a.diagnostics_sha256,
                                  problem_warnings_path=a.problem_warnings,
                                  expected_problem_warnings_sha256=a.problem_warnings_sha256)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        result = dict(schema_version=SCHEMA, production_enabled=False, evaluation_status="failed", eligibility_state="not_assessed",
                      blocking_reasons=["input_or_artifact_contract_failure"], error=str(exc), would_conditionally_assign=None,
                      near_species_boundary=None,
                      assignment=None, purity_established=False, residual_sequence_explained=False)
    with a.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return 0 if result["evaluation_status"] != "failed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
