"""Additive, policy-free reports from saved assembly evidence and pinned references.

No ANI/AF defaults, taxonomic assignments, quality assessments or tree runs are
introduced here. A development guide annotates an observation only. The existing
scope adapter must apply its approved unresolved-residual guard before reporting.
"""
from __future__ import annotations

import copy
import csv
import json
import re
from pathlib import Path

from .evidence import collect_evidence, explicit_percent, read_tsv, sha256
from .reference_package import canonical_json
from .package_runtime import load_runtime_package
from .scope import TAXON_FIELDS, query_digest, POLICY_ID

SCHEMA = "1.0.0-dev-uncertainty-report"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _pinned_file(record: dict) -> Path:
    _require(isinstance(record, dict) and set(record) == {"path", "sha256"},
             "report artifact requires exactly path and sha256")
    path = Path(record["path"])
    _require(path.is_absolute() and path.is_file(), "report artifact must be an existing absolute path")
    _require(sha256(path) == record["sha256"], "report artifact SHA256 mismatch")
    return path


def _check_assignment_scope(assignment: dict, scope: dict | None) -> None:
    """Shared output invariant; never infer a positive assignment or scope cause."""
    _require(isinstance(assignment, dict), "assignment must be an object")
    _require(assignment.get("status") in {"assigned", "unresolved", "not_evaluated"}, "unsupported assignment state")
    _require(isinstance(assignment.get("reason"), str) and bool(assignment["reason"].strip()), "assignment reason is required")
    if assignment["status"] != "assigned":
        _require(all(assignment.get(k) is None for k in TAXON_FIELDS), "unassigned state contains a positive taxon")
    else:
        _require(bool(assignment.get("species_id")) and bool(assignment.get("species_name")), "assigned state requires supplied species ID and name")
    if scope is not None:
        _require(isinstance(scope, dict), "scope assessment must be an object")
        if scope.get("outcome") == "unexplained_after_checks":
            _require(assignment["status"] == "unresolved" and assignment["reason"] == "unexplained_residual_sequence",
                     "apply the approved unresolved-residual guard before reporting")


class ReportingContext:
    """Verify one existing local package once; use its primary component only."""

    def __init__(self, manifest_path: Path, *, expected_sha256: str):
        path = _pinned_file({"path": str(manifest_path.absolute()), "sha256": expected_sha256})
        self._package_view = load_runtime_package(path, expected_sha256=expected_sha256)
        self._manifest = self._package_view.manifest
        self.validation = self._package_view.validation
        self.package_sha256 = expected_sha256
        self.package_path = str(path)
        self._concepts = {c["concept_id"]: c for c in self._manifest["concepts"]}
        self._references = {r["reference_id"]: r for r in self._manifest["references"]
                            if r["component_id"] == "primary_whole"}
        self._component = next(c for c in self._manifest["components"] if c["component_id"] == "primary_whole")
        table = self._manifest["recipe"]["source_tables"]["primary_references"]
        self._source_rows = {r["accession"]: r for r in read_tsv(Path(table["path"]), {"accession"})}
        _require(set(self._source_rows) == {r["accession"] for r in self._references.values()},
                 "primary source inventory does not match packaged references")

    def evidence_references(self) -> list[dict]:
        """Explicit adapter to the existing evidence schema, with pinned IDs."""
        rows = []
        for ref in self._references.values():
            species = self._concepts[ref["species_id"]]
            sub = self._concepts.get(ref["subspecies_id"])
            rows.append({"reference_id": ref["reference_id"], "search_path": ref["search_label"],
                         "species_id": species["concept_id"], "species_name": species["display_name"],
                         "subspecies_id": sub["concept_id"] if sub else "",
                         "subspecies_name": sub["display_name"] if sub else "", "lineage_id": "",
                         "name_status": (sub or species)["nomenclature_status"],
                         "reference_role": "type_anchor" if ref["type_anchor"] else "diversity"})
        return rows

    def reference_source_sha256(self, reference_id: str) -> str:
        """Original whole-reference digest recorded by the verified package."""
        return self._references[reference_id]["source"]["sha256"]

    def representation(self, concept_id: str) -> dict:
        """Count assemblies and recorded blocks, never independent isolates."""
        c = self._concepts[concept_id]
        _require(c["rank"] in {"species", "subspecies"}, "representation requires a named species/subspecies")
        refs = [r for r in self._references.values() if r[c["rank"] + "_id"] == concept_id]
        blocks = {}
        for field in ("clone_block", "split_block"):
            values = [self._source_rows[r["accession"]].get(field, "") for r in refs]
            blocks[field] = {"distinct_recorded": len(set(values) - {""}),
                             "assemblies_without_recorded_block": sum(not x for x in values),
                             "interpretation": "recorded correlation blocks; not isolate identity"}
        n = len(refs)
        caution = ("No primary reference assembly represents this concept." if not n else
                   "Within-group diversity is not characterised by one reference assembly." if n == 1 else
                   "Multiple development assemblies do not establish independent sampling or a characterised species boundary.")
        return {"concept_id": concept_id, "rank": c["rank"], "name": c["display_name"],
                "nomenclature_status": c["nomenclature_status"], "identifier_status": c["identifier_status"],
                "component_id": "primary_whole", "reference_assemblies": n,
                "reference_ids": sorted(r["reference_id"] for r in refs),
                "type_anchor_assemblies": sum(r["type_anchor"] for r in refs),
                "recorded_blocks": blocks, "sampling_limitation": caution,
                "source_qualifications": copy.deepcopy(c.get("source_annotation", {}))}

    def report(self, evidence_path: Path, *, development_guides: dict | None = None,
               tree_follow_up: dict | None = None) -> dict:
        """Read and verify saved evidence, preserving its assignment byte-values.

        Evidence provenance must pin this package/component and the query file.
        Historical execution remains attested, not rerun by the formatter.
        """
        evidence = json.loads(evidence_path.read_text())
        _require(evidence.get("record_type") == "development_evidence_only", "unsupported evidence record type")
        _require(evidence.get("input_kind") == "assembly", "read reporting adapter is not implemented")
        provenance = evidence.get("provenance", {})
        _require(provenance.get("reference_package_sha256") == self.package_sha256,
                 "evidence/reference package SHA256 mismatch")
        _require(provenance.get("reference_component_id") == "primary_whole", "wrong evidence reference component")
        _require(provenance.get("engine") == self._component["engine"] and
                 provenance.get("engine_version") == self._component["engine_version"], "ANI engine/version mismatch")
        query = Path(evidence.get("query_search_path", ""))
        _require(query.is_absolute() and query.is_file() and sha256(query) == query_digest(evidence),
                 "query path/hash mismatch; relocation is not supported")
        refs = {r["reference_id"]: r for r in self.evidence_references()}
        ranked = evidence.get("ranked_reference_evidence", [])
        _require(isinstance(ranked, list), "ranked evidence must be a list")
        raw = []
        for h in ranked:
            _require(h.get("reference_id") in refs, "evidence reference is outside primary component")
            ref = refs[h["reference_id"]]
            for field in ("species_id", "species_name", "subspecies_id", "subspecies_name", "lineage_id", "name_status", "reference_role"):
                _require(h.get(field) == ref[field], f"evidence/reference {field} mismatch")
            raw.append({"Query_file": str(query), "Ref_file": ref["search_path"], "ANI": h["ani_percent"],
                        "Align_fraction_query": h["aligned_fraction_query_percent"],
                        "Align_fraction_ref": h["aligned_fraction_reference_percent"]})
        checked = collect_evidence(raw, list(refs.values()), query_search_path=str(query), provenance=provenance)
        for key in ("reference_count", "unique_hit_count", "not_returned_reference_count", "nearest_reference_ids",
                    "ranked_reference_evidence", "evidence_status", "units"):
            _require(evidence.get(key) == checked[key], f"evidence {key} mismatch")
        assignment = evidence.get("assignment", {})
        _check_assignment_scope(assignment, evidence.get("scope_assessment"))
        status = assignment["status"]
        if status == "assigned":
            for rank in ("species", "subspecies"):
                cid, name = assignment.get(rank + "_id"), assignment.get(rank + "_name")
                _require(bool(cid) == bool(name), "assigned taxon ID/name incomplete")
                if cid:
                    c = self._concepts.get(cid, {})
                    _require(c.get("rank") == rank and c.get("display_name") == name, "assigned taxon registry mismatch")
                    if rank == "subspecies":
                        _require(c["parent_concept_id"] == assignment["species_id"], "assigned subspecies parent mismatch")
            # No genomovar caller exists in this interface; do not give an
            # unsupported, positive finer-rank assertion a validated appearance.
            _require(not assignment.get("genomovar_id") and not assignment.get("lineage_id"),
                     "positive lineage/genomovar reporting is not implemented")
        scope = evidence.get("scope_assessment")
        if scope is not None:
            _require(scope.get("query_sha256") == query_digest(evidence), "scope/query hash mismatch")
            _require(scope.get("outcome") in {"not_assessed", "explained_after_checks", "unexplained_after_checks"},
                     "unsupported scope state")
            _require(scope.get("schema_version") == "1.0.0-scope-assessment" and
                     isinstance(scope.get("reason"), str) and bool(scope["reason"].strip()), "scope provenance/reason incomplete")
            checks = scope.get("checks")
            _require(isinstance(checks, list) and bool(checks) == (scope["outcome"] != "not_assessed"), "scope check state mismatch")
            _require(all(isinstance(c, dict) for c in checks), "scope checks must be objects")
            _require(len({c.get("check_id") for c in checks}) == len(checks), "duplicate scope checks")
            for check in checks:
                _require(check.get("status") == "completed" and all(isinstance(check.get(k), str) and check[k].strip() for k in ("check_id", "method", "method_version")), "incomplete scope checks")
            if scope["outcome"] == "unexplained_after_checks":
                _require(evidence.get("reporting_policy", {}).get("id") == POLICY_ID, "unexplained scope reporting policy missing")
                _require(status == "unresolved" and assignment["reason"] == "unexplained_residual_sequence",
                         "apply the approved unresolved-residual guard before reporting")
            _require(isinstance(scope.get("evidence_sha256"), str) and
                     re.fullmatch(r"[0-9a-f]{64}", scope["evidence_sha256"]) is not None,
                     "scope original-evidence digest must be a lowercase SHA256")
            verified_checks = []
            for check in checks:
                path = _pinned_file({"path": check["artifact_path"], "sha256": check["sha256"]})
                verified_checks.append({k: check[k] for k in ("check_id", "method", "method_version", "status", "sha256")})
                verified_checks[-1]["artifact_path"] = str(path)
            assessment_path = scope.get("assessment_source_path")
            assessment_digest = scope.get("assessment_file_sha256")
            _require(bool(assessment_path) == bool(assessment_digest), "incomplete scope assessment-file pin")
            if assessment_path:
                _pinned_file({"path": assessment_path, "sha256": assessment_digest})
            # Interpretations are producer assertions. Verifying artifact bytes
            # cannot authenticate biology or bind a missing original evidence file.
            scope = {k: scope[k] for k in ("schema_version", "query_sha256", "evidence_sha256", "outcome", "reason")}
            scope.update(checks=verified_checks, artifact_hashes_verified=len(verified_checks),
                         assessment_source_path=assessment_path, assessment_file_sha256=assessment_digest,
                         verification="artifact_hashes_checked_by_reporting; producer_attested_interpretation; not_a_validated_scope_detector",
                         original_evidence_binding="digest_recorded_not_independently_bound_without_original_evidence_path")
        else:
            scope = {"outcome": "not_assessed", "reason": "No scope assessment supplied; coverage alone does not explain extra sequence."}
        nearest = [h for h in ranked if h["reference_id"] in evidence["nearest_reference_ids"]]
        competing = []
        species_seen = set()
        for h in ranked:
            if h["species_id"] not in species_seen:
                species_seen.add(h["species_id"])
                competing.append({"species_id": h["species_id"], "species_name": h["species_name"],
                                  "nomenclature_status": self._concepts[h["species_id"]]["nomenclature_status"],
                                  "best_reference_evidence": copy.deepcopy(h),
                                  "ani_difference_from_nearest_percentage_points": round(nearest[0]["ani_percent"] - h["ani_percent"], 12)})
        relevant = {h[k] for h in nearest for k in ("species_id", "subspecies_id") if h[k]}
        relevant.update(assignment[k] for k in ("species_id", "subspecies_id") if assignment.get(k))
        investigation = self._investigation(nearest, development_guides)
        follow_up = self._follow_up(tree_follow_up)
        report = {"schema_version": SCHEMA, "report_type": "development_uncertainty_report_no_calling_policy",
                "assignment": copy.deepcopy(assignment), "nearest_candidates": copy.deepcopy(nearest),
                "candidate_competing_taxa": competing,
                "ranked_reference_evidence": copy.deepcopy(ranked),
                "reference_representation": [self.representation(cid) for cid in sorted(relevant)],
                "evidence_status": evidence["evidence_status"],
                "missing_reference_estimates": evidence["not_returned_reference_count"],
                "scope_assessment": copy.deepcopy(scope),
                "subspecies_reason": assignment.get("subspecies_reason", "supplied_assignment" if assignment.get("subspecies_id") else "No subspecies assignment supplied; nearest subspecies evidence is not a rank decision."),
                "qc": copy.deepcopy(evidence.get("qc", {"status": "not_assessed"})),
                "mlst": copy.deepcopy(evidence.get("mlst", {"status": "not_run"})),
                "development_investigation": investigation, "optional_tree_follow_up": follow_up,
                "provenance": {"evidence_path": str(evidence_path.absolute()), "evidence_sha256": sha256(evidence_path),
                               "query_sha256": query_digest(evidence), "input_kind": "assembly",
                               "ani_engine": provenance["engine"], "ani_engine_version": provenance["engine_version"],
                               "query_sequence_scope": provenance.get("sequence_scope", "not_recorded"),
                               "reference_package_sha256": self.package_sha256, "reference_package_path": self.package_path,
                               "reference_package_id": self._manifest["package_id"], "reference_package_version": self._manifest["package_version"],
                               "reference_release_state": self._manifest["release_state"], "reference_component_id": "primary_whole",
                               "reference_sequence_scope": self._component["sequence_scope"], "reference_assemblies": len(refs),
                               "source_evidence_provenance": copy.deepcopy(provenance)},
                "limitations": ["This formatter preserves supplied assignments and never validates or creates a positive call.",
                                "Nearest candidates, source names and supplementary ST do not establish taxonomic assignments.",
                                "Genomovar nomination requires approved parent/scope and all-member compatibility; no positive genomovar reporting is implemented.",
                                "Absent estimates are unmeasured, not zero ANI; missing competitors cannot establish identity.",
                                "Assembly/block counts are not counts of independent isolates or confidence intervals.",
                                "Coverage can reflect query/reference content, divergence, method behaviour or quality; cause is not inferred here."]}

        report["uncertainty_messages"] = _uncertainty_messages(report)
        return report

    @staticmethod
    def _investigation(nearest: list[dict], guides: dict | None) -> dict:
        if guides is None:
            return {"status": "not_configured", "calling_policy_enabled": False, "observations": []}
        _require(set(guides) == {"boundary_ani_low_percent", "boundary_ani_high_percent", "coverage_af_percent"},
                 "development guides require explicit ANI interval and AF comparison; no defaults")
        g = {k: explicit_percent(v, k) for k, v in guides.items()}
        _require(g["boundary_ani_low_percent"] < g["boundary_ani_high_percent"], "invalid investigation interval")
        observations = []
        for h in nearest:
            if g["boundary_ani_low_percent"] <= h["ani_percent"] < g["boundary_ani_high_percent"]:
                observations.append({"reference_id": h["reference_id"], "kind": "configured_development_ANI_boundary",
                                     "message": "This nearest match lies in the configured investigation interval; it does not establish a confident species assignment."})
            directions = [d for d in ("query", "reference") if h[f"aligned_fraction_{d}_percent"] < g["coverage_af_percent"]]
            if directions:
                observations.append({"reference_id": h["reference_id"], "kind": "configured_development_AF_comparison",
                                     "limiting_directions": directions,
                                     "message": "The named direction is below the configured development AF comparison; this is not a production coverage requirement or a diagnosis of poor assembly quality."})
        return {"status": "explicitly_configured_development_only", "calling_policy_enabled": False,
                "guides": g, "interval_semantics": "lower_inclusive_upper_exclusive", "observations": observations}

    def _follow_up(self, pointer: dict | None) -> dict:
        if pointer is None:
            return {"status": "not_configured"}
        _require(set(pointer) == {"workflow", "species_reference_package_sha256", "reference_package", "sequence_scope", "paper_tree_relationship"},
                 "tree follow-up needs explicit workflow, package, sequence scope and paper-tree relationship")
        _require(pointer["species_reference_package_sha256"] == self.package_sha256, "tree pointer uses a different species package")
        _pinned_file(pointer["workflow"])
        _require(all(isinstance(pointer[k], str) and pointer[k].strip() for k in ("reference_package", "sequence_scope", "paper_tree_relationship")),
                 "tree follow-up scope fields are required")
        return {"status": "optional_investigation_not_run", **copy.deepcopy(pointer),
                "caution": "Inspect the named comparisons first. Placement weights are conditional on this tree/model; placement does not guarantee a formal name or resolve a rank question. Follow-up cannot promote a candidate or bypass unexplained scope."}


def _uncertainty_messages(report: dict) -> list[str]:
    assignment = report["assignment"]
    messages = []
    if assignment["status"] != "assigned":
        messages.append(f"Species assignment {assignment['status']}: {assignment['reason']}.")
    if not assignment.get("subspecies_id"):
        messages.append(report["subspecies_reason"])
    messages.append(report["scope_assessment"]["reason"])
    if len(report["nearest_candidates"]) > 1:
        messages.append(f"{len(report['nearest_candidates'])} references tie at nearest ANI; no tie winner is selected.")
    for observation in report["development_investigation"]["observations"]:
        directions = ", ".join(observation.get("limiting_directions", []))
        messages.append(f"{observation['reference_id']}: {observation['message']}" +
                        (f" Limiting direction(s): {directions}." if directions else ""))
    for representation in report["reference_representation"]:
        messages.append(f"{representation['name']}: {representation['sampling_limitation']} "
                        f"Name qualification: {representation['nomenclature_status']}.")
    return messages


def _supplementary_columns(report: dict) -> dict:
    """Expose typing without promoting a compatibility candidate to a call."""
    mlst = report.get("mlst", {"status": "not_run"})
    _require(isinstance(mlst, dict), "supplementary MLST must be an object")
    if mlst.get("status") == "not_run":
        _require(mlst in ({"status": "not_run"}, {"status": "not_run", "sequence_type": None}), "unexecuted MLST cannot contain typing claims")
        st = dict(mlst_state="not_run", reported_ST="", profile_compatible_candidate_ST="",
                  mlst_submission_cutoff="", mlst_database_limitation="")
    else:
        from .mlst import LIMITATION, CUTOFF
        _require(mlst.get("record_type") == "supplementary_mlst_evidence"
                 and mlst.get("query_sha256") == report["provenance"]["query_sha256"]
                 and mlst.get("taxonomic_assignment_performed") is False
                 and mlst.get("candidate_is_reported_ST") is False
                 and mlst.get("submission_cutoff") == CUTOFF and mlst.get("limitation") == LIMITATION,
                 "unsupported supplementary MLST binding or role")
        st_value = mlst.get("sequence_type")
        _require(mlst.get("state") in {"ambiguous_alleles", "incomplete_profile", "inexact_alleles", "exact_profile_match", "unmatched_profile"},
                 "unsupported supplementary MLST state")
        candidate = mlst.get("profile_compatible_candidate_ST")
        _require(candidate is None or (mlst.get("state") == "ambiguous_alleles" and st_value is None
                 and isinstance(candidate, str) and candidate.isdigit()), "candidate ST must remain ambiguous and unreported")
        _require(mlst.get("state") != "exact_profile_match" or st_value is not None, "exact profile match lacks ST")
        _require(st_value is None or (mlst.get("state") == "exact_profile_match"
                 and isinstance(st_value, str) and st_value.isdigit()), "only exact MLST profiles have reported STs")
        st = dict(mlst_state=mlst["state"], reported_ST=st_value or "",
                  profile_compatible_candidate_ST=mlst.get("profile_compatible_candidate_ST") or "",
                  mlst_submission_cutoff=mlst["submission_cutoff"], mlst_database_limitation=mlst["limitation"])
    gv = report.get("genomovar_compatibility")
    if gv is None:
        st.update(genomovar_compatibility_state="not_requested", compatible_genomovar_ids="", assigned_genomovar_id="")
    else:
        _require(isinstance(gv, dict) and gv.get("assigned_genomovar_id") is None
                 and gv.get("membership_policy_enabled") is False, "genomovar compatibility cannot contain an assignment")
        if gv.get("status") == "not_evaluated":
            _require(not gv.get("compatible_concept_ids") and gv.get("compatibility_result") is None, "unexecuted genomovar analysis contains compatibility claims")
        else:
            _require(gv.get("query_sha256") == report["provenance"]["query_sha256"]
                     and gv.get("reference_package_sha256") == report["provenance"]["reference_package_sha256"]
                     and gv.get("input_assignment") == report["assignment"], "genomovar upstream binding mismatch")
        st.update(genomovar_compatibility_state=gv.get("compatibility_result", gv.get("status")),
                  compatible_genomovar_ids=";".join(gv.get("compatible_concept_ids", [])), assigned_genomovar_id="")
    return st


def _eligibility_columns(report: dict) -> dict:
    """Hypothetical eligibility stays in its own columns; assignment columns are unchanged."""
    record = report.get("experimental_eligibility")
    if record is None:
        return dict(experimental_eligibility_status="not_requested", experimental_eligibility_state="",
                    experimental_eligibility_reasons="", experimental_scope_screen_state="",
                    hypothetical_species_id="", hypothetical_species="",
                    hypothetical_subspecies_id="", hypothetical_subspecies="",
                    near_boundary_species_id="", near_boundary_species="", near_boundary_ani_percent="", near_boundary_note="")
    from .eligibility import ELIGIBLE, SCHEMA as ELIGIBILITY_SCHEMA
    _require(isinstance(record, dict) and record.get("schema_version") == ELIGIBILITY_SCHEMA
             and record.get("production_enabled") is False and record.get("purity_established") is False
             and record.get("residual_sequence_explained") is False and record.get("assignment") == report["assignment"],
             "experimental eligibility must stay disabled and preserve the assignment")
    _require(isinstance(record.get("blocking_reasons"), list), "experimental eligibility reasons must be a list")
    rank = report.get("experimental_rank_candidate_evidence")
    _require(rank is None or "rank_evidence" not in record or record["rank_evidence"] == rank,
             "experimental eligibility and rank evidence disagree")
    hypothetical = record.get("would_conditionally_assign")
    if hypothetical is None:
        _require(record.get("eligibility_state") != ELIGIBLE, "eligible state requires separate hypothetical IDs")
        hypothetical = {}
    else:
        candidate = record["rank_evidence"]
        _require(isinstance(hypothetical, dict) and bool(hypothetical.get("species_id"))
                 and record.get("evaluation_status") == "completed" and record.get("eligibility_state") == ELIGIBLE
                 and not record["blocking_reasons"]
                 and hypothetical.get("species_id") == candidate["species"]["candidate_id"]
                 and hypothetical.get("subspecies_id") == candidate["subspecies"]["candidate_id"],
                 "hypothetical IDs require a completed unblocked screen matching the rank candidate")
    near = record.get("near_species_boundary")
    expected = None
    if record.get("rank_evidence") is not None:
        # Re-derive the note from the record's own rank evidence and rule, with names from this report.
        from .eligibility import FROZEN_NEAR_BOUNDARY_POINT, near_species_boundary
        rank = record["rank_evidence"]
        try:
            expected = near_species_boundary(rank["species"]["raw_nearest"], rank["rule"], FROZEN_NEAR_BOUNDARY_POINT,
                                             {h["reference_id"]: h["species_name"] for h in report["nearest_candidates"]})
        except ValueError:
            expected = "unverifiable"
    _require(near == expected and (near is None or (record.get("would_conditionally_assign") is None
                                                     and record["rank_evidence"]["species"]["candidate_id"] is None)),
             "a near-boundary note is report-only and must match the rank evidence")
    near = near or {}
    return dict(experimental_eligibility_status=record["evaluation_status"],
                experimental_eligibility_state=record["eligibility_state"],
                experimental_eligibility_reasons=";".join(record["blocking_reasons"]),
                experimental_scope_screen_state=record.get("scope_screen_state", ""),
                hypothetical_species_id=hypothetical.get("species_id") or "",
                hypothetical_species=hypothetical.get("species_name") or "",
                hypothetical_subspecies_id=hypothetical.get("subspecies_id") or "",
                hypothetical_subspecies=hypothetical.get("subspecies_name") or "",
                near_boundary_species_id=near.get("species_id", ""), near_boundary_species=near.get("species_name", ""),
                near_boundary_ani_percent=near.get("ani_percent", ""), near_boundary_note=near.get("message", ""))


def write_reports(reports: list[dict], *, json_path: Path, tsv_path: Path) -> None:
    """Same full report in JSON and TSV details; distinct visible label columns."""
    _require(not json_path.exists() and not tsv_path.exists() and json_path.absolute() != tsv_path.absolute(),
             "refusing to overwrite report output")
    rows = []
    for report in reports:
        _require(report.get("schema_version") == SCHEMA, "unsupported report schema")
        assignment = report["assignment"]
        _check_assignment_scope(assignment, report.get("scope_assessment"))
        nearest = report["nearest_candidates"]
        single = nearest[0] if len(nearest) == 1 else None
        representation = {r["concept_id"]: r for r in report["reference_representation"]}
        rows.append({**_supplementary_columns(report), **_eligibility_columns(report),
                     "assignment_status": assignment["status"], "assignment_reason": assignment["reason"],
                     "assigned_species_id": assignment.get("species_id") or "", "assigned_species": assignment.get("species_name") or "",
                     "assigned_subspecies_id": assignment.get("subspecies_id") or "", "assigned_subspecies": assignment.get("subspecies_name") or "",
                     "nearest_species_candidates": ";".join(sorted({h["species_name"] for h in report["nearest_candidates"]})),
                     "nearest_reference_ids": ";".join(h["reference_id"] for h in report["nearest_candidates"]),
                     "nearest_reference_count": len(nearest), "nearest_reference_tie": len(nearest) > 1,
                     "nearest_ani_percent": nearest[0]["ani_percent"] if nearest else "",
                     "single_nearest_af_query_percent": single["aligned_fraction_query_percent"] if single else "",
                     "single_nearest_af_reference_percent": single["aligned_fraction_reference_percent"] if single else "",
                     "single_nearest_species_reference_assemblies": representation[single["species_id"]]["reference_assemblies"] if single else "",
                     "single_nearest_subspecies_reference_assemblies": representation[single["subspecies_id"]]["reference_assemblies"] if single and single["subspecies_id"] else "",
                     "scope_outcome": report["scope_assessment"]["outcome"],
                     "uncertainty_messages": " | ".join(report["uncertainty_messages"]),
                     "query_sha256": report["provenance"]["query_sha256"], "ani_engine_version": report["provenance"]["ani_engine_version"],
                     "details_json": json.dumps(report, sort_keys=True, ensure_ascii=False, allow_nan=False)})
    fields = ["assignment_status", "assignment_reason", "assigned_species_id", "assigned_species", "assigned_subspecies_id",
              "assigned_subspecies", "nearest_species_candidates", "nearest_reference_ids", "nearest_reference_count", "nearest_reference_tie",
              "nearest_ani_percent", "single_nearest_af_query_percent", "single_nearest_af_reference_percent",
              "single_nearest_species_reference_assemblies", "single_nearest_subspecies_reference_assemblies",
              "scope_outcome", "uncertainty_messages", "query_sha256", "ani_engine_version",
              "mlst_state", "reported_ST", "profile_compatible_candidate_ST", "mlst_submission_cutoff",
              "mlst_database_limitation", "genomovar_compatibility_state", "compatible_genomovar_ids", "assigned_genomovar_id",
              "experimental_eligibility_status", "experimental_eligibility_state", "experimental_eligibility_reasons",
              "experimental_scope_screen_state", "hypothetical_species_id", "hypothetical_species",
              "hypothetical_subspecies_id", "hypothetical_subspecies", "near_boundary_species_id", "near_boundary_species",
              "near_boundary_ani_percent", "near_boundary_note", "details_json"]
    payload = canonical_json(reports)
    with json_path.open("x") as handle:
        handle.write(payload)
    with tsv_path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', required=True, type=Path)
    parser.add_argument('--reference-package', required=True, type=Path)
    parser.add_argument('--expected-package-sha256', required=True)
    parser.add_argument('--development-guides', type=Path, help='Optional explicit investigation values; never calling gates')
    parser.add_argument('--tree-follow-up', type=Path, help='Optional hash-pinned workflow pointer; does not run placement')
    parser.add_argument('--json-output', required=True, type=Path)
    parser.add_argument('--tsv-output', required=True, type=Path)
    args = parser.parse_args()
    try:
        context = ReportingContext(args.reference_package, expected_sha256=args.expected_package_sha256)
        report = context.report(args.evidence,
                                development_guides=json.loads(args.development_guides.read_text()) if args.development_guides else None,
                                tree_follow_up=json.loads(args.tree_follow_up.read_text()) if args.tree_follow_up else None)
        write_reports([report], json_path=args.json_output, tsv_path=args.tsv_output)
    except (ValueError, OSError, TypeError, KeyError) as exc:
        parser.exit(2, f'uncertainty reporting error: {exc}\n')


if __name__ == '__main__':
    main()
