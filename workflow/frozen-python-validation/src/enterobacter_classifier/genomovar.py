"""Compatibility with frozen ANI genomovars, never new clustering or assignment.

The historical chromosome reference definition is retained for whole queries.
Missing comparisons and multiple compatible groups remain explicit. A medoid
match is not a substitute for comparisons with every definition member.
"""
from __future__ import annotations

import copy
from decimal import Decimal, InvalidOperation
import json
from pathlib import Path

from .evidence import explicit_percent, sha256
from .package_runtime import load_runtime_package
from .reporting import _check_assignment_scope
from .scope import query_digest

SCHEMA = "1.0.0-dev-genomovar-compatibility"
DEFINITION_RULE = "skani_chromosome_ANI_ge_99.5_complete_linkage"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _ani(value: object) -> Decimal:
    _require(type(value) in (str, int, float), "ANI must be an explicit percent")
    try:
        number = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError("invalid ANI") from exc
    _require(number.is_finite() and 0 <= number <= 100 and number * 100 == (number * 100).to_integral_value(),
             "expected finite ANI percent reported to at most two decimals")
    return number


def assess_groups(groups: dict[str, list[str]], hits: dict[str, dict]) -> dict:
    """Pure all-member99.5% compatibility over explicitly required reference IDs.

    A caller may supply an explicitly reduced group for a recorded development
    diagnostic. This helper never infers self-exclusion from a query ID. The
    package-facing interface always supplies the complete frozen membership.
    """
    _require(isinstance(groups, dict) and bool(groups), "explicit definition groups required")
    _require(all(isinstance(k, str) and k and isinstance(v, list)
                 and all(isinstance(x, str) and x for x in v) for k, v in groups.items()),
             "invalid definition membership")
    members = [r for rs in groups.values() for r in rs]
    _require(len(members) == len(set(members)), "definition members must be unique across groups")
    _require(isinstance(hits, dict) and set(hits) <= set(members), "unexpected definition-member evidence")
    numeric = {}
    for rid, h in hits.items():
        _require(isinstance(h, dict), "member evidence must be an object")
        aq = explicit_percent(h.get("af_query_percent"), "query AF")
        ar = explicit_percent(h.get("af_reference_percent"), "reference AF")
        _require(aq > 0 and ar > 0, "returned ANI requires positive directional AF")
        numeric[rid] = (_ani(h.get("ani_percent")), aq, ar)
    rows = []
    for cid, members in sorted(groups.items()):
        observed = sorted(set(members) & set(numeric))
        missing = sorted(set(members) - set(numeric))
        failing = [r for r in observed if numeric[r][0] < Decimal("99.5")]
        state = ("unassessable_no_other_member" if not members else
                 "incompatible_observed_ANI" if failing else
                 "incomplete_missing_pairs" if missing else "compatible_ANI_only")
        rows.append(dict(concept_id=cid, required_member_count=len(members), observed_member_count=len(observed),
            missing_reference_ids=missing, below_threshold_reference_ids=failing, compatibility_state=state,
            minimum_ani_percent=float(min(numeric[r][0] for r in observed)) if observed else None,
            minimum_query_af_percent=min((numeric[r][1] for r in observed), default=None),
            minimum_reference_af_percent=min((numeric[r][2] for r in observed), default=None)))
    compatible = [r["concept_id"] for r in rows if r["compatibility_state"] == "compatible_ANI_only"]
    incomplete = [r["concept_id"] for r in rows if r["compatibility_state"] == "incomplete_missing_pairs"]
    unassessable = [r["concept_id"] for r in rows if r["compatibility_state"] == "unassessable_no_other_member"]
    result = ("unresolved_unassessable_definition" if unassessable else
              "unresolved_incomplete_evidence" if incomplete else "multiple_compatible_groups" if len(compatible) > 1
              else "one_compatible_group_ANI_only" if len(compatible) == 1 else "no_compatible_group")
    return dict(groups=rows, compatible_concept_ids=compatible, incomplete_concept_ids=incomplete,
                unassessable_concept_ids=unassessable, compatibility_result=result, ani_threshold_percent=99.5,
                threshold_inclusive=True, AF_gate="not_defined", self_exclusion="none_inferred")


class GenomovarContext:
    """Verify the local package once; retain its historical definition component."""

    def __init__(self, manifest_path: Path, *, expected_sha256: str):
        path = Path(manifest_path).absolute()
        _require(path.is_file() and sha256(path) == expected_sha256, "genomovar package SHA256 mismatch")
        view = load_runtime_package(path, expected_sha256=expected_sha256)
        self.manifest = view.manifest
        self.validation = view.validation
        self.package_sha256 = expected_sha256
        self.refs = {r["search_label"]: r for r in self.manifest["references"] if r["component_id"] == "genomovar_chromosome"}
        self.concepts = {c["concept_id"]: c for c in self.manifest["concepts"] if c["rank"] == "genomovar"}
        self.groups = {k: list(c["definition_member_ids"]) for k, c in self.concepts.items()}
        _require(bool(self.groups) and all(c["historical_definition_rule"] == DEFINITION_RULE
                 and c["definition_sequence_scope"] == "historical_chromosome_record_set"
                 and c["production_membership_policy_enabled"] is False for c in self.concepts.values()),
                 "unsupported genomovar definition or enabled policy")

    def evaluate(self, raw_rows: list[dict], *, upstream_evidence: dict, provenance: dict) -> dict:
        """Inspect a completed whole-query/historical-reference search.

        Upstream evidence must already have passed its reporting/package join.
        This function verifies query/package binding but does not authenticate
        the historical execution or independently interpret supplied QC/scope.
        No positive genomovar policy exists, even with one compatible group.
        """
        _require(isinstance(upstream_evidence, dict) and upstream_evidence.get("input_kind") == "assembly",
                 "assembly upstream evidence required")
        assignment = upstream_evidence.get("assignment", {})
        _check_assignment_scope(assignment, upstream_evidence.get("scope_assessment"))
        query = upstream_evidence.get("query_search_path")
        _require(isinstance(query, str) and Path(query).is_absolute() and Path(query).is_file(), "absolute assembly query required")
        query_hash = query_digest(upstream_evidence)
        _require(sha256(Path(query)) == query_hash, "genomovar query hash mismatch")
        _require(upstream_evidence.get("provenance", {}).get("reference_package_sha256") == self.package_sha256,
                 "upstream reference package differs")
        _require(isinstance(provenance, dict) and provenance.get("input_kind") == "assembly"
                 and provenance.get("engine") == "skani" and provenance.get("engine_version") == "0.3.2"
                 and type(provenance.get("search_exit_code")) is int and provenance["search_exit_code"] == 0
                 and provenance.get("reference_component_id") == "genomovar_chromosome"
                 and provenance.get("reference_package_sha256") == self.package_sha256
                 and provenance.get("query_sha256") == query_hash
                 and isinstance(provenance.get("search_parameters_record"), str) and provenance["search_parameters_record"],
                 "completed pinned genomovar search provenance required")
        _require(isinstance(raw_rows, list), "raw genomovar search rows required")
        hits = {}
        for r in raw_rows:
            _require(isinstance(r, dict) and r.get("Query_file") == query and r.get("Ref_file") in self.refs,
                     "genomovar query/reference path mismatch")
            rid = self.refs[r["Ref_file"]]["reference_id"]
            _require(rid not in hits, "duplicate genomovar-member comparison")
            hits[rid] = dict(ani_percent=r.get("ANI"), af_query_percent=r.get("Align_fraction_query"),
                             af_reference_percent=r.get("Align_fraction_ref"))
        result = assess_groups(self.groups, hits)
        for row in result["groups"]:
            c = self.concepts[row["concept_id"]]
            row.update(display_name=c["display_name"], definition_version=c["definition_version"],
                       parent_concept_id=c["parent_concept_id"], sampling_state=c["sampling_state"],
                       definition_rule=c["historical_definition_rule"])
        result.update(schema_version=SCHEMA, record_type="development_genomovar_compatibility_only",
            query_search_path=query, query_sha256=query_hash, reference_package_sha256=self.package_sha256,
            definition_reference_scope="historical_chromosome_record_set", query_scope="whole_assembly_as_supplied",
            reference_count=len(self.refs), returned_reference_count=len(hits),
            provenance=copy.deepcopy(provenance), input_assignment=copy.deepcopy(assignment),
            assigned_genomovar_id=None, membership_policy_enabled=False,
            upstream_eligibility="separate_not_established_here",
            supplied_scope_outcome=upstream_evidence.get("scope_assessment", {}).get("outcome", "not_assessed"),
            supplied_qc_status=upstream_evidence.get("qc", {}).get("status", "not_assessed"),
            limitations=["Compatibility is not an assignment or independent validation.",
                "Every frozen member is required; no self-exclusion is inferred from an accession string.",
                "No coverage gate is adopted here; ANI compatibility does not establish adequate coverage or purity.",
                "Missing alternative-group evidence prevents a unique compatibility conclusion.",
                "A singleton definition cannot establish its population boundary from one reference."])
        return result
