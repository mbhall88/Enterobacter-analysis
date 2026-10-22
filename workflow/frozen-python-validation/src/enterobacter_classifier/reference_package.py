"""Pinned, local development reference packages; no classification policy.

The recipe identifies reviewed tables and pre-existing skani databases by hash.
Building and validating both check all packaged FASTAs, database files and index
membership. This is an integrity/contract check, not authenticated provenance or
an independent recomputation of sketches, taxonomy or biological scope.
"""
from __future__ import annotations

import argparse
import copy
import json
import re
import struct
from collections import Counter
from pathlib import Path

from .evidence import read_tsv, sha256

SCHEMA = "1.0.0-dev-reference-package"
INDEX_FORMAT = "skani-0.3.2-consolidated-bincode1"
COMPONENTS = {
    "primary_whole": ("primary_named_reference_evidence", "whole_assembly_as_deposited"),
    "genomovar_chromosome": ("genomovar_definition_member", "historical_chromosome_record_set"),
}
TABLE_KEYS = {
    "named_registry": "concept_id", "lineage_registry": "concept_id",
    "primary_references": "accession", "definition_references": "accession",
    "definition_memberships": "accession", "definition_groups": "genomovar_id",
    "genomovar_crosswalk": "genomovar_id", "development_context": "accession",
    "cohort_coverage": "accession", "historical_types": "historical_taxon_name",
    "assignments": "accession",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def canonical_json(value: dict) -> str:
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n"


def file_record(path: str | Path) -> dict:
    """Create an explicit pin; callers must supply the intended source path."""
    p = Path(path).absolute()
    require(p.is_file(), f"missing source file: {p}")
    return {"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size}


def unique(rows: list[dict], key: str, label: str) -> dict[str, dict]:
    values = [r.get(key, "") for r in rows]
    require(all(isinstance(x, str) and x and x == x.strip() for x in values),
            f"{label}: empty/invalid {key}")
    require(len(set(values)) == len(values), f"{label}: duplicate {key}")
    return dict(zip(values, rows))


def read_skani_index(path: Path, sketches_bytes: int) -> list[dict]:
    """Read only the pinned v0.3.2 format (Vec<String,u64,u64>, bincode1).

    The source writer appends every sketch contiguously; all offsets/lengths and
    trailing bytes are checked. Sketch payloads and markers are pinned by hash,
    not decoded or biologically revalidated here.
    """
    data, cursor = path.read_bytes(), 0

    def number() -> int:
        nonlocal cursor
        require(cursor + 8 <= len(data), "truncated skani index")
        value = struct.unpack_from("<Q", data, cursor)[0]
        cursor += 8
        return value

    count = number()
    require(0 < count <= len(data) // 24, "invalid skani index count")
    entries, expected_offset = [], 0
    for _ in range(count):
        length = number()
        require(0 < length <= len(data) - cursor, "invalid skani index filename length")
        try:
            name = data[cursor:cursor + length].decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError("invalid skani index filename encoding") from exc
        cursor += length
        offset, size = number(), number()
        require(Path(name).is_absolute() and "\x00" not in name, "index path must be absolute")
        require(offset == expected_offset and size > 0, "noncontiguous/empty skani index segment")
        expected_offset = offset + size
        require(expected_offset <= sketches_bytes, "skani index exceeds sketches file")
        entries.append({"path": name, "offset": offset, "bytes": size})
    require(cursor == len(data), "trailing skani index data")
    require(expected_offset == sketches_bytes, "unindexed bytes in sketches file")
    unique(entries, "path", "skani index")
    return entries


class _Builder:
    def __init__(self, recipe: dict):
        self.recipe = copy.deepcopy(recipe)
        self.checked: dict[str, dict] = {}

    def artifact(self, item: dict) -> dict:
        require(isinstance(item, dict) and set(item) == {"path", "sha256", "bytes"},
                "artifact requires exactly path, sha256, bytes")
        path, digest, size = item["path"], item["sha256"], item["bytes"]
        require(isinstance(path, str) and Path(path).is_absolute(), "artifact path must be absolute")
        require(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
                "artifact requires lowercase SHA-256")
        require(type(size) is int and size >= 0, "artifact bytes must be a nonnegative integer")
        if path in self.checked:
            require(self.checked[path] == item, "conflicting pins for one artifact path")
        else:
            require(Path(path).is_file(), f"missing artifact: {path}")
            require(Path(path).stat().st_size == size, f"artifact size mismatch: {path}")
            require(sha256(Path(path)) == digest, f"artifact SHA-256 mismatch: {path}")
            self.checked[path] = dict(item)
        return dict(item)

    def sequence(self, path: str, digest: str) -> dict:
        require(Path(path).is_absolute() and Path(path).is_file(), f"missing sequence source: {path}")
        return self.artifact({"path": path, "sha256": digest, "bytes": Path(path).stat().st_size})

    def build(self) -> dict:
        recipe = self.recipe
        require(set(recipe) == {"package_id", "package_version", "source_tables", "components",
                                "provenance_files", "expected_counts"}, "unknown/missing recipe fields")
        for key in ("package_id", "package_version"):
            require(isinstance(recipe[key], str) and re.fullmatch(r"[A-Za-z0-9._-]+", recipe[key]) is not None,
                    f"invalid {key}")
        require(set(recipe["source_tables"]) == set(TABLE_KEYS), "missing/unknown source tables")
        require(set(recipe["components"]) == set(COMPONENTS), "missing/unknown reference components")
        tables = {}
        for name, key in TABLE_KEYS.items():
            source = self.artifact(recipe["source_tables"][name])
            rows = read_tsv(Path(source["path"]), {key})
            tables[name] = unique(rows, key, name)
        require(bool(recipe["provenance_files"]), "source provenance records required")
        for artifact in recipe["provenance_files"]:
            self.artifact(artifact)

        concepts = []
        for name in ("named_registry", "lineage_registry"):
            for row in tables[name].values():
                rank = row["rank"]
                require(rank in ({"species", "subspecies"} if name == "named_registry" else {"lineage"}),
                        "undefined or misplaced concept rank")
                require(row["validated_supported_caller_rank"] == "False", "positive caller claim in proposal")
                n = int(row["primary140_reference_count"])
                require(n >= 0 and row["represented_reference_evidence"] == str(n > 0),
                        "inconsistent reference representation state")
                if name == "lineage_registry":
                    require(n == 0 and row["genomovar173_reference_count"] == "0", "context lineage promoted")
                concepts.append({
                    "concept_id": row["concept_id"], "rank": rank,
                    "parent_concept_id": row.get("parent_concept_id", "") or None,
                    "display_name": row.get("project_name", row.get("project_label")),
                    "identifier_status": "stable_id_proposal_not_adopted_registry",
                    "nomenclature_status": row["nomenclature_status"],
                    "reference_state": "represented" if n else "context_only" if rank == "lineage" else "unsupported_no_reference",
                    "sampling_state": "absent" if not n else "single_reference" if n == 1 else "multiple_development_references",
                    "validated_supported_caller_rank": False,
                    "source_annotation": dict(row),
                })
        by_id = unique(concepts, "concept_id", "concept registry")
        by_name = unique(concepts, "display_name", "concept display names")
        for c in concepts:
            parent = c["parent_concept_id"]
            if c["rank"] == "subspecies":
                require(parent in by_id and by_id[parent]["rank"] == "species", "undefined/non-species parent")
            else:
                require(parent is None, "species/context lineage cannot have an inferred parent")

        primary = tables["primary_references"]
        definitions = tables["definition_references"]
        memberships = tables["definition_memberships"]
        groups = tables["definition_groups"]
        crosswalk = tables["genomovar_crosswalk"]
        require(set(definitions) == set(memberships), "silent missing/extra definition member")
        require(set(groups) == set(crosswalk) == {r["genomovar_id"] for r in definitions.values()},
                "missing/extra genomovar definition group")
        refs = []
        for component, rows in (("primary_whole", primary), ("genomovar_chromosome", definitions)):
            role, scope = COMPONENTS[component]
            for accession, row in rows.items():
                require(row["partition"] == "development", "non-development reference in package")
                require(row["group"] in by_name, "undefined reference concept")
                c = by_name[row["group"]]
                require(c["rank"] in {"species", "subspecies"}, "context reference promoted to named evidence")
                species = c if c["rank"] == "species" else by_id[c["parent_concept_id"]]
                overlay = tables["assignments"].get(accession)
                require(overlay is not None and overlay["project_species"] == species["display_name"],
                        "reference species disagrees with reviewed overlay")
                subspecies = c["display_name"] if c["rank"] == "subspecies" else ""
                expected_subspecies = f'{species["display_name"]} subsp. {overlay["project_subspecies"]}' if overlay["project_subspecies"] else ""
                require(subspecies == expected_subspecies, "reference subspecies disagrees with reviewed overlay")
                require("no_named_classifier_reference" not in overlay["usage_restrictions"], "restricted context reference promoted")
                if component == "primary_whole":
                    require(row["project_species"] == species["display_name"] and row["named_target"] == "1",
                            "primary reference is not its stated named target")
                    require(row["project_subspecies"] == overlay["project_subspecies"], "primary subspecies join mismatch")
                    source = self.sequence(row["whole_fasta"], row["whole_sha256"])
                    is_type = row["type_anchor"]
                else:
                    require(c["rank"] == "subspecies", "genomovar definition requires explicit subspecies parent")
                    source = self.sequence(row["chromosome_fasta"], row["chromosome_sha256"])
                    is_type = row["is_type"]
                    m = memberships[accession]
                    require(m["project_species"] == species["display_name"] and m["project_subspecies_proposal"] == subspecies,
                            "definition membership taxonomy mismatch")
                    for manifest_field, member_field in (("genomovar_id", "genomovar_id"),
                                                          ("definition_version", "genomovar_version"),
                                                          ("definition_rule", "genomovar_rule")):
                        require(row[manifest_field] == m[member_field] == overlay[member_field],
                                f"definition/overlay {manifest_field} mismatch")
                require(is_type in {"0", "1"} and row["qc_pass"] in {"0", "1"}, "invalid reference QC/type flag")
                refs.append({
                    "reference_id": f"{component}:{accession}", "accession": accession,
                    "component_id": component, "reference_role": role, "sequence_scope": scope,
                    "species_id": species["concept_id"], "subspecies_id": c["concept_id"] if subspecies else None,
                    "stratum_concept_id": c["concept_id"], "source": source,
                    "search_label": source["path"],
                    "type_anchor": is_type == "1", "supplementary_qc_pass": row["qc_pass"] == "1",
                    "quality_flags": row["quality_flags"],
                    "quality_state": "type_QC_exception" if is_type == "1" and row["qc_pass"] == "0" else "source_QC_pass" if row["qc_pass"] == "1" else "source_QC_fail",
                    "source_selection_reason": row.get("selection_reason", "frozen_genomovar_definition_member"),
                    "historical_origin_path": row.get("chromosome_original_source"),
                    "historical_origin_validation": "frozen_metadata_only_not_an_additional_packaged_source",
                    "historical_source_basis": row.get("chromosome_source_basis"),
                    "assignment_status": overlay["assignment_status"],
                })
        unique(refs, "reference_id", "reference registry")
        for c in concepts:
            if c["rank"] == "lineage":
                continue
            represented = [r for r in refs if r["component_id"] == "primary_whole" and
                           c["concept_id"] == r[c["rank"] + "_id"]]
            row = c["source_annotation"]
            require(len(represented) == int(row["primary140_reference_count"]), "concept reference count mismatch")
            require(sum(r["type_anchor"] for r in represented) == int(row["primary140_type_count"]), "type-anchor count mismatch")
            for field, predicate in (
                ("all_represented_historical_type_accessions", lambda r: r["type_anchor"]),
                ("supplementary_reference_accessions", lambda r: not r["type_anchor"]),
                ("type_QC_exception_accessions", lambda r: r["quality_state"] == "type_QC_exception"),
            ):
                expected = set(filter(None, row[field].split(";")))
                require({r["accession"] for r in represented if predicate(r)} == expected, f"{field} membership mismatch")

        definition_groups = []
        for gid, g in groups.items():
            rows = [r for r in definitions.values() if r["genomovar_id"] == gid]
            x = crosswalk[gid]
            parents = {by_name[r["group"]]["concept_id"] for r in rows}
            versions = {r["definition_version"] for r in rows}
            rules = {r["definition_rule"] for r in rows}
            require(parents == {x["parent_concept_id"]} and by_id[x["parent_concept_id"]]["rank"] == "subspecies",
                    "undefined/wrong genomovar parent")
            require(versions == {x["version"]} and rules == {x["historical_rule"]}, "mixed definition version/rule")
            require(len(rows) == int(g["genomes"]) == int(x["definition_member_count"]), "silent missing definition group member")
            accessions = {r["accession"] for r in rows}
            require(g["representative_accession"] == x["representative_accession"] and g["representative_accession"] in accessions,
                    "representative not a definition member")
            cid = f'ECTX-GV-{x["version"]}-{gid}'
            definition_groups.append({
                "concept_id": cid, "rank": "genomovar", "parent_concept_id": x["parent_concept_id"],
                "display_name": gid, "identifier_status": "stable_id_proposal_not_adopted_registry",
                "nomenclature_status": "informal_operational_ANI_group", "definition_version": x["version"],
                "historical_definition_rule": x["historical_rule"],
                "historical_definition_engine": "skani 0.3.1", "definition_sequence_scope": COMPONENTS["genomovar_chromosome"][1],
                "definition_member_ids": sorted(f"genomovar_chromosome:{a}" for a in accessions),
                "representative_reference_id": f'genomovar_chromosome:{g["representative_accession"]}',
                "representative_role": "nomination_only_not_membership_shortcut",
                "sampling_state": "singleton_no_independent_own_group_test" if len(rows) == 1 else "multiple_definition_members",
                "reference_state": "represented_definition", "validated_supported_caller_rank": False,
                "production_membership_policy_enabled": False,
            })
        unique(concepts + definition_groups, "concept_id", "full concept registry")

        primary_by_accession = {r["accession"]: r for r in refs if r["component_id"] == "primary_whole"}
        type_accessions = set()
        for row in tables["historical_types"].values():
            represented = set(filter(None, row["represented_name_bearing_accessions"].split(";")))
            for accession in represented:
                require(accession in primary_by_accession and primary_by_accession[accession]["type_anchor"],
                        "historical type crosswalk points outside primary type anchors")
                ref = primary_by_accession[accession]
                require(row["current_project_species"] == by_id[ref["species_id"]]["display_name"] and
                        row["current_project_stratum"] == by_id[ref["stratum_concept_id"]]["display_name"],
                        "historical type crosswalk concept mismatch")
            type_accessions.update(represented)
        require(type_accessions == {a for a, r in primary_by_accession.items() if r["type_anchor"]},
                "historical type crosswalk silently omits a type anchor")

        context = tables["development_context"]
        require(not (set(context) & (set(primary) | set(definitions))), "context-only assembly in packaged references")
        for row in context.values():
            require(row["reference_in_proposed_package"] == "False", "context-only package claim changed")
        coverage = tables["cohort_coverage"]
        for a, row in coverage.items():
            require(row["validated_automatic_assignment"] == "False", "coverage audit implies validated caller")
            require(row["in_primary140"] == str(a in primary) and row["in_genomovar173"] == str(a in definitions),
                    "cohort/reference membership mismatch")
            for rank in ("species", "subspecies"):
                cid = row[rank + "_concept_id"]
                require(not cid or (cid in by_id and by_id[cid]["rank"] == rank), "undefined coverage concept/rank")
                present = bool(cid) and by_id[cid]["reference_state"] == "represented"
                require(row[rank + "_reference_evidence_present"] == str(present), "cohort reference coverage mismatch")

        components = []
        for name, (role, scope) in COMPONENTS.items():
            configuration = recipe["components"][name]
            require(set(configuration) == {"source_list", "database_files", "build_records", "skani_binary"}, "component fields missing/unknown")
            source_list = self.artifact(configuration["source_list"])
            paths = Path(source_list["path"]).read_text().splitlines()
            declared = [r["search_label"] for r in refs if r["component_id"] == name]
            require(len(paths) == len(set(paths)) == len(declared) == len(set(declared)) and set(paths) == set(declared),
                    "duplicate/missing/mixed-scope source-list member")
            files = [self.artifact(a) for a in configuration["database_files"]]
            file_by_name = unique([{**f, "name": Path(f["path"]).name} for f in files], "name", "database inventory")
            require(set(file_by_name) == {"index.db", "markers.bin", "sketches.db"}, "database inventory incomplete/unknown")
            directories = {Path(f["path"]).parent for f in files}
            require(len(directories) == 1, "mixed database directories")
            directory = directories.pop()
            require({p.name for p in directory.iterdir()} == set(file_by_name), "database directory inventory mismatch")
            index = read_skani_index(Path(file_by_name["index.db"]["path"]), file_by_name["sketches.db"]["bytes"])
            require({e["path"] for e in index} == set(declared) and len(index) == len(declared), "database index reference/scope mismatch")
            require(bool(configuration["build_records"]), "archived build records required")
            for a in configuration["build_records"]:
                self.artifact(a)
            self.artifact(configuration["skani_binary"])
            components.append({
                "component_id": name, "reference_role": role, "sequence_scope": scope,
                "sequence_scope_basis": "frozen_source_manifest_not_new_replicon_detection",
                "reference_ids": sorted(r["reference_id"] for r in refs if r["component_id"] == name),
                "engine": "skani", "engine_version": "0.3.2", "database_format": INDEX_FORMAT,
                "database_files": sorted(files, key=lambda f: f["path"]),
                "database_index_entries": sorted(index, key=lambda e: e["path"]),
                "search_label_semantics": "exact_archived_index_label_not_a_deployed_FASTA_requirement",
                "provenance_state": "archived_build_records_hash_checked_not_reexecuted",
                "sketch_payload_recomputed": False,
            })

        counts = {
            "primary_whole_references": len(primary), "historical_chromosome_definition_references": len(definitions),
            "unique_accessions": len(set(primary) | set(definitions)), "reference_scope_records": len(refs),
            "named_concepts": len(tables["named_registry"]), "context_lineage_concepts": len(tables["lineage_registry"]),
            "genomovar_concepts": len(groups), "development_context_records_not_packaged": len(context),
            "eligible_cohort_records": len(coverage),
            "eligible_species_reference_represented": sum(r["species_reference_evidence_present"] == "True" for r in coverage.values()),
            "primary_historical_type_accessions": sum(r["type_anchor"] for r in refs if r["component_id"] == "primary_whole"),
        }
        require(recipe["expected_counts"] == counts and all(type(v) is int for v in recipe["expected_counts"].values()),
                "frozen package count contract mismatch")
        return {
            "schema_version": SCHEMA, "package_id": recipe["package_id"], "package_version": recipe["package_version"],
            "release_state": "development_proposal_not_distributable", "positive_calling_policy_enabled": False,
            "path_semantics": "absolute_existing_local_sources_no_relocation_or_sequence_copy",
            "integrity_scope": "all_pins_and_packaged_sequence_hashes_exact_source_joins_and_database_index_membership",
            "provenance_limit": "hashes_pin_reviewed_local_records_not_authenticated_origin_or_recomputed_sketches",
            "recipe": recipe, "counts": counts,
            "concepts": sorted(concepts + definition_groups, key=lambda c: c["concept_id"]),
            "references": sorted(refs, key=lambda r: r["reference_id"]),
            "components": sorted(components, key=lambda c: c["component_id"]),
            "historical_type_crosswalk": sorted(tables["historical_types"].values(), key=lambda r: r["historical_taxon_name"]),
            "context_not_packaged": [{k: v for k, v in row.items() if k not in {"whole_fasta", "whole_sha256"}}
                                     for _, row in sorted(context.items())],
            "cohort_coverage_summary": dict(sorted(Counter(r["coverage_category"] for r in coverage.values()).items())),
            "verified_artifacts": sorted(self.checked.values(), key=lambda a: a["path"]),
        }


def build_package(recipe: dict) -> dict:
    """Reconstruct only the frozen reviewed proposal, verifying external files."""
    try:
        return _Builder(recipe).build()
    except (KeyError, TypeError, AttributeError) as exc:
        raise ValueError(f"malformed package recipe/source table: {exc}") from exc


def validate_package(manifest: dict) -> dict:
    """Rebuild from pinned sources and require exact, type-sensitive equality."""
    require(isinstance(manifest, dict) and isinstance(manifest.get("recipe"), dict), "manifest recipe missing")
    expected = build_package(manifest["recipe"])
    require(canonical_json(manifest) == canonical_json(expected), "manifest differs from checked frozen source joins")
    return {"valid": True, "schema_version": SCHEMA, "counts": expected["counts"],
            "verified_artifacts": len(expected["verified_artifacts"]),
            "database_bytes": sum(f["bytes"] for c in expected["components"] for f in c["database_files"]),
            "positive_calling_policy_enabled": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    build = sub.add_parser("build")
    build.add_argument("--recipe", required=True, type=Path)
    build.add_argument("--output", required=True, type=Path)
    validate = sub.add_parser("validate")
    validate.add_argument("--manifest", required=True, type=Path)
    validate.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    output = args.output if args.command == "build" else args.report
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    try:
        value = (build_package(json.loads(args.recipe.read_text())) if args.command == "build"
                 else validate_package(json.loads(args.manifest.read_text())))
        with output.open("x") as handle:
            handle.write(canonical_json(value))
    except (ValueError, OSError) as exc:
        parser.exit(2, f"reference package error: {exc}\n")


if __name__ == "__main__":
    main()
