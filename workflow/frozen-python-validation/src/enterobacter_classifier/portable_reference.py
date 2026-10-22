"""Relocate frozen reference sketches and metadata without claiming a release.

Export verifies the full local derivation contract. Loading instead trusts an
externally supplied package digest, checks copied artifacts and index mappings,
and never dereferences archived FASTA or executable paths. Database bytes and
opaque index labels are unchanged; no engine or biological analysis runs here.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path, PurePosixPath
import re
import shutil

from .evidence import read_tsv, sha256
from .reference_package import (COMPONENTS, INDEX_FORMAT, SCHEMA as LOCAL_SCHEMA,
                                TABLE_KEYS, canonical_json, file_record,
                                read_skani_index, require, unique, validate_package)

SCHEMA = "1.0.0-dev-portable-reference-package"
TRUST = "expected_package_digest_and_export_attestation_not_authenticated_release"
RELEASE = "development_proposal_not_distributable"


def _record(record: dict) -> None:
    require(isinstance(record, dict) and set(record) == {"path", "sha256", "bytes"}, "invalid artifact record")
    require(isinstance(record["path"], str) and bool(record["path"]), "artifact path required")
    require(isinstance(record["sha256"], str) and re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is not None,
            "explicit lowercase SHA256 required")
    require(type(record["bytes"]) is int and record["bytes"] >= 0, "invalid artifact byte count")


def _contained(root: Path, name: str) -> Path:
    require(isinstance(name, str) and bool(name) and "\\" not in name and "\x00" not in name,
            "unsafe portable relative path")
    parts = PurePosixPath(name)
    require(not parts.is_absolute() and parts.as_posix() == name
            and all(p not in {".", ".."} for p in parts.parts), "unsafe portable relative path")
    path = root
    for part in parts.parts:
        path = path / part
        require(not path.is_symlink(), "portable assets cannot be symlinks")
    require(path.resolve().is_relative_to(root.resolve()), "portable path escapes package")
    return path


def _check_file(path: Path, record: dict) -> None:
    _record(record)
    require(path.is_file() and path.stat().st_size == record["bytes"]
            and sha256(path) == record["sha256"], "portable artifact missing or hash/size mismatch")


def _asset_specs(source: dict) -> list[dict]:
    """Exact recipe allowlist; no directory globs, FASTAs or binary copying."""
    recipe, records = source["recipe"], {}
    require(set(recipe["source_tables"]) == set(TABLE_KEYS), "source table contract mismatch")
    require(set(recipe["components"]) == set(COMPONENTS), "source component contract mismatch")
    def add(record, target, role):
        _record(record)
        require(Path(record["path"]).is_absolute(), "archived source path must be absolute")
        existing = records.get(record["path"])
        if existing is not None:
            require(existing["sha256"] == record["sha256"] and existing["bytes"] == record["bytes"], "conflicting source pins")
            existing["roles"].append(role)
        else:
            records[record["path"]] = dict(path=target, source_path=record["path"],
                sha256=record["sha256"], bytes=record["bytes"], roles=[role])
    for name in sorted(TABLE_KEYS):
        add(recipe["source_tables"][name], f"metadata/tables/{name}.tsv", f"table:{name}")
    for i, record in enumerate(recipe["provenance_files"]):
        add(record, f"metadata/provenance/{i:03d}.artifact", f"provenance:{i}")
    for cid in sorted(COMPONENTS):
        c = recipe["components"][cid]
        add(c["source_list"], f"metadata/lists/{cid}.txt", f"source_list:{cid}")
        for i, record in enumerate(c["build_records"]):
            add(record, f"metadata/build/{cid}/{i:03d}.artifact", f"build:{cid}:{i}")
        names = [Path(f["path"]).name for f in c["database_files"]]
        require(len(names) == 3 and set(names) == {"index.db", "sketches.db", "markers.bin"}, "database file contract mismatch")
        for record in c["database_files"]:
            name = Path(record["path"]).name
            add(record, f"databases/{cid}/{name}", f"database:{cid}:{name}")
    excluded = {r["source"]["path"] for r in source["references"]}
    excluded.update(c["skani_binary"]["path"] for c in recipe["components"].values())
    require(not excluded.intersection(records), "sequence or executable leaked into metadata allowlist")
    require(len({r["path"] for r in records.values()}) == len(records), "duplicate portable asset path")
    for row in records.values(): row["roles"].sort()
    return sorted(records.values(), key=lambda r: r["path"])


def _model(source: dict) -> None:
    """Check frozen metadata invariants without reopening original sequences."""
    require(source.get("schema_version") == LOCAL_SCHEMA and source.get("release_state") == RELEASE
            and source.get("positive_calling_policy_enabled") is False, "unsupported or enabled source package")
    refs = unique(source["references"], "reference_id", "references")
    concepts = unique(source["concepts"], "concept_id", "concepts")
    components = unique(source["components"], "component_id", "components")
    require(set(components) == set(COMPONENTS), "wrong component inventory")
    require(source["counts"] == source["recipe"]["expected_counts"], "frozen count contract mismatch")
    for c in concepts.values():
        require(c["rank"] in {"species", "subspecies", "lineage", "genomovar"}
                and c["validated_supported_caller_rank"] is False, "invalid concept rank/caller state")
        parent = concepts.get(c.get("parent_concept_id"))
        if c["rank"] in {"subspecies", "genomovar"}:
            require(parent is not None and parent["rank"] == ("species" if c["rank"] == "subspecies" else "subspecies"),
                    "undefined/wrong concept parent")
        else:
            require(c.get("parent_concept_id") is None, "unexpected concept parent")
        if c["rank"] == "genomovar":
            require(c["production_membership_policy_enabled"] is False
                    and c["definition_sequence_scope"] == COMPONENTS["genomovar_chromosome"][1], "wrong genomovar scope/policy")
            members = c["definition_member_ids"]
            require(members and len(set(members)) == len(members) and set(members) <= set(refs), "invalid definition membership")
            require(c["representative_reference_id"] in members, "representative outside definition")
            require(all(refs[r]["component_id"] == "genomovar_chromosome"
                        and refs[r]["subspecies_id"] == c["parent_concept_id"] for r in members), "definition member scope/parent mismatch")
    memberships = [r for c in concepts.values() if c["rank"] == "genomovar" for r in c["definition_member_ids"]]
    definition_ids = {r for r, v in refs.items() if v["component_id"] == "genomovar_chromosome"}
    require(len(memberships) == len(set(memberships)) and set(memberships) == definition_ids, "definition members omitted/duplicated")
    for rid, r in refs.items():
        cid = r["component_id"]
        require(cid in COMPONENTS and rid == f"{cid}:{r['accession']}", "unstable/wrong-scope reference ID")
        role, scope = COMPONENTS[cid]
        require((r["reference_role"], r["sequence_scope"]) == (role, scope), "reference scope/role mismatch")
        _record(r["source"])
        require(Path(r["source"]["path"]).is_absolute() and r["search_label"] == r["source"]["path"], "archived index/source label mismatch")
        require(concepts.get(r["species_id"], {}).get("rank") == "species", "unknown reference species")
        sub = r["subspecies_id"]
        require(sub is None or (concepts.get(sub, {}).get("rank") == "subspecies"
                and concepts[sub]["parent_concept_id"] == r["species_id"]), "wrong reference subspecies parent")
        require(r["stratum_concept_id"] == (sub or r["species_id"]), "wrong reference stratum")
        require(type(r["type_anchor"]) is bool and type(r["supplementary_qc_pass"]) is bool, "invalid type/QC field types")
        quality = "type_QC_exception" if r["type_anchor"] and not r["supplementary_qc_pass"] else "source_QC_pass" if r["supplementary_qc_pass"] else "source_QC_fail"
        require(r["quality_state"] == quality, "reference QC annotation mismatch")
    for cid, c in components.items():
        require((c["reference_role"], c["sequence_scope"]) == COMPONENTS[cid], "component scope/role mismatch")
        require(c["engine"] == "skani" and c["engine_version"] == "0.3.2" and c["database_format"] == INDEX_FORMAT
                and c["sketch_payload_recomputed"] is False, "unsupported database engine/derivation")
        members = {k for k, r in refs.items() if r["component_id"] == cid}
        require(len(c["reference_ids"]) == len(members) and set(c["reference_ids"]) == members, "component reference membership mismatch")
        unique([refs[k] for k in members], "search_label", "component index labels")
        require(canonical_json({"files": c["database_files"]}) == canonical_json({"files": sorted(source["recipe"]["components"][cid]["database_files"], key=lambda x:x["path"])}), "component/recipe database mismatch")
    counts = source["counts"]
    expected = dict(primary_whole_references=sum(r["component_id"] == "primary_whole" for r in refs.values()),
        historical_chromosome_definition_references=len(definition_ids), unique_accessions=len({r["accession"] for r in refs.values()}),
        reference_scope_records=len(refs), named_concepts=sum(c["rank"] in {"species", "subspecies"} for c in concepts.values()),
        context_lineage_concepts=sum(c["rank"] == "lineage" for c in concepts.values()),
        genomovar_concepts=sum(c["rank"] == "genomovar" for c in concepts.values()),
        development_context_records_not_packaged=len(source["context_not_packaged"]),
        primary_historical_type_accessions=sum(r["component_id"] == "primary_whole" and r["type_anchor"] for r in refs.values()))
    require(all(type(counts.get(k)) is int and counts[k] == v for k,v in expected.items()), "reference/concept count mismatch")
    require(not ({r["accession"] for r in source["context_not_packaged"]} & {r["accession"] for r in refs.values()}), "context promoted into references")


class PortableReference:
    """A checked runtime view; archived source paths remain provenance only."""
    def __init__(self, path: Path, digest: str, runtime: dict, validation: dict):
        self.package_path = str(path)
        self.package_sha256 = digest
        self._manifest = runtime
        self.validation = validation

    @property
    def manifest(self) -> dict:
        return copy.deepcopy(self._manifest)

    def engine_record(self, component_id: str, binary_path: Path) -> dict:
        require(component_id in COMPONENTS, "unknown reference component")
        require(binary_path is not None, "portable reference requires an explicit skani binary")
        original = self._manifest["recipe"]["components"][component_id]["skani_binary"]
        actual = Path(binary_path).absolute()
        _check_file(actual, original)
        return {**original, "path": str(actual)}


def load_portable(manifest_path: Path, *, expected_sha256: str) -> PortableReference:
    """Load only copied artifacts; the expected digest supplies release trust.

    Full source-FASTA derivation is attested by export, not rerun at deployment.
    This API neither authenticates a publisher nor enables reference release.
    """
    try:
        path = Path(manifest_path).absolute()
        require(path.name == "portable.json" and not path.is_symlink(), "expected regular portable.json")
        require(re.fullmatch(r"[0-9a-f]{64}", expected_sha256 or "") is not None
                and path.is_file() and sha256(path) == expected_sha256, "portable manifest SHA256 mismatch")
        root = path.parent.resolve()
        p = json.loads(path.read_text())
        require(set(p) == {"schema_version", "package_id", "package_version", "release_state", "trust_basis",
                          "positive_calling_policy_enabled", "source_manifest", "export_validation", "assets"}, "unknown/missing portable fields")
        require(p["schema_version"] == SCHEMA and p["release_state"] == RELEASE and p["trust_basis"] == TRUST
                and p["positive_calling_policy_enabled"] is False, "unsupported portable policy/schema")
        for key in ("source_manifest", "export_validation"):
            _check_file(_contained(root, p[key]["path"]), p[key])
        original = json.loads(_contained(root, p["source_manifest"]["path"]).read_text())
        _model(original)
        require((p["package_id"], p["package_version"]) == (original["package_id"], original["package_version"]), "package identity/version changed")
        expected = _asset_specs(original)
        require(canonical_json({"assets": p["assets"]}) == canonical_json({"assets": expected}), "portable inventory differs from recipe allowlist")
        copied = {}
        for a in expected:
            file = _contained(root, a["path"])
            _check_file(file, {k:a[k] for k in ("path", "sha256", "bytes")})
            copied[a["source_path"]] = str(file)
        allowed = {a["path"] for a in expected} | {p["source_manifest"]["path"], p["export_validation"]["path"], "portable.json"}
        require(len(allowed) == len(expected) + 3, "colliding portable paths")
        observed = set()
        for item in root.rglob("*"):
            require(not item.is_symlink(), "portable package contains a symlink")
            if item.is_file(): observed.add(item.relative_to(root).as_posix())
        require(observed == allowed, "portable directory inventory mismatch")
        attestation = json.loads(_contained(root, p["export_validation"]["path"]).read_text())
        require(attestation["source_manifest_sha256"] == p["source_manifest"]["sha256"]
                and attestation["local_validation"]["valid"] is True
                and attestation["local_validation"]["counts"] == original["counts"]
                and attestation["scope"] == "full_local_source_validation_at_export_not_reexecuted_at_load", "invalid export derivation attestation")
        runtime = copy.deepcopy(original)
        for table in runtime["recipe"]["source_tables"].values(): table["path"] = copied[table["path"]]
        for record in runtime["recipe"]["provenance_files"]: record["path"] = copied[record["path"]]
        index_count = 0
        for c in runtime["components"]:
            cid = c["component_id"]
            files = {Path(f["path"]).name: f for f in c["database_files"]}
            index = read_skani_index(Path(copied[files["index.db"]["path"]]), files["sketches.db"]["bytes"])
            require(sorted(index, key=lambda x:x["path"]) == c["database_index_entries"], "copied/archived database index mismatch")
            labels = {r["search_label"] for r in runtime["references"] if r["component_id"] == cid}
            require(len(index) == len(labels) and {e["path"] for e in index} == labels, "database index reference/scope mismatch")
            configuration = runtime["recipe"]["components"][cid]
            source_list = Path(copied[configuration["source_list"]["path"]]).read_text().splitlines()
            require(len(source_list) == len(labels) and set(source_list) == labels, "source-list index membership mismatch")
            for record in c["database_files"]: record["path"] = copied[record["path"]]
            for record in [configuration["source_list"], *configuration["database_files"], *configuration["build_records"]]:
                record["path"] = copied[record["path"]]
            index_count += len(index)
        primary = runtime["recipe"]["source_tables"]["primary_references"]
        rows = read_tsv(Path(primary["path"]), {"accession"})
        require(set(unique(rows, "accession", "primary source rows")) == {r["accession"] for r in runtime["references"] if r["component_id"] == "primary_whole"}, "primary source/reference membership mismatch")
        runtime["path_semantics"] = "relocated_assets_with_archived_labels_and_sequence_paths_as_opaque_provenance"
        runtime["integrity_scope"] = "copied_artifacts_index_mapping_and_frozen_metadata_invariants;source_derivation_attested_at_export"
        validation = dict(valid=True, schema_version=SCHEMA, counts=original["counts"],
            copied_artifacts=len(expected)+2, database_index_entries=index_count,
            database_bytes=sum(f["bytes"] for c in runtime["components"] for f in c["database_files"]),
            original_fastas_required_at_load=False, sketched_sequences_recomputed=False,
            trust_basis=TRUST, source_derivation="full_local_validation_attested_at_export_only",
            positive_calling_policy_enabled=False, source_manifest_sha256=p["source_manifest"]["sha256"])
        require(sha256(path) == expected_sha256, "portable manifest changed during load")
        return PortableReference(path, expected_sha256, runtime, validation)
    except (KeyError, TypeError, AttributeError) as exc:
        raise ValueError(f"malformed portable package: {exc}") from exc


def export_portable(local_manifest: Path, destination: Path, *, expected_sha256: str) -> Path:
    """Validate full local sources once, then copy the exact explicit allowlist."""
    source_path = Path(local_manifest).absolute()
    require(source_path.is_file() and sha256(source_path) == expected_sha256, "local manifest SHA256 mismatch")
    source = json.loads(source_path.read_text())
    validation = validate_package(source)
    _model(source)
    specs = _asset_specs(source)
    destination = Path(destination).absolute()
    require(not destination.exists() and not destination.is_symlink(), "refusing to overwrite portable package")
    destination.mkdir(parents=True)
    for spec in specs:
        output = _contained(destination, spec["path"])
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(spec["source_path"], output)
        _check_file(output, {k:spec[k] for k in ("path", "sha256", "bytes")})
    original = destination / "metadata/source-manifest.json"
    shutil.copyfile(source_path, original)
    require(sha256(original) == expected_sha256 == sha256(source_path), "local manifest changed during export")
    attestation = dict(scope="full_local_source_validation_at_export_not_reexecuted_at_load",
        source_manifest_sha256=expected_sha256, local_validation=validation,
        exporter_source=file_record(Path(__file__)), local_validator_source=file_record(Path(__file__).with_name("reference_package.py")),
        no_sequence_or_sketch_recomputation=True, no_external_distribution=True)
    attestation_path = destination / "metadata/export-validation.json"
    attestation_path.write_text(canonical_json(attestation))
    def relative_pin(file):return {**file_record(file), "path": file.relative_to(destination).as_posix()}
    manifest = dict(schema_version=SCHEMA, package_id=source["package_id"], package_version=source["package_version"],
        release_state=RELEASE, trust_basis=TRUST, positive_calling_policy_enabled=False,
        source_manifest=relative_pin(original), export_validation=relative_pin(attestation_path), assets=specs)
    result = destination / "portable.json"
    result.write_text(canonical_json(manifest))
    load_portable(result, expected_sha256=sha256(result))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    export = sub.add_parser("export")
    export.add_argument("--local-manifest", type=Path, required=True)
    export.add_argument("--destination", type=Path, required=True)
    export.add_argument("--expected-sha256", required=True)
    validate = sub.add_parser("validate")
    validate.add_argument("--manifest", type=Path, required=True)
    validate.add_argument("--expected-sha256", required=True)
    args = parser.parse_args()
    try:
        if args.command == "export":
            path = export_portable(args.local_manifest, args.destination, expected_sha256=args.expected_sha256)
            print(canonical_json({"manifest": str(path), "sha256": sha256(path)}), end="")
        else:
            print(canonical_json(load_portable(args.manifest, expected_sha256=args.expected_sha256).validation), end="")
    except (ValueError, OSError) as exc:
        parser.exit(2, f"portable reference error: {exc}\n")


if __name__ == "__main__":
    main()
