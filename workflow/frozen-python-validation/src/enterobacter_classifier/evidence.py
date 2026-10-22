"""Development evidence contract; deliberately contains no assignment policy.

Run with ``python -m enterobacter_classifier.evidence --help``. The historical
Task13 classifier stays separate. Only completed, assembly-derived skani output
is accepted here; future read adapters must define their own evidence semantics.
"""
from __future__ import annotations

import argparse
import copy
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

SCHEMA_VERSION = "2.0.0-dev-evidence"
REFERENCE_FIELDS = {
    "reference_id", "search_path", "species_id", "species_name", "subspecies_id",
    "subspecies_name", "lineage_id", "name_status", "reference_role",
}
SKANI_FIELDS = {"Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query"}
ROLES = {"type_anchor", "diversity", "boundary", "context_only"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_assembly(path: Path) -> None:
    """Check nucleotide FASTA syntax, without equating syntax with good QC."""
    opener = gzip.open if path.suffix == ".gz" else open
    seen = False
    length = 0
    identifiers: set[str] = set()
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                identifier = line[1:].split()[0] if line[1:].split() else ""
                if not identifier or identifier in identifiers or (seen and not length):
                    raise ValueError("assembly FASTA has empty/duplicate record ID or empty record")
                identifiers.add(identifier)
                seen, length = True, 0
            else:
                if not seen or not set(line.upper()) <= set("ACGTRYSWKMBDHVN"):
                    raise ValueError("input must be nucleotide assembly FASTA, not raw reads")
                length += len(line)
    if not seen or not length:
        raise ValueError("assembly FASTA is empty or ends with an empty record")


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not required <= set(fields):
            raise ValueError(f"{path}: missing or duplicated required TSV columns")
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(f"{path}: ragged TSV row")
    return rows


def explicit_percent(value: object, field: str) -> float:
    """Input is percent, never infer units from its magnitude (1 means 1%)."""
    if isinstance(value, bool):
        raise ValueError(f"{field}: boolean is not a percentage")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field}: missing/non-numeric percentage") from exc
    if not math.isfinite(number) or not 0 <= number <= 100:
        raise ValueError(f"{field}: percentage must be finite and in [0, 100]")
    return number


def validate_references(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Require explicit file mapping and stable concept IDs separate from names."""
    by_path: dict[str, dict[str, str]] = {}
    ids: set[str] = set()
    concepts: dict[tuple[str, str], tuple[str, str]] = {}
    if not rows:
        raise ValueError("reference inventory is empty")
    for row in rows:
        if not REFERENCE_FIELDS <= row.keys():
            raise ValueError("reference inventory columns missing")
        if any(v != v.strip() for v in row.values()):
            raise ValueError("reference inventory contains surrounding whitespace")
        rid, source = row["reference_id"], row["search_path"]
        if not rid or not source or rid in ids or source in by_path:
            raise ValueError("reference identifiers and search paths must be nonempty and unique")
        if row["reference_role"] not in ROLES or not row["name_status"]:
            raise ValueError("reference role/name status missing or invalid")
        for rank in ("species", "subspecies"):
            cid, name = row[f"{rank}_id"], row[f"{rank}_name"]
            if bool(cid) != bool(name):
                raise ValueError(f"{rid}: {rank} ID/name must both be present or absent")
            if not cid:
                continue
            parent = row["species_id"] if rank == "subspecies" else ""
            if rank == "subspecies" and not parent:
                raise ValueError(f"{rid}: subspecies requires a parent species")
            key, interpretation = (rank, cid), (parent, name)
            if key in concepts and concepts[key] != interpretation:
                raise ValueError(f"{rid}: conflicting concept parent/display name")
            concepts[key] = interpretation
        if not row["species_id"] and not row["lineage_id"] and row["reference_role"] != "context_only":
            raise ValueError(f"{rid}: assignable reference needs a named concept or explicit lineage")
        if row["lineage_id"]:
            key = ("lineage", row["lineage_id"])
            interpretation = (row["species_id"], row["subspecies_id"])
            if key in concepts and concepts[key] != interpretation:
                raise ValueError(f"{rid}: conflicting lineage parent")
            concepts[key] = interpretation
        ids.add(rid)
        by_path[source] = dict(row)
    return by_path


def collect_evidence(
    raw_hits: list[dict[str, str]], reference_rows: list[dict[str, str]], *,
    query_search_path: str, provenance: dict,
) -> dict:
    """Validate one query's returned estimates without filtering on ANI or AF.

    Missing rows are censored/not-returned, not zero ANI. Identical duplicate
    rows carry no extra support; conflicting duplicate estimates are an error.
    Processing failure must not be represented as an empty successful search.
    """
    references = validate_references(reference_rows)
    if not query_search_path or query_search_path != query_search_path.strip():
        raise ValueError("query search path must be explicit and nonempty")
    if type(provenance.get("search_exit_code")) is not int or provenance["search_exit_code"] != 0:
        raise ValueError("only a successfully completed search may produce an evidence record")
    if provenance.get("input_kind") != "assembly":
        raise ValueError("read evidence adapter is not implemented; do not use assembly ANI for reads")
    if not provenance.get("engine_version") or not provenance.get("search_parameters_record"):
        raise ValueError("engine version and search-parameter provenance are required")
    hits: dict[str, dict] = {}
    duplicates = 0
    for raw in raw_hits:
        if not SKANI_FIELDS <= raw.keys():
            raise ValueError("skani output missing required columns")
        if raw["Query_file"] != query_search_path:
            raise ValueError("output contains a different query; explicit input-path mapping required")
        if raw["Ref_file"] not in references:
            raise ValueError("output reference missing from inventory; do not infer ID from basename")
        ref = references[raw["Ref_file"]]
        item = {k: ref[k] for k in REFERENCE_FIELDS if k != "search_path"}
        item.update(
            ani_percent=explicit_percent(raw["ANI"], "ANI"),
            aligned_fraction_query_percent=explicit_percent(raw["Align_fraction_query"], "AF query"),
            aligned_fraction_reference_percent=explicit_percent(raw["Align_fraction_ref"], "AF reference"),
        )
        if item["aligned_fraction_query_percent"] == 0 or item["aligned_fraction_reference_percent"] == 0:
            raise ValueError("returned ANI cannot be supported by zero aligned fraction")
        item["minimum_directional_af_percent"] = min(
            item["aligned_fraction_query_percent"], item["aligned_fraction_reference_percent"])
        rid = ref["reference_id"]
        if rid in hits:
            if item != hits[rid]:
                raise ValueError("conflicting duplicate estimates for one query/reference")
            duplicates += 1
        hits[rid] = item
    ranked = sorted(hits.values(), key=lambda h: (-h["ani_percent"], h["reference_id"]))
    tied = [h["reference_id"] for h in ranked if h["ani_percent"] == ranked[0]["ani_percent"]] if ranked else []
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": "development_evidence_only",
        "input_kind": "assembly", "query_search_path": query_search_path,
        "provenance": copy.deepcopy(provenance),
        "ranking": "highest_returned_ani_no_coverage_filter; ties_retained; not_a_call",
        "units": {"ani": "percent", "aligned_fraction": "percent"},
        "evidence_status": "estimates_returned" if ranked else "no_estimate_returned",
        "raw_hit_count": len(raw_hits), "unique_hit_count": len(ranked),
        "duplicate_rows_collapsed": duplicates,
        "reference_count": len(references),
        "not_returned_reference_count": len(references) - len(ranked),
        "nearest_reference_ids": tied, "ranked_reference_evidence": ranked,
        "assignment": {
            "status": "not_evaluated", "reason": "calling_policy_pending",
            "species_id": None, "species_name": None,
            "subspecies_id": None, "subspecies_name": None,
            "lineage_id": None, "genomovar_id": None,
        },
        "qc": {"status": "not_assessed", "flags": []},
        "mlst": {"status": "not_run", "sequence_type": None},
        "limitations": [
            "Nearest reference is evidence, not an assignment or confidence estimate.",
            "A not-returned pair has no ANI estimate; search censoring depends on recorded parameters.",
            "No ANI/AF gate, taxonomic calling rule, QC rule or read classifier is enabled.",
            "This inventory may lack relevant taxa; absence of competitors does not establish identity.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--references", required=True, type=Path)
    parser.add_argument("--query", required=True, type=Path, help="Original assembly file, hashed as supplied")
    parser.add_argument("--query-search-path", required=True, help="Exact Query_file in saved skani output")
    parser.add_argument("--engine-version", required=True)
    parser.add_argument("--search-parameters-record", required=True, type=Path)
    parser.add_argument("--search-exit-code", required=True, type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite output")
    try:
        if args.query.resolve() != Path(args.query_search_path).resolve():
            raise ValueError("query path differs from saved search path; relocation mapping is not supported")
        validate_assembly(args.query)
        provenance = {
            "input_kind": "assembly", "sequence_scope": "assembly_as_supplied",
            "compartment_scope": "not_verified",
            "execution_provenance": "caller_attested_not_reexecuted",
            "engine": "skani", "engine_version": args.engine_version,
            "search_exit_code": args.search_exit_code,
            "search_parameters_record": str(args.search_parameters_record),
            "input_paths": {name: str(path.resolve()) for name, path in {
                "query": args.query, "hits": args.hits, "references": args.references,
                "search_parameters_record": args.search_parameters_record,
            }.items()},
            "sha256": {name: sha256(path) for name, path in {
                "query": args.query, "hits": args.hits, "references": args.references,
                "search_parameters_record": args.search_parameters_record,
            }.items()},
        }
        # This CLI audits supplied output. Hashes/path checks do not prove that
        # the stated historical search command/version was actually executed.
        result = collect_evidence(read_tsv(args.hits, SKANI_FIELDS),
                                  read_tsv(args.references, REFERENCE_FIELDS),
                                  query_search_path=args.query_search_path, provenance=provenance)
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
    except (OSError, ValueError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
