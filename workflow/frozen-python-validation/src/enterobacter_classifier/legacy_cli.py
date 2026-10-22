from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

from . import __version__
from .core import Reference, classify_hits


REQUIRED_MANIFEST_COLUMNS = {"reference_id", "fasta_path", "species", "database_version", "nomenclature_version"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_db(args: argparse.Namespace) -> None:
    manifest = Path(args.manifest).resolve()
    output = Path(args.output).resolve()
    if output.exists():
        raise SystemExit(f"refusing to overwrite existing database: {output}")
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not REQUIRED_MANIFEST_COLUMNS.issubset(rows[0]):
        raise SystemExit("manifest requires: " + ", ".join(sorted(REQUIRED_MANIFEST_COLUMNS)))
    manifest_fields = list(rows[0])
    if "sha256" in manifest_fields:
        raise SystemExit("manifest must not supply sha256; build-db derives it from copied FASTAs")
    identifiers = [row["reference_id"].strip() for row in rows]
    versions = {(row["database_version"].strip(), row["nomenclature_version"].strip()) for row in rows}
    if len(identifiers) != len(set(identifiers)) or not all(identifiers) or len(versions) != 1 or not all(next(iter(versions))):
        raise SystemExit("reference IDs must be unique; database and nomenclature versions must be one non-empty pair")
    output.mkdir(parents=True)
    copied = []
    for row in sorted(rows, key=lambda value: value["reference_id"]):
        source = Path(row["fasta_path"]).resolve()
        if not source.is_file():
            raise SystemExit(f"reference FASTA missing: {source}")
        suffix = source.suffix.lower() if source.suffix.lower() in {".fa", ".fna", ".fasta", ".fas"} else ".fna"
        target = output / "reference-inputs" / f"{row['reference_id']}{suffix}"
        target.parent.mkdir(exist_ok=True)
        shutil.copyfile(source, target)
        row["fasta_path"] = str(target.relative_to(output))
        row["sha256"] = sha256(target)
        copied.append(row)
    fields = manifest_fields + ["sha256"]
    with (output / "references.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader(); writer.writerows(copied)
    database_version, nomenclature_version = next(iter(versions))
    metadata = {"database_schema_version": "1.0.0", "builder_version": __version__,
                "database_version": database_version, "nomenclature_version": nomenclature_version,
                "reference_count": len(copied), "references_sha256": sha256(output / "references.tsv")}
    (output / "database.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metadata, sort_keys=True))


def assign(args: argparse.Namespace) -> None:
    database = Path(args.database).resolve(); query = Path(args.query).resolve(); output = Path(args.output).resolve()
    if not query.is_file() or not (database / "database.json").is_file(): raise SystemExit("query or database missing")
    if output.exists(): raise SystemExit(f"refusing to overwrite existing output: {output}")
    config = json.loads(Path(args.config).read_text())
    rows = list(csv.DictReader((database / "references.tsv").open(), delimiter="\t"))
    references = {row["reference_id"]: Reference(row["reference_id"], row["species"], row.get("subspecies", ""), row.get("operational_cluster_id", "")) for row in rows}
    output.mkdir(parents=True)
    raw_hits = output / "skani.tsv"
    if args.hits:
        shutil.copyfile(args.hits, raw_hits)
    else:
        command = [args.skani, "dist", "-q", str(query), "-r"] + [str(database / row["fasta_path"]) for row in rows] + ["-t", str(args.threads), "-o", str(raw_hits)]
        (output / "command.json").write_text(json.dumps(command) + "\n")
        subprocess.run(command, check=True)
    with raw_hits.open(newline="") as handle:
        parsed = list(csv.DictReader(handle, delimiter="\t"))
    hits = [{"reference_id": Path(row.get("Ref_file", row.get("reference_id", ""))).stem,
             "ani": row.get("ANI", row.get("ani")), "af_ref": row.get("Align_fraction_ref", row.get("af_ref")),
             "af_query": row.get("Align_fraction_query", row.get("af_query"))} for row in parsed]
    result = classify_hits(hits, references, config)
    result["database"] = json.loads((database / "database.json").read_text())
    result["query_sha256"] = sha256(query)
    (output / "assignment.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser(prog="enterobacter_classifier")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    build = sub.add_parser("build-db"); build.add_argument("--manifest", required=True); build.add_argument("--output", required=True); build.set_defaults(func=build_db)
    assigner = sub.add_parser("assign")
    for name in ("query", "database", "config", "output"): assigner.add_argument(f"--{name}", required=True)
    assigner.add_argument("--threads", type=int, default=4); assigner.add_argument("--skani", default="skani"); assigner.add_argument("--hits")
    assigner.set_defaults(func=assign)
    args = parser.parse_args(); args.func(args)


if __name__ == "__main__": main()
