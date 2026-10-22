"""Supplementary assembly MLST against the pinned pre-2025-submission package.

No tool is executed here. A checked command plan and a caller-attested successful
result are separate records. ST is supplementary typing, never taxonomic support.
"""
from __future__ import annotations

from collections import Counter
import copy
import csv
import hashlib
import gzip
import json
from pathlib import Path
import re

from .evidence import explicit_percent, sha256

SCHEMA_VERSION = "1.0.0-dev-supplementary-mlst"
LOCI = ("dnaA", "fusA", "gyrB", "leuS", "pyrG", "rplB", "rpoB")
CUTOFF = "2024-12-31"
LIMITATION = (
    "Bundled MLST database: limited to PubMLST records submitted on or before "
    "31 December 2024. Alleles, ST profiles and schemes submitted from "
    "1 January 2025 onward are not included. Unmatched profiles may "
    "therefore correspond to STs described after this cutoff."
)
DATABASE_SUFFIXES = ("", ".ndb", ".nhr", ".nin", ".njs", ".nog", ".nos", ".not", ".nsq", ".ntf", ".nto")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _text(value: object) -> bool:
    return isinstance(value, str) and bool(value) and value == value.strip()


def _digest(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _canonical_sha(value: dict) -> str:
    try:
        data = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    except (TypeError, ValueError) as exc:
        raise ValueError("record must contain finite JSON values") from exc
    return hashlib.sha256(data).hexdigest()


def _fasta(path: Path, *, allele_locus: str | None = None) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier, parts = None, []
    def finish():
        if identifier is None:
            return
        sequence = "".join(parts).upper()
        _require(bool(sequence), "empty FASTA record")
        if allele_locus:
            match = re.fullmatch(re.escape(allele_locus) + r"_([0-9]+)", identifier)
            _require(match is not None and set(sequence) <= set("ACGT"), "invalid allele record/sequence")
            key = match[1]
        else:
            _require(set(sequence) <= set("ACGTRYSWKMBDHVN"), "invalid nucleotide assembly sequence")
            key = identifier
        _require(key not in records, "duplicate FASTA identifier")
        records[key] = sequence
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii", newline="") as handle:
        for number, raw in enumerate(handle, 1):
            _require(not raw.endswith("\r"), "bare CR line ending is not supported; use LF or CRLF")
            _require(number != 1 or raw.startswith(">"), "assembly/allele FASTA must start with a header at its first byte")
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                finish()
                fields = line[1:].split()
                _require(bool(fields), "empty FASTA identifier")
                identifier, parts = fields[0], []
            else:
                _require(identifier is not None, "sequence before FASTA header")
                parts.append(line)
    finish()
    _require(bool(records), "empty FASTA file")
    return records


class MlstContext:
    """Load checked data once; prepare and import explicit assembly-only calls.

    Version strings describe the tested caller contract. File hashes are checked;
    this class does not execute version probes or prove historical execution.
    Package and runtime pins are rechecked before preparation and interpretation.
    """

    def __init__(self, manifest_path: Path, *, expected_sha256: str):
        self._path = Path(manifest_path).resolve(strict=True)
        _require(_digest(expected_sha256) and sha256(self._path) == expected_sha256, "MLST manifest hash mismatch")
        self._manifest_sha256 = expected_sha256
        self._root = self._path.parent
        manifest = json.loads(self._path.read_text())
        _require(isinstance(manifest, dict) and type(manifest.get("schema_version")) is int and
                 manifest["schema_version"] == 1 and manifest.get("scheme") == "ecloacae", "unsupported MLST package")
        _require(manifest.get("loci") == list(LOCI) and manifest.get("submission_cutoff") == CUTOFF,
                 "unsupported locus order or submission cutoff")
        _require(manifest.get("limitation") == LIMITATION and _text(manifest.get("package_id")) and
                 _text(manifest.get("retrieved_utc")) and _text(manifest.get("sequence_state")), "missing cutoff/retrieval provenance")
        _require(manifest.get("complete_through_cutoff_independently_established") is False,
                 "this package does not establish exhaustive historical completeness")
        expected = {f"pubmlst/ecloacae/{l}.tfa" for l in LOCI} | {"pubmlst/ecloacae/ecloacae.txt"}
        hashes = manifest.get("data_sha256")
        _require(isinstance(hashes, dict) and set(hashes) == expected, "unexpected package data inventory")
        self._pins = {str(self._path): expected_sha256}
        for relative, digest in hashes.items():
            path = (self._root / relative).resolve(strict=True)
            _require(path.is_relative_to(self._root) and _digest(digest) and sha256(path) == digest,
                     "allele/profile data hash mismatch or external package path")
            self._pins[str(path)] = digest
        self._alleles = {l: _fasta(self._root / f"pubmlst/ecloacae/{l}.tfa", allele_locus=l) for l in LOCI}
        profiles = {}
        with (self._root / "pubmlst/ecloacae/ecloacae.txt").open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            _require(reader.fieldnames == ["ST", *LOCI, "clonal_complex"], "unexpected MLST profile columns")
            seen = set()
            for row in reader:
                _require(None not in row and all(v is not None for v in row.values()), "ragged profile row")
                st, signature = row["ST"], tuple(row[l] for l in LOCI)
                _require(re.fullmatch(r"[0-9]+", st) is not None and st not in profiles and signature not in seen,
                         "invalid/duplicate ST or duplicate profile")
                _require(all(row[l] in self._alleles[l] for l in LOCI), "profile references an absent allele")
                profiles[st] = signature
                seen.add(signature)
        _require(bool(profiles) and type(manifest.get("profile_count")) is int and len(profiles) == manifest["profile_count"] and
                 manifest.get("allele_counts") == {l: len(self._alleles[l]) for l in LOCI}, "package record count mismatch")
        self._profiles = profiles
        self._profile_to_st = {profile: st for st, profile in profiles.items()}
        self._manifest = copy.deepcopy(manifest)
        # This is the exact private FASTA convention tested with BLAST 2.16.0+.
        database = "".join(f">ecloacae.{l}_{a}\n{seq}\n" for l in LOCI for a, seq in self._alleles[l].items())
        self._database_fasta_sha256 = hashlib.sha256(database.encode()).hexdigest()
        self._check_pins(self._pins)

    @staticmethod
    def _check_pins(pins: dict[str, str]) -> None:
        for name, digest in pins.items():
            path = Path(name)
            _require(path.is_absolute() and path.is_file() and _digest(digest) and sha256(path) == digest,
                     f"MLST artifact hash mismatch: {name}")

    def _runtime(self, runtime: dict) -> dict[str, str]:
        _require(isinstance(runtime, dict) and set(runtime) == {"executables", "database"}, "explicit pinned MLST runtime required")
        executables = runtime["executables"]
        _require(isinstance(executables, dict) and set(executables) == {"mlst", "blastn", "any2fasta"}, "incomplete MLST toolchain")
        pins = {}
        for name, version in (("mlst", "2.32.2"), ("blastn", "2.16.0+"), ("any2fasta", None)):
            item = executables[name]
            keys = {"path", "sha256", "version"} if version else {"path", "sha256"}
            _require(isinstance(item, dict) and set(item) == keys and _text(item.get("path")), "malformed executable pin")
            _require(version is None or item["version"] == version, "untested MLST/BLAST version contract")
            _require(item["path"] not in pins, "toolchain executable paths must be distinct")
            pins[item["path"]] = item["sha256"]
        db = runtime["database"]
        _require(isinstance(db, dict) and set(db) == {"prefix", "sha256"} and _text(db["prefix"]), "explicit private BLAST database required")
        files = db["sha256"]
        expected = {db["prefix"] + suffix for suffix in DATABASE_SUFFIXES}
        _require(isinstance(files, dict) and set(files) == expected, "incomplete/unexpected tested BLAST index inventory")
        _require(files[db["prefix"]] == self._database_fasta_sha256, "BLAST source FASTA differs from verified package alleles")
        _require(not set(files) & set(pins), "database/toolchain path collision")
        pins.update(files)
        self._check_pins(pins)
        return pins

    def build_command(self, *, runtime: dict, query_path: Path, query_id: str, output_path: Path,
                      threads: int, min_identity_percent: object, min_coverage_percent: object) -> dict:
        """Prepare explicit argv/pins, without executing or claiming success.

        The caller must put the pinned executable directory on PATH for mlst's
        BLAST/any2fasta children. ``child_search_path`` identifies that directory.
        Matcher parameters are explicit experiment settings, not product defaults.
        Existing outputs are allowed for saved-result validation; an executing
        driver remains responsible for refusing an unintended overwrite.
        """
        self._check_pins(self._pins)
        pins = self._runtime(runtime)
        _require(_text(query_id) and type(threads) is int and threads > 0, "query identifier and positive thread count required")
        _require(not any(c in query_id for c in "\t\r\n"), "invalid query identifier")
        query = Path(query_path).resolve(strict=True)
        output = Path(output_path).absolute()
        _require(query != output and not output.is_dir(), "MLST output must be a file path distinct from query")
        _fasta(query)
        mi, mc = explicit_percent(min_identity_percent, "MLST minid"), explicit_percent(min_coverage_percent, "MLST mincov")
        _require(mi > 0 and mc > 0, "positive explicit matcher percentages required")
        directories = {str(Path(x["path"]).parent) for x in runtime["executables"].values()}
        _require(len(directories) == 1, "pinned MLST/BLAST/any2fasta executables must share a child PATH directory")
        argv = [runtime["executables"]["mlst"]["path"], "--scheme", "ecloacae", "--datadir", str(self._root / "pubmlst"),
                "--blastdb", runtime["database"]["prefix"], "--threads", str(threads), "--minid", f"{mi:g}",
                "--mincov", f"{mc:g}", "--full", "--debug", "--label", query_id, "--json", str(output), str(query)]
        return {"schema_version": SCHEMA_VERSION, "record_type": "prepared_command_not_executed", "input_kind": "assembly",
                "query_id": query_id, "query_path": str(query), "query_sha256": sha256(query), "output_path": str(output),
                "package_manifest_sha256": self._manifest_sha256, "runtime": copy.deepcopy(runtime),
                "runtime_sha256": pins, "child_search_path": directories.pop(), "argv": argv,
                "parameters": {"threads": threads, "min_identity_percent": mi, "min_coverage_percent": mc}}

    def _check_command(self, command: dict) -> None:
        _require(isinstance(command, dict) and command.get("schema_version") == SCHEMA_VERSION and
                 command.get("record_type") == "prepared_command_not_executed" and command.get("input_kind") == "assembly",
                 "explicit prepared assembly command required; reads are not implemented")
        _require(command.get("package_manifest_sha256") == self._manifest_sha256, "command/package digest mismatch")
        params = command.get("parameters")
        _require(isinstance(params, dict) and set(params) == {"threads", "min_identity_percent", "min_coverage_percent"},
                 "explicit MLST parameters required")
        # Reconstruct the fixed argv without touching the output; the driver
        # separately owns execution and output-creation safety.
        _require(_text(command.get("query_id")) and _text(command.get("query_path")) and _text(command.get("output_path")),
                 "command query/output fields missing")
        query = Path(command["query_path"])
        _require(query.is_absolute() and _digest(command.get("query_sha256")) and sha256(query) == command["query_sha256"],
                 "query hash mismatch")
        output = Path(command["output_path"])
        _require(output.is_absolute() and output != query and not output.is_dir(),
                 "invalid assembly/output path contract")
        _require(not any(c in command["query_id"] for c in "\t\r\n"), "invalid query identifier")
        _fasta(query)
        _require(sha256(query) == command["query_sha256"], "query changed during FASTA validation")
        _require(type(params["threads"]) is int and params["threads"] > 0, "invalid thread count")
        mi = explicit_percent(params["min_identity_percent"], "MLST minid")
        mc = explicit_percent(params["min_coverage_percent"], "MLST mincov")
        _require(mi > 0 and mc > 0, "positive matcher percentages required")
        self._check_pins(self._pins)
        pins = self._runtime(command.get("runtime"))
        _require(command.get("runtime_sha256") == pins, "command runtime inventory mismatch")
        runtime = command["runtime"]
        directories = {str(Path(x["path"]).parent) for x in runtime["executables"].values()}
        _require(directories == {command.get("child_search_path")}, "command child PATH mismatch")
        expected = [runtime["executables"]["mlst"]["path"], "--scheme", "ecloacae", "--datadir", str(self._root / "pubmlst"),
                    "--blastdb", runtime["database"]["prefix"], "--threads", str(params["threads"]), "--minid", f"{mi:g}",
                    "--mincov", f"{mc:g}", "--full", "--debug", "--label", command["query_id"], "--json", command["output_path"], str(query)]
        _require(command.get("argv") == expected, "MLST argv differs from explicit pinned command contract")

    def _locations(self, rows: list[dict] | None, calls: dict, command: dict) -> dict[str, list[dict]]:
        loci = {l: [] for l in LOCI}
        if rows is None:
            return loci
        _require(isinstance(rows, list), "exact locations must be a list or absent")
        genome = _fasta(Path(command["query_path"]))
        seen = set()
        for row in rows:
            _require(isinstance(row, dict), "location must be an object")
            locus, allele, contig = row.get("locus"), row.get("allele"), row.get("contig")
            _require(row.get("accession") == command["query_id"] and row.get("source_fasta_sha256") == command["query_sha256"],
                     "location/query identity or sequence hash mismatch")
            _require(isinstance(locus, str) and locus in LOCI and isinstance(allele, str) and allele in self._alleles[locus],
                     "unknown location locus/allele")
            _require(isinstance(contig, str) and contig in genome, "unknown location record")
            _require(re.fullmatch(r"[0-9]+(?:,[0-9]+)*", calls[locus]) is not None and allele in calls[locus].split(","),
                     "location not supported by an exact caller token")
            coords = []
            for key in ("start", "end"):
                value = row.get(key)
                _require(type(value) in (int, str) and re.fullmatch(r"[0-9]+", str(value)) is not None, "invalid exact location coordinate")
                coords.append(int(value))
            start, end = coords
            _require(1 <= start <= end <= len(genome[contig]), "exact location is out of bounds")
            key = locus, allele, contig, start, end
            _require(key not in seen, "duplicate exact location")
            seen.add(key)
            target = self._alleles[locus][allele]
            segment = genome[contig][start-1:end]
            _require(segment == target or segment == target.translate(str.maketrans("ACGT", "TGCA"))[::-1],
                     "exact location does not match the pinned allele sequence")
            loci[locus].append({"allele_id": allele, "record_id": contig, "start_1based": start, "end_1based_inclusive": end})
        for locus in LOCI:
            code = calls[locus]
            if re.fullmatch(r"[0-9]+(?:,[0-9]+)*", code):
                _require(Counter(code.split(",")) == Counter(r["allele_id"] for r in loci[locus]),
                         "exact location multiplicity disagrees with caller tokens")
        _require(sha256(Path(command["query_path"])) == command["query_sha256"], "query changed during location validation")
        return loci

    def interpret(self, record: dict, *, command: dict, search_exit_code: int,
                  exact_locations: list[dict] | None = None) -> dict:
        """Import one explicit successful result; never force a nearest ST.

        Optional location rows are checked against the hashed original sequence.
        Without them, repeated hits remain unconfirmed and cannot yield a copy-
        supported candidate ST. Supplying rows requires all reported exact hits.
        """
        _require(type(search_exit_code) is int and search_exit_code == 0, "failed search cannot become a negative MLST result")
        self._check_command(command)
        _require(isinstance(record, dict) and set(record) == {"id", "filename", "scheme", "sequence_type", "alleles"},
                 "unsupported or incomplete MLST JSON record")
        _require(record["id"] == command["query_id"] and record["filename"] == command["query_path"] and
                 record["scheme"] == "ecloacae", "result query/path/scheme mismatch")
        calls = record["alleles"]
        _require(isinstance(calls, dict) and set(calls) == set(LOCI), "missing/unexpected locus fields")
        for locus, code in calls.items():
            _require(isinstance(code, str) and re.fullmatch(r"-|[0-9]+(?:,[0-9]+)*|~[0-9]+|[0-9]+\?", code) is not None,
                     "unsupported allele token")
            _require(all(a in self._alleles[locus] for a in re.findall(r"[0-9]+", code)), "allele absent from pinned package")
        claimed = record["sequence_type"]
        _require(isinstance(claimed, str) and re.fullmatch(r"-|[0-9]+", claimed) is not None, "unsupported raw ST token; no novelty claims")
        codes = [calls[l] for l in LOCI]
        st = None
        if any("," in c for c in codes):
            state = "ambiguous_alleles"
        elif "-" in codes:
            state = "incomplete_profile"
        elif any("~" in c or "?" in c for c in codes):
            state = "inexact_alleles"
        else:
            st = self._profile_to_st.get(tuple(codes))
            state = "exact_profile_match" if st else "unmatched_profile"
        _require(claimed == (st or "-"), "raw ST conflicts with frozen exact-profile evidence")
        locations = self._locations(exact_locations, calls, command)
        repeated, discordant, same_sequence_ids, unconfirmed, signature, per_locus = [], [], [], [], [], []
        for locus in LOCI:
            code = calls[locus]
            numbers = re.findall(r"[0-9]+", code)
            exact = re.fullmatch(r"[0-9]+(?:,[0-9]+)*", code) is not None
            ids = set(numbers) if exact else set()
            signature.append(next(iter(ids)) if len(ids) == 1 else None)
            locus_state = "multiple_exact_alleles" if "," in code else "no_reported_match" if code == "-" else \
                          "inexact_full_length" if code.startswith("~") else "partial_match" if code.endswith("?") else "single_exact_allele"
            if "," in code:
                sequences = {self._alleles[locus][a] for a in ids}
                if len(sequences) > 1:
                    discordant.append(locus)
                elif len(ids) > 1:
                    same_sequence_ids.append(locus)
                elif len({(x["record_id"], x["start_1based"], x["end_1based_inclusive"]) for x in locations[locus]}) > 1:
                    repeated.append(locus)
                else:
                    unconfirmed.append(locus)
            per_locus.append({"locus": locus, "raw_token": code, "state": locus_state,
                              "reported_allele_ids": numbers, "reported_ids_are_exact": exact,
                              "sequence_verified_exact_locations": locations[locus]})
        candidate = self._profile_to_st.get(tuple(signature)) if repeated and not unconfirmed and all(x is not None for x in signature) else None
        reasons = [reason for values, reason in [(repeated, "repeated_identical_allele_copy"),
                   (discordant, "discordant_exact_allele_sequences"), (same_sequence_ids, "multiple_allele_ids_same_sequence"),
                   (unconfirmed, "repeated_hit_without_confirmed_multiple_locations")] if values]
        return {
            "schema_version": SCHEMA_VERSION, "record_type": "supplementary_mlst_evidence", "input_kind": "assembly",
            "query_id": command["query_id"], "query_path": command["query_path"], "query_sha256": command["query_sha256"],
            "scheme": "ecloacae", "state": state, "sequence_type": st, "alleles": copy.deepcopy(calls), "loci": per_locus,
            "profile_compatible_candidate_ST": candidate, "candidate_is_reported_ST": False,
            "copy_diagnostic_reasons": reasons, "repeated_same_allele_loci": repeated,
            "discordant_exact_allele_loci": discordant, "multiple_ids_same_sequence_loci": same_sequence_ids,
            "unconfirmed_repeated_hit_loci": unconfirmed,
            "database_package": self._manifest["package_id"], "package_manifest_sha256": self._manifest_sha256,
            "submission_cutoff": CUTOFF, "retrieved_utc": self._manifest["retrieved_utc"],
            "sequence_state": self._manifest["sequence_state"], "complete_through_cutoff_independently_established": False,
            "limitation": LIMITATION, "typing_role": "supplementary; does not establish genus, species, subspecies, lineage or genomovar",
            "taxonomic_assignment_performed": False, "scheme_routing": "explicitly_selected; applicability_not_established_by_ST",
            "execution": {"status": "caller_attested_successful_search", "search_exit_code": search_exit_code,
                          "command": copy.deepcopy(command), "raw_record_canonical_sha256": _canonical_sha(record),
                          "verification": "pinned_file_bytes_and_command_contract; historical_execution_not_reexecuted"},
            "location_verification": "supplied_exact_locations_checked_against_hashed_sequence" if exact_locations is not None else "not_supplied",
        }
