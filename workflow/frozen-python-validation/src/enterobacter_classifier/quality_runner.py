"""Execute explicitly pinned CheckM2 and retain evidence, without a QC gate.

Only the tested CheckM2 1.1.0/database3 contract is supported. The original
assembly is never rewritten. A separately recorded execution copy changes
only universal line endings and a missing terminal newline.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile
import time

from .evidence import read_tsv, sha256
from .quality import assess_assembly
from .quality_evidence import attach_checkm2, DATABASE_SHA256, NATIVE_STRUCTURE

CONFIG_SCHEMA = "1.0.0-dev-checkm2-runner-config"
NORMALIZATION_SCHEMA = "1.0.0-dev-checkm2-normalization"
NORMALIZATION_MODE = "ascii_universal_newlines_ensure_terminal_lf"
PACKAGE_NAMES = {"scikit-learn", "tensorflow", "lightgbm", "numpy", "pandas", "scipy", "keras"}


def _write(path: Path, data: dict | list) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _pinned(record: dict) -> Path:
    _require(isinstance(record, dict), "pinned file record required")
    path, digest = record.get("path"), record.get("sha256")
    _require(isinstance(path, str) and Path(path).is_absolute(), "pinned path must be absolute")
    _require(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
             "explicit SHA256 required")
    result = Path(path).resolve(strict=True)
    _require(result.is_file() and sha256(result) == digest, "pinned file SHA256 mismatch")
    return result


def _native_path(path: Path) -> None:
    # Upstream CheckM2 constructs unquoted Prodigal/DIAMOND shell commands.
    # Original query paths are never passed to those commands.
    _require(re.fullmatch(r"[A-Za-z0-9_./:-]+", str(path)) is not None,
             "CheckM2 native paths must contain no whitespace or shell metacharacters")


def normalize_assembly(source: Path, destination: Path) -> dict:
    """Create an exclusive copy with LF line endings and one terminal LF.

    Header text, sequence wrapping, case, bases and blank lines are preserved.
    ASCII universal-newline reading also supports an original .gz FASTA. This
    deterministic byte transformation can be independently repeated by the
    evidence importer; the equivalence is not merely a producer assertion.
    """
    source = Path(source).resolve(strict=True)
    destination = Path(destination).absolute()
    _require(source != destination, "normalization cannot replace the original assembly")
    before = assess_assembly(source)
    opener = gzip.open if source.suffix == ".gz" else open
    digest = hashlib.sha256()
    last_ended_with_lf = False
    with opener(source, "rt", encoding="ascii", newline=None) as inp, destination.open("xb") as out:
        for line in inp:
            data = line.encode("ascii")
            out.write(data)
            digest.update(data)
            last_ended_with_lf = line.endswith("\n")
        if not last_ended_with_lf:
            out.write(b"\n")
            digest.update(b"\n")
    _require(sha256(source) == before["query_sha256"], "original assembly changed during normalization")
    after = assess_assembly(destination)
    _require(after["query_sha256"] == digest.hexdigest(), "normalized copy changed")
    _require(before["structural_metrics"] == after["structural_metrics"]
             and before["record_metrics"] == after["record_metrics"], "normalization changed assembly structure")
    return {
        "schema_version": NORMALIZATION_SCHEMA, "mode": NORMALIZATION_MODE, "mode_version": "1",
        "original_path": str(source), "original_sha256": before["query_sha256"],
        "normalized_path": str(destination), "normalized_sha256": after["query_sha256"],
        "record_count": before["structural_metrics"]["record_count"],
        "total_bp": before["structural_metrics"]["total_bp"],
    }


def _manifest(path: Path, query: Path, digest: str) -> None:
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["input_id", "path", "sha256"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({"input_id": "query", "path": str(query), "sha256": digest})


def _command(argv: list[str], root: Path, name: str, env: dict, timeout: int) -> dict:
    started = time.monotonic()
    stdout, stderr = root / f"{name}.stdout", root / f"{name}.stderr"
    timed_out = False
    with stdout.open("x") as out, stderr.open("x") as err:
        process = subprocess.Popen(argv, stdout=out, stderr=err, env=env, start_new_session=True)
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        except BaseException:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise
    return {"argv": argv, "exit_status": process.returncode, "timed_out": timed_out,
            "elapsed_seconds": time.monotonic() - started,
            "stdout_path": str(stdout), "stdout_sha256": sha256(stdout),
            "stderr_path": str(stderr), "stderr_sha256": sha256(stderr)}


def run_checkm2(assembly: Path, config_path: Path, output_dir: Path, *,
               expected_config_sha256: str, threads: int = 2) -> Path:
    """Return one completed measured-quality envelope; failures retain logs.

    Require a new output directory and a byte-pinned configuration. No model or
    database is downloaded or changed. No numerical quality gate, purity claim,
    residual-sequence explanation or taxonomic assignment is produced.
    """
    _require(type(threads) is int and threads > 0, "threads must be a positive integer")
    assembly = Path(assembly).resolve(strict=True)
    config_path = Path(config_path).resolve(strict=True)
    output = Path(output_dir).resolve()
    _native_path(output)
    output.mkdir(parents=True, exist_ok=False)
    logs = output / "logs"
    logs.mkdir()
    state_path = output / "execution.json"
    started = time.monotonic()
    state = {"schema_version": "1.0.0-dev-checkm2-execution", "status": "running", "commands": [],
             "started_utc": datetime.now(timezone.utc).isoformat(), "threads": threads,
             "config_path": str(config_path), "config_sha256": expected_config_sha256}
    temp = None
    try:
        _pinned({"path": str(config_path), "sha256": expected_config_sha256})
        config = json.loads(config_path.read_text())
        _require(config.get("schema_version") == CONFIG_SCHEMA, "unsupported CheckM2 configuration")
        tool = _pinned(config["tool"])
        _require(config["tool"].get("version") == "1.1.0", "unsupported CheckM2 version")
        database = _pinned(config["database"])
        _require(config["database"].get("version") == "3" and config["database"]["sha256"] == DATABASE_SHA256,
                 "unsupported CheckM2 database")
        inventory = _pinned(config["installed_files_inventory"])
        installed = json.loads(inventory.read_text())
        _require(isinstance(installed, list) and installed and len({r["path"] for r in installed}) == len(installed),
                 "installed CheckM2 inventory must be nonempty and unique")
        for record in installed:
            _pinned(record)
        lock = _pinned(config["environment_lock"])
        dependencies = config["dependencies"]
        _require(set(dependencies) == {"python", "diamond", "prodigal"}, "explicit dependency executables required")
        executables = {name: _pinned(record) for name, record in dependencies.items()}
        _require(all(p.parent == tool.parent for p in executables.values()), "dependency executables must share the pinned environment")
        source_dir = output / "source"
        source_dir.mkdir()
        source_records = []
        for name in ("quality_runner.py", "quality.py", "quality_evidence.py", "evidence.py"):
            source = Path(__file__).resolve().parent / name
            snapshot = source_dir / name
            digest = sha256(source)
            shutil.copyfile(source, snapshot)
            _require(sha256(snapshot) == digest, "runner source changed during snapshot")
            source_records.append({"path": str(source), "sha256": digest,
                                   "snapshot_path": str(snapshot)})
        state["source_files"] = source_records
        packages = config["python_packages"]
        _require(isinstance(packages, dict) and set(packages) == PACKAGE_NAMES
                 and all(isinstance(v, str) and v for v in packages.values()), "explicit model dependency versions required")
        timeout = config["timeout_seconds"]
        _require(type(timeout) is int and timeout > 0, "positive timeout_seconds required")
        for path in (tool.parent, database):
            _native_path(path)
        temp_root = Path(config["temp_root"])
        _require(temp_root.is_absolute() and temp_root.is_dir(), "existing absolute temp_root required")
        _native_path(temp_root)
        _require(len(os.fsencode(str(temp_root))) <= 60, "temp_root is too long for CheckM2 Unix sockets")
        temp = Path(tempfile.mkdtemp(prefix="q", dir=temp_root))
        structural = assess_assembly(assembly)
        _write(output / "structural.json", structural)
        inputs = output / "execution-input"
        inputs.mkdir()
        normalized = normalize_assembly(assembly, inputs / "query.fna")
        normalization_path = output / "normalization.json"
        _write(normalization_path, normalized)
        manifest = output / "input-manifest.tsv"
        execution_manifest = output / "execution-input-manifest.tsv"
        _manifest(manifest, assembly, structural["query_sha256"])
        _manifest(execution_manifest, Path(normalized["normalized_path"]), normalized["normalized_sha256"])
        shutil.copyfile(config_path, output / "config.json")
        _require(sha256(output / "config.json") == expected_config_sha256, "configuration copy changed")
        state.update(query_sha256=structural["query_sha256"], input_manifest_sha256=sha256(manifest),
                     execution_input_manifest_sha256=sha256(execution_manifest), normalization_record_sha256=sha256(normalization_path),
                     database_sha256=DATABASE_SHA256, installed_files_inventory_sha256=sha256(inventory),
                     environment_lock_sha256=sha256(lock))
        env = {k: v for k, v in os.environ.items() if k not in {"PYTHONPATH", "PYTHONHOME"}}
        env.update(PATH=str(tool.parent) + os.pathsep + os.defpath, PYTHONNOUSERSITE="1", TMPDIR=str(temp),
                   OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                   TF_NUM_INTRAOP_THREADS=str(threads), TF_NUM_INTEROP_THREADS="1",
                   NUMEXPR_MAX_THREADS=str(threads), MPLCONFIGDIR=str(output / "cache/matplotlib"),
                   XDG_CACHE_HOME=str(output / "cache"))
        state["runtime_environment"] = {k: env[k] for k in ("PATH", "PYTHONNOUSERSITE", "TMPDIR", "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "TF_NUM_INTRAOP_THREADS", "TF_NUM_INTEROP_THREADS", "NUMEXPR_MAX_THREADS")}

        def command(name, args):
            record = _command(args, logs, name, env, timeout)
            state["commands"].append(record)
            _write(state_path, state)
            _require(record["exit_status"] == 0 and not record["timed_out"], f"CheckM2 step failed: {name}")
            return Path(record["stdout_path"]).read_text(), Path(record["stderr_path"]).read_text()

        out, _ = command("checkm2-version", [str(tool), "--version"])
        _require(out.strip() == "1.1.0", "CheckM2 runtime version mismatch")
        out, _ = command("diamond-version", [str(executables["diamond"]), "version"])
        _require(out.strip() == "diamond version " + dependencies["diamond"]["version"], "DIAMOND runtime version mismatch")
        out, err = command("prodigal-version", [str(executables["prodigal"]), "-v"])
        _require(any(line.startswith("Prodigal V" + dependencies["prodigal"]["version"] + ":")
                     for line in (out + err).splitlines()), "Prodigal runtime version mismatch")
        out, _ = command("python-version", [str(executables["python"]), "--version"])
        _require(out.strip() == "Python " + dependencies["python"]["version"], "Python runtime version mismatch")
        code = "import importlib.metadata,json;print(json.dumps({n:importlib.metadata.version(n) for n in " + repr(sorted(PACKAGE_NAMES)) + "},sort_keys=True))"
        out, _ = command("model-dependency-versions", [str(executables["python"]), "-c", code])
        _require(json.loads(out) == packages, "model dependency runtime versions mismatch")
        native_dir = output / "native"
        command("predict", [str(tool), "predict", "--input", str(inputs), "--extension", ".fna",
            "--output-directory", str(native_dir), "--threads", str(threads), "--lowmem", "--debug",
            "--database_path", str(database)])
        report = native_dir / "quality_report.tsv"
        rows = read_tsv(report, {"Name", "Completeness", "Contamination", "Completeness_Model_Used", *NATIVE_STRUCTURE})
        _require(len(rows) == 1 and rows[0]["Name"] == "query", "native quality output must contain exactly query")
        row = rows[0]
        complete, contamination = float(row["Completeness"]), float(row["Contamination"])
        _require(math.isfinite(complete) and 0 <= complete <= 100 and math.isfinite(contamination) and contamination >= 0,
                 "invalid native quality estimates")
        for key, metric in NATIVE_STRUCTURE.items():
            _require(re.fullmatch(r"[0-9]+", row[key]) is not None and int(row[key]) == structural["structural_metrics"][metric],
                     f"native/structural sequence mismatch: {key}")
        _require(sha256(assembly) == structural["query_sha256"] and sha256(Path(normalized["normalized_path"])) == normalized["normalized_sha256"],
                 "original or normalized query changed during execution")
        _require(sha256(config_path) == expected_config_sha256, "configuration changed during execution")
        for record in (config["tool"], config["database"], config["installed_files_inventory"],
                       config["environment_lock"], *dependencies.values(), *installed, *source_records):
            _pinned(record)
        for record in source_records:
            _pinned({"path": record["snapshot_path"], "sha256": record["sha256"]})
        state.update(status="passed", rows=1, report_sha256=sha256(report),
                     finished_utc=datetime.now(timezone.utc).isoformat(), elapsed_seconds=time.monotonic()-started)
        _write(state_path, state)
        envelope = dict(schema_version="1.0.0-dev-measured-quality-evidence", input_id="query", query_sha256=structural["query_sha256"],
            evidence_origin="local_calculation", status="estimated", tool="CheckM2", tool_version="1.1.0",
            completeness_percent=complete, contamination_percent=contamination, completeness_model=row["Completeness_Model_Used"],
            database_version="3", database_sha256=DATABASE_SHA256, native_output_path=str(report), native_output_sha256=sha256(report),
            native_row_key="query", native_row=row, input_manifest_path=str(manifest), input_manifest_sha256=sha256(manifest),
            input_manifest_fields={"identifier": "input_id", "path": "path", "sha256": "sha256"},
            installed_files_inventory_path=str(inventory), installed_files_inventory_sha256=sha256(inventory),
            execution_record_path=str(state_path), execution_record_sha256=sha256(state_path),
            normalization_record_path=str(normalization_path), normalization_record_sha256=sha256(normalization_path),
            execution_input_manifest_path=str(execution_manifest), execution_input_manifest_sha256=sha256(execution_manifest),
            config_path=str(output / "config.json"), config_sha256=expected_config_sha256,
            purity_claim=False, residual_sequence_explained=False, biological_quality_gate="not_defined", taxonomic_assignment="not_evaluated")
        pending = output / "quality-evidence.pending.json"
        _write(pending, envelope)
        attach_checkm2(structural, pending, expected_sha256=sha256(pending))
        result = output / "quality-evidence.json"
        pending.rename(result)
        return result
    except BaseException as exc:
        state.update(status="failed", error=f"{type(exc).__name__}: {exc}",
                     finished_utc=datetime.now(timezone.utc).isoformat(), elapsed_seconds=time.monotonic()-started)
        _write(state_path, state)
        raise
    finally:
        if temp is not None:
            shutil.rmtree(temp)
