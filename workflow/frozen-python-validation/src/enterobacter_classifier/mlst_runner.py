"""Execute the explicitly pinned supplementary MLST contract on one assembly.

The caller selects this typing scheme; neither ST nor successful execution proves
scheme applicability, genus, species, biological quality or a taxonomic rank.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time

from .evidence import sha256
from .mlst import LOCI, MlstContext

CONFIG_SCHEMA = "1.0.0-dev-mlst-runner"
EXECUTION_SCHEMA = "1.0.0-dev-mlst-execution"
_HIT = re.compile(r"^\[(\d+)\] (\S+):(\d+)-(\d+)\((plus|minus)\) \| ecloacae (\w+) (\d+) \| id=(\d+)/(\d+) \| cov=(\d+)/(\d+) \| seq=(\S+)$")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _write(path: Path, value: object) -> None:
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def parse_exact_locations(debug_text: str, *, command: dict) -> list[dict]:
    """Parse the tested mlst2.32.2 debug grammar, without asserting execution.

    This extracts all full exact hits, including duplicate allele tokens. It
    rejects malformed hit rows and duplicate locations; it does not pick a hit
    subset to make an ST match. MlstContext.interpret must then check query and
    runtime pins, reported token multiplicities and original allele sequences.
    """
    _require(isinstance(debug_text, str) and isinstance(command, dict), "debug text and prepared command required")
    _require(command.get("record_type") == "prepared_command_not_executed" and
             command.get("runtime", {}).get("executables", {}).get("mlst", {}).get("version") == "2.32.2" and
             command.get("runtime", {}).get("executables", {}).get("blastn", {}).get("version") == "2.16.0+" and
             "--debug" in command.get("argv", []), "unsupported native debug contract")
    query_id, digest = command.get("query_id"), command.get("query_sha256")
    _require(isinstance(query_id, str) and bool(query_id) and isinstance(digest, str) and
             re.fullmatch(r"[0-9a-f]{64}", digest) is not None, "debug query identity/hash required")
    rows, seen, hit_numbers = [], set(), set()
    for line in debug_text.splitlines():
        match = _HIT.fullmatch(line)
        if not match:
            _require(not (re.match(r"^\[\d+\]", line) or "| id=" in line or "| cov=" in line),
                     "malformed native MLST alignment record")
            continue
        number, contig, start, end, strand, locus, allele, nident, alen, covlen, hlen, sequence = match.groups()
        number, start, end, nident, alen, covlen, hlen = map(int, [number, start, end, nident, alen, covlen, hlen])
        _require(number > 0 and number not in hit_numbers, "duplicate/invalid native hit number")
        hit_numbers.add(number)
        _require(locus in LOCI and 0 < nident <= alen <= hlen and alen == covlen == len(sequence) and
                 abs(end-start)+1 == alen and min(start, end) > 0 and set(sequence.upper()) <= set("ACGTRYSWKMBDHVN"),
                 "inconsistent native MLST alignment fields")
        if nident != alen or alen != hlen:
            continue
        start, end = sorted([start, end])
        key = locus, allele, contig, start, end
        _require(key not in seen, "duplicate exact native location")
        seen.add(key)
        rows.append(dict(accession=query_id, locus=locus, allele=allele, contig=contig,
                         start=start, end=end, source_fasta_sha256=digest))
    return rows


def _command(argv: list[str], *, output: Path, name: str, env: dict) -> dict:
    started, clock = datetime.now(timezone.utc).isoformat(), time.monotonic()
    stdout, stderr = output / f"{name}.stdout", output / f"{name}.stderr"
    record = dict(argv=argv, cwd=str(Path.cwd()), started_utc=started,
                  child_PATH_prefix=env["PATH"].split(os.pathsep)[0], temp_root=env.get("TMPDIR"),
                  stdout_path=str(stdout), stderr_path=str(stderr), exit_code=None)
    try:
        with stdout.open("x") as out, stderr.open("x") as err:
            completed = subprocess.run(argv, stdout=out, stderr=err, env=env, check=False)
        record["exit_code"] = completed.returncode
    except BaseException as exc:
        record["launch_error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        record["elapsed_seconds"] = time.monotonic() - clock
        for label, path in (("stdout", stdout), ("stderr", stderr)):
            if path.is_file():
                record[label + "_sha256"] = sha256(path)
        _write(output / f"{name}.command.json", record)
    return record


def run_mlst(query: Path, config_path: Path, output_dir: Path, *,
             expected_config_sha256: str, threads: int) -> dict:
    """Execute one assembly, returning checked supplementary typing evidence.

    A new output directory is required. Failure logs and execution state remain
    inspectable and no result is produced after a failed command or import. The
    query filename is an execution label only, never a biological identity.
    """
    _require(type(threads) is int and threads > 0, "positive explicit thread count required")
    output = Path(output_dir).absolute()
    _require(not output.exists() and not output.is_symlink(), "refusing to overwrite MLST output directory")
    _require(output.parent.is_dir(), "MLST output parent must exist")
    output.mkdir(exist_ok=False)
    logs = output / "logs"
    logs.mkdir()
    started = time.monotonic()
    state = dict(schema_version=EXECUTION_SCHEMA, status="running", commands=[],
                 started_utc=datetime.now(timezone.utc).isoformat(), threads=threads,
                 config_path=str(Path(config_path).absolute()), expected_config_sha256=expected_config_sha256)
    try:
        query, config_path = Path(query).resolve(strict=True), Path(config_path).resolve(strict=True)
        _require(isinstance(expected_config_sha256, str) and
                 re.fullmatch(r"[0-9a-f]{64}", expected_config_sha256) is not None and
                 sha256(config_path) == expected_config_sha256, "MLST configuration hash mismatch")
        config = json.loads(config_path.read_text())
        _require(isinstance(config, dict) and set(config) == {"schema_version", "package", "runtime", "matcher"} and
                 config["schema_version"] == CONFIG_SCHEMA, "unsupported MLST runner configuration")
        package = config["package"]
        _require(isinstance(package, dict) and set(package) == {"path", "sha256"} and
                 isinstance(package["path"], str) and Path(package["path"]).is_absolute(), "explicit package pin required")
        matcher = config["matcher"]
        _require(isinstance(matcher, dict) and set(matcher) == {"min_identity_percent", "min_coverage_percent"},
                 "explicit MLST matcher settings required")
        context = MlstContext(Path(package["path"]), expected_sha256=package["sha256"])
        plan = context.build_command(runtime=config["runtime"], query_path=query, query_id=query.name,
                                     output_path=output / "raw.json", threads=threads, **matcher)
        _require(sha256(config_path) == expected_config_sha256, "MLST configuration changed during preflight")
        shutil.copyfile(config_path, output / "config.json")
        _require(sha256(output / "config.json") == expected_config_sha256, "MLST config snapshot differs from original bytes")
        _write(output / "command-plan.json", plan)
        state.update(query_path=str(query), query_sha256=plan["query_sha256"], package_manifest_sha256=package["sha256"])
        source_dir = output / "source"
        source_dir.mkdir()
        state["source_files"] = []
        for name in ("mlst_runner.py", "mlst.py", "evidence.py"):
            source = Path(__file__).resolve().parent / name
            snapshot = source_dir / name
            digest = sha256(source)
            shutil.copyfile(source, snapshot)
            _require(sha256(snapshot) == digest, "MLST source changed during snapshot")
            state["source_files"].append(dict(path=str(source), sha256=digest, snapshot_path=str(snapshot)))
        env = dict(os.environ)
        env["PATH"] = plan["child_search_path"] + os.pathsep + env.get("PATH", "")
        for name, item in config["runtime"]["executables"].items():
            _require(Path(item["path"]).name == name and
                     Path(shutil.which(name, path=env["PATH"]) or "").resolve() == Path(item["path"]).resolve(),
                     "child executable resolution differs from pinned MLST runtime")
        # Upstream writes its temporary FASTA/BLAST paths unquoted in shell
        # strings. A private safe temporary root avoids inheriting unsafe TMPDIR.
        with tempfile.TemporaryDirectory(prefix="enterobacter-mlst-", dir="/tmp") as temp:
            env.update(TMPDIR=temp, TMP=temp, TEMP=temp)
            for name, args, expected in (("mlst", ["--version"], "mlst 2.32.2"),
                                         ("blastn", ["-version"], "blastn: 2.16.0+")):
                item = _command([config["runtime"]["executables"][name]["path"], *args],
                                output=logs, name=f"{name}-version", env=env)
                state["commands"].append(item)
                _require(item["exit_code"] == 0 and Path(item["stdout_path"]).read_text().splitlines()[:1] == [expected],
                         f"{name} version probe failed or differs from pinned contract")
            # A version process must not alter any pinned query/package/runtime.
            context._check_command(plan)
            item = _command(plan["argv"], output=logs, name="mlst", env=env)
            state["commands"].append(item)
            _require(item["exit_code"] == 0, f"MLST search failed with exit code {item['exit_code']}; see retained logs")
        raw_path = output / "raw.json"
        _require(raw_path.is_file(), "successful MLST command did not produce raw JSON")
        raw_digest = sha256(raw_path)
        debug_path = logs / "mlst.stderr"
        _require(sha256(debug_path) == item["stderr_sha256"], "MLST debug log changed after execution")
        raw = json.loads(raw_path.read_text())
        _require(isinstance(raw, list) and len(raw) == 1, "single-query MLST JSON must contain exactly one record")
        locations = parse_exact_locations(debug_path.read_text(), command=plan)
        result = context.interpret(raw[0], command=plan, search_exit_code=item["exit_code"], exact_locations=locations)
        _require(sha256(config_path) == expected_config_sha256, "MLST configuration changed during execution")
        _require(sha256(raw_path) == raw_digest and sha256(debug_path) == item["stderr_sha256"],
                 "MLST raw/debug output changed during interpretation")
        _write(output / "exact-locations.json", locations)
        result["execution"].update(status="executed_and_checked", verification="runner_observed_exit_zero_and_pinned_artifacts; exact_locations_checked_against_original_sequence",
                                   runner_execution_path=str(output / "execution.json"), config_sha256=expected_config_sha256,
                                   raw_json_path=str(raw_path), raw_json_sha256=sha256(raw_path),
                                   debug_stderr_path=str(logs / "mlst.stderr"), debug_stderr_sha256=sha256(logs / "mlst.stderr"))
        _write(output / "result.json", result)
        state.update(status="complete", result_path=str(output / "result.json"), result_sha256=sha256(output / "result.json"),
                     exact_locations_path=str(output / "exact-locations.json"), exact_locations_sha256=sha256(output / "exact-locations.json"),
                     taxonomic_assignment_performed=False)
        return result
    except BaseException as exc:
        state.update(status="failed", error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        state["elapsed_seconds"] = time.monotonic() - started
        state["finished_utc"] = datetime.now(timezone.utc).isoformat()
        _write(output / "execution.json", state)
