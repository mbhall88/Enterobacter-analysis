"""Run the unified development assembly workflow; no positive policy is enabled.

One command validates an assembly, runs the package-pinned skani binary, retains
raw evidence and produces the existing uncertainty report. Structural validity
does not establish biological quality. Reads and final assignment rules are not
implemented. Original assembly records are neither selected nor discarded.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import resource
import shutil
import subprocess
import time
from pathlib import Path

from .evidence import SKANI_FIELDS, collect_evidence, read_tsv, sha256
from .quality import assess_assembly
from .quality_evidence import attach_checkm2
from .quality_runner import run_checkm2
from .rank import evaluate_rank, validate_rule
from .eligibility import evaluate_eligibility, primary_inventory, validate_config as validate_eligibility_config
from .genomovar import GenomovarContext
from .reference_package import canonical_json
from .reporting import ReportingContext, write_reports

SCHEMA = "1.0.0-dev-unified-workflow"


def _write(path: Path, data: dict) -> None:
    with path.open("x") as handle:
        handle.write(canonical_json(data))


def _run(argv: list[str], output: Path, name: str) -> dict:
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    clock = time.monotonic()
    with (output / f"{name}.stdout").open("x") as stdout, (output / f"{name}.stderr").open("x") as stderr:
        try:
            completed = subprocess.run(argv, stdout=stdout, stderr=stderr, check=False)
            code = completed.returncode
        except OSError as exc:
            _write(output / f"{name}.command.json", dict(argv=argv, cwd=str(Path.cwd()),
                   started_utc=started, elapsed_seconds=time.monotonic()-clock,
                   exit_code=None, launch_error=str(exc)))
            raise
    record = dict(argv=argv, cwd=str(Path.cwd()), started_utc=started,
                  elapsed_seconds=time.monotonic()-clock, exit_code=code,
                  stdout_path=str(output / f"{name}.stdout"), stderr_path=str(output / f"{name}.stderr"))
    _write(output / f"{name}.command.json", record)
    if code != 0:
        raise ValueError(f"{name} failed with exit code {code}; inspect retained logs")
    return record


def run_assembly(*, query: Path, reference_package: Path,
                 expected_package_sha256: str, output: Path, threads: int = 1,
                 development_guides: dict | None = None,
                 quality_evidence_path: Path | None = None,
                 expected_quality_sha256: str | None = None,
                 experimental_rank_rule: dict | None = None,
                 checkm2_config: Path | None = None,
                 expected_checkm2_config_sha256: str | None = None,
                 experimental_genomovars: bool = False,
                 mlst_config: Path | None = None,
                 expected_mlst_config_sha256: str | None = None,
                 experimental_eligibility_config: Path | None = None,
                 expected_eligibility_config_sha256: str | None = None,
                 skani_binary: Path | None = None) -> dict:
    """Execute a new local run, rejecting output replacement and unpinned inputs.

    The search reporting threshold (15% on either side) is inherited from the
    frozen skani comparison, not a taxonomic AF requirement. No sequence input
    enters this function as a final independent-test sample implicitly.

    An explicit eligibility configuration adds the disabled built-in screen of
    the frozen candidate policy (rank evidence and the two-parent alert; any QC
    estimate is reported only). Optional ANIm and
    warning envelopes bind the saved evidence digest, so they belong to the
    post-analysis ``eligibility.EligibilityContext`` rather than this run.
    """
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    if (quality_evidence_path is None) != (expected_quality_sha256 is None):
        raise ValueError("quality evidence path and expected SHA256 must be supplied together")
    if (checkm2_config is None) != (expected_checkm2_config_sha256 is None):
        raise ValueError("automatic CheckM2 configuration and expected SHA256 must be supplied together")
    if checkm2_config is not None and quality_evidence_path is not None:
        raise ValueError("select automatic quality execution or a pinned completed estimate, not both")
    if (mlst_config is None) != (expected_mlst_config_sha256 is None):
        raise ValueError("MLST configuration and expected SHA256 must be supplied together")
    if type(experimental_genomovars) is not bool:
        raise ValueError("experimental genomovar selection must be boolean")
    if experimental_genomovars and experimental_rank_rule is None:
        raise ValueError("targeted genomovar investigation requires an explicit experimental rank rule")
    if (experimental_eligibility_config is None) != (expected_eligibility_config_sha256 is None):
        raise ValueError("experimental eligibility configuration and expected SHA256 must be supplied together")
    if experimental_eligibility_config is not None and experimental_rank_rule is None:
        raise ValueError("experimental eligibility requires an explicit experimental rank rule")
    if experimental_rank_rule is not None:
        validate_rule(experimental_rank_rule)
        experimental_rank_rule = json.loads(canonical_json(experimental_rank_rule))
    eligibility_config = None
    if experimental_eligibility_config is not None:
        experimental_eligibility_config = Path(experimental_eligibility_config).absolute()
        if sha256(experimental_eligibility_config) != expected_eligibility_config_sha256:
            raise ValueError("experimental eligibility configuration SHA256 mismatch")
        eligibility_config = json.loads(experimental_eligibility_config.read_text())
        validate_eligibility_config(eligibility_config, experimental_rank_rule)
    query, reference_package, output = query.absolute(), reference_package.absolute(), output.absolute()
    if output.exists() or output.is_symlink():
        raise ValueError("refusing to overwrite output directory")
    if not output.parent.is_dir():
        raise ValueError("output parent must already exist")
    structural = assess_assembly(query)
    quality = (attach_checkm2(structural, quality_evidence_path, expected_sha256=expected_quality_sha256)
               if quality_evidence_path is not None else structural)
    # Conservative workflow input guard, informed by pinned skani0.3.2:
    # records shorter than500bp are skipped; DNA screening uses21-mers.
    # Ambiguity handling differs by CPU path; this is not exact sketch emulation
    # or a guarantee that enough informative anchors exist.
    eligible = [r for r in quality["record_metrics"]
                if r["length_bp"] >= 500 and r["max_unambiguous_acgt_run_bp"] >= 21]
    if not eligible:
        raise ValueError("unsupported by development input guard: no record is at least500bp with a contiguous21-base ACGT run; no taxonomic search performed")
    length_eligible = [r for r in quality["record_metrics"] if r["length_bp"] >= 500]
    query_engine_scope = dict(engine="skani", engine_version="0.3.2",
        records_passed_to_engine=quality["structural_metrics"]["record_count"],
        skani_length_eligible_record_count=len(length_eligible),
        skani_length_ignored_record_count=quality["structural_metrics"]["record_count"]-len(length_eligible),
        skani_length_denominator_bp=sum(r["length_bp"] for r in length_eligible),
        skani_length_ignored_bp=sum(r["length_bp"] for r in quality["record_metrics"] if r["length_bp"] < 500),
        interpretation="Input unchanged; skani skips records shorter than500bp. Length eligibility does not ensure informative sketch anchors.")
    context = ReportingContext(reference_package, expected_sha256=expected_package_sha256)
    gv_context = GenomovarContext(reference_package, expected_sha256=expected_package_sha256) if experimental_genomovars else None
    # Validate the optional display configuration before running an executable.
    context._investigation([], development_guides)
    manifest = context._manifest
    is_portable = context.validation.get("original_fastas_required_at_load") is False
    if is_portable and skani_binary is None:
        raise ValueError("portable package requires an explicit local skani binary path")
    if sha256(reference_package) != expected_package_sha256:
        raise ValueError("reference package changed during preflight")
    component = next(c for c in manifest["components"] if c["component_id"] == "primary_whole")
    binary_record = context._package_view.engine_record("primary_whole", skani_binary)
    binary = Path(binary_record["path"])
    files = component["database_files"]
    dirs = {Path(f["path"]).parent for f in files}
    if len(dirs) != 1 or component["engine"] != "skani" or component["engine_version"] != "0.3.2":
        raise ValueError("workflow requires the frozen skani0.3.2 primary component")
    database = dirs.pop()
    pins = {str(query): quality["query_sha256"], str(reference_package): expected_package_sha256,
            str(binary): binary_record["sha256"]}
    pins.update({f["path"]: f["sha256"] for f in files})
    gv_database = None
    if gv_context is not None:
        gv_component = next(c for c in manifest["components"] if c["component_id"] == "genomovar_chromosome")
        gv_binary = context._package_view.engine_record("genomovar_chromosome", skani_binary)
        if gv_binary != binary_record or gv_component["engine_version"] != "0.3.2":
            raise ValueError("genomovar workflow requires the same pinned skani executable")
        gv_dirs = {Path(f["path"]).parent for f in gv_component["database_files"]}
        if len(gv_dirs) != 1:
            raise ValueError("genomovar database files must share a directory")
        gv_database = gv_dirs.pop()
        pins.update({f["path"]: f["sha256"] for f in gv_component["database_files"]})
    if is_portable:
        # The wrapper digest pins the inventory; keep copied metadata as well
        # as sketches stable across the entire execution window.
        portable = json.loads(reference_package.read_text())
        for item in [portable["source_manifest"], portable["export_validation"], *portable["assets"]]:
            pins[str(reference_package.parent / item["path"])] = item["sha256"]
    if checkm2_config is not None:
        checkm2_config = Path(checkm2_config).absolute()
        pins[str(checkm2_config)] = expected_checkm2_config_sha256
    if mlst_config is not None:
        mlst_config = Path(mlst_config).absolute()
        pins[str(mlst_config)] = expected_mlst_config_sha256
    if experimental_eligibility_config is not None:
        pins[str(experimental_eligibility_config)] = expected_eligibility_config_sha256
    for path, digest in pins.items():
        if sha256(Path(path)) != digest:
            raise ValueError(f"input changed during preflight: {path}")
    output.mkdir()
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    clock = time.monotonic()
    try:
        _write(output / "input-pins.json", pins)
        if experimental_rank_rule is not None:
            _write(output / "experimental-rank-rule.json", experimental_rank_rule)
        if eligibility_config is not None:
            _write(output / "experimental-eligibility-config.json", eligibility_config)
        _write(output / "query-engine-scope.json", query_engine_scope)
        source_pins = {str(p): sha256(p) for p in sorted(Path(__file__).parent.glob("*.py"))}
        _write(output / "source-pins.json", source_pins)
        (output / "source").mkdir()
        for path, digest in source_pins.items():
            snapshot = output / "source" / Path(path).name
            shutil.copyfile(path, snapshot)
            if sha256(snapshot) != digest:
                raise ValueError("workflow source changed during provenance capture")
        _write(output / "configuration.json", dict(positive_calling_policy_enabled=False,
            biological_scope_detection_enabled=False, input_kind="assembly", threads=threads,
            development_display_guides=development_guides, search_reporting_min_af_percent=15,
            imported_quality_evidence_sha256=expected_quality_sha256,
            automatic_checkm2_config_sha256=expected_checkm2_config_sha256,
            quality_execution_mode="automatic_checkm2" if checkm2_config is not None else
                                   "imported_estimate" if quality_evidence_path is not None else "structural_only",
            quality_estimate_interpretation="no_quality_gate_or_purity_inference",
            experimental_rank_evaluation_enabled=experimental_rank_rule is not None,
            experimental_genomovar_evaluation_enabled=experimental_genomovars,
            experimental_eligibility_evaluation_enabled=eligibility_config is not None,
            experimental_eligibility_config_sha256=expected_eligibility_config_sha256,
            supplementary_mlst_config_sha256=expected_mlst_config_sha256,
            approximate_search_screen_percent=80, query_records="all_input_records_unchanged",
            records_passing_development_input_guard=len(eligible),
            development_input_guard="conservative_not_exact_skani_emulation; record_length_ge500_and_ACGT_run_ge21; no_guarantee_of_informative_sketch"))
        if checkm2_config is not None:
            quality_evidence_path = run_checkm2(query, checkm2_config, output / "quality-run",
                expected_config_sha256=expected_checkm2_config_sha256, threads=threads)
            expected_quality_sha256 = sha256(quality_evidence_path)
            quality = attach_checkm2(structural, quality_evidence_path, expected_sha256=expected_quality_sha256)
        _write(output / "quality.json", quality)
        if quality_evidence_path is not None:
            shutil.copyfile(quality_evidence_path, output / "imported-quality-evidence.json")
            if sha256(output / "imported-quality-evidence.json") != expected_quality_sha256:
                raise ValueError("quality evidence changed during provenance capture")
        _run([str(binary), "--version"], output, "skani-version")
        if (output / "skani-version.stdout").read_text().strip() != "skani 0.3.2":
            raise ValueError("pinned binary reports an unexpected skani version")
        raw = output / "hits.tsv"
        argv = [str(binary), "search", "-q", str(query), "-d", str(database),
                "-o", str(raw), "--min-af", "15", "-s", "80", "-t", str(threads)]
        command = _run(argv, output, "skani-search")
        for path, digest in pins.items():
            if sha256(Path(path)) != digest:
                raise ValueError(f"input changed during execution: {path}")
        if quality_evidence_path is not None:
            checked_quality = attach_checkm2(structural, quality_evidence_path, expected_sha256=expected_quality_sha256)
            if checked_quality != quality:
                raise ValueError("quality evidence changed during execution")
        provenance = dict(input_kind="assembly", sequence_scope="whole_assembly_as_supplied",
            compartment_scope="not_verified", engine="skani", engine_version="0.3.2",
            search_exit_code=command["exit_code"], query_sha256=quality["query_sha256"],
            search_parameters_record=str(output / "skani-search.command.json"),
            search_parameters_sha256=sha256(output / "skani-search.command.json"),
            raw_hits_sha256=sha256(raw), reference_package_sha256=expected_package_sha256,
            reference_component_id="primary_whole", execution_provenance="executed_by_development_workflow",
            quality_evidence_path=str(output / "quality.json"), quality_evidence_sha256=sha256(output / "quality.json"))
        provenance["query_engine_scope"] = query_engine_scope
        evidence = collect_evidence(read_tsv(raw, SKANI_FIELDS), context.evidence_references(),
                                    query_search_path=str(query), provenance=provenance)
        evidence["qc"] = quality
        if mlst_config is not None:
            # Optional module import keeps non-ST execution independent of its
            # native toolchain; the pinned runner owns version and copy checks.
            from .mlst_runner import run_mlst
            mlst_result = run_mlst(query, mlst_config, output / "mlst-run",
                expected_config_sha256=expected_mlst_config_sha256, threads=threads)
            if (mlst_result.get("query_sha256") != quality["query_sha256"]
                    or mlst_result.get("taxonomic_assignment_performed") is not False):
                raise ValueError("MLST result has mismatched query or taxonomic role")
            evidence["mlst"] = mlst_result
            _write(output / "mlst.json", mlst_result)
            for path, digest in pins.items():
                if sha256(Path(path)) != digest:
                    raise ValueError(f"input changed during MLST execution: {path}")
        ep = output / "evidence.json"
        _write(ep, evidence)
        report = context.report(ep, development_guides=development_guides)
        if experimental_rank_rule is not None:
            hierarchy = evaluate_rank(evidence, experimental_rank_rule)
            if hierarchy["assignment"] != report["assignment"]:
                raise ValueError("experimental rank evaluator changed assignment state")
            _write(output / "experimental-rank.json", hierarchy)
            report["experimental_rank_candidate_evidence"] = hierarchy
        if eligibility_config is not None:
            # Built-in screen only. The evidence and QC were produced and
            # pin-checked by this run, so the pure evaluator binds its file digest.
            evidence_digest = sha256(ep)
            eligibility = evaluate_eligibility(evidence, experimental_rank_rule, eligibility_config,
                reference_inventory=primary_inventory(context), evidence_sha256=evidence_digest)
            if eligibility["assignment"] != report["assignment"] or eligibility.get("rank_evidence", hierarchy) != hierarchy:
                raise ValueError("experimental eligibility changed assignment state or rank evidence")
            eligibility["workflow_verification"] = dict(
                evidence_path=str(ep), evidence_sha256=evidence_digest, produced_by="this_development_workflow_run",
                optional_envelopes="not_accepted_by_analyze; bind_saved_evidence_sha256_with_post_analysis_EligibilityContext")
            _write(output / "experimental-eligibility.json", eligibility)
            report["experimental_eligibility"] = eligibility
        if gv_context is not None:
            parents = {c["parent_concept_id"] for c in gv_context.concepts.values()}
            if hierarchy["subspecies"]["candidate_id"] in parents:
                gv_raw = output / "genomovar-hits.tsv"
                gv_command = _run([str(binary), "search", "-q", str(query), "-d", str(gv_database),
                    "-o", str(gv_raw), "--min-af", "0", "-s", "80", "-t", str(threads)], output, "genomovar-search")
                for path, digest in pins.items():
                    if sha256(Path(path)) != digest:
                        raise ValueError(f"input changed during genomovar execution: {path}")
                gv_provenance = dict(input_kind="assembly", engine="skani", engine_version="0.3.2",
                    search_exit_code=gv_command["exit_code"], query_sha256=quality["query_sha256"],
                    reference_component_id="genomovar_chromosome", reference_package_sha256=expected_package_sha256,
                    search_parameters_record=str(output / "genomovar-search.command.json"),
                    search_parameters_sha256=sha256(output / "genomovar-search.command.json"),
                    raw_hits_sha256=sha256(gv_raw), execution_provenance="executed_by_development_workflow")
                gv_result = gv_context.evaluate(read_tsv(gv_raw, SKANI_FIELDS), upstream_evidence=evidence, provenance=gv_provenance)
            else:
                gv_result = dict(status="not_evaluated", reason="no_experimental_candidate_in_target_subspecies",
                    target_parent_concept_ids=sorted(parents), assigned_genomovar_id=None,
                    membership_policy_enabled=False, biological_interpretation="not_a_negative_genomovar_result")
            _write(output / "genomovar-compatibility.json", gv_result)
            report["genomovar_compatibility"] = gv_result
        # Supplementary searches extend the execution window. Recheck original
        # inputs and imported estimates immediately before publishing success.
        for path, digest in pins.items():
            if sha256(Path(path)) != digest:
                raise ValueError(f"input changed before final reporting: {path}")
        if quality_evidence_path is not None:
            if attach_checkm2(structural, quality_evidence_path, expected_sha256=expected_quality_sha256) != quality:
                raise ValueError("quality evidence changed before final reporting")
        write_reports([report], json_path=output / "report.json", tsv_path=output / "report.tsv")
        artifacts = {str(p.relative_to(output)): sha256(p) for p in sorted(output.rglob("*")) if p.is_file()}
        record = dict(schema_version=SCHEMA, status="complete", started_utc=started,
            elapsed_seconds=time.monotonic()-clock, child_peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
            child_peak_rss_scope="cumulative_process_children_linux_kib_not_isolated_per_command",
            assignment_status=evidence["assignment"]["status"], search_estimates=evidence["unique_hit_count"],
            missing_reference_estimates=evidence["not_returned_reference_count"],
            positive_calling_policy_enabled=False, artifacts_sha256=artifacts)
        _write(output / "run.json", record)
        return record
    except (ValueError, OSError, KeyError, TypeError, subprocess.SubprocessError) as exc:
        _write(output / "failure.json", dict(schema_version=SCHEMA, status="failed",
            started_utc=started, elapsed_seconds=time.monotonic()-clock, error=str(exc),
            interpretation="execution_failure_not_negative_taxonomic_evidence"))
        raise


def main(argv: list[str] | None = None, *, prog: str | None = None) -> None:
    p = argparse.ArgumentParser(prog=prog, description=__doc__)
    p.add_argument("--query", required=True, type=Path)
    p.add_argument("--reference-package", required=True, type=Path)
    p.add_argument("--expected-package-sha256", required=True)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--development-guides", type=Path, help="Explicit display annotations; never assignment gates")
    p.add_argument("--quality-evidence", type=Path, help="One pinned local CheckM2 evidence JSON object; imports estimates without applying a QC gate")
    p.add_argument("--expected-quality-sha256", help="Required with --quality-evidence")
    p.add_argument("--experimental-rank-rule", type=Path, help="Explicit disabled development rule; candidate hierarchy only, never positive assignments")
    p.add_argument("--checkm2-config", type=Path, help="Pinned configuration to run CheckM2 within this analysis; no QC gate is implied")
    p.add_argument("--expected-checkm2-config-sha256", help="Required with --checkm2-config")
    p.add_argument("--experimental-genomovars", action="store_true", help="Investigate frozen all-member compatibility when the experimental rank result nominates a represented target subspecies")
    p.add_argument("--mlst-config", type=Path, help="Pinned historical MLST package/runtime configuration; supplementary typing only")
    p.add_argument("--expected-mlst-config-sha256", help="Required with --mlst-config")
    p.add_argument("--experimental-eligibility-config", type=Path, help="Pinned disabled eligibility configuration (frozen candidate policy) bound to --experimental-rank-rule; needs no QC estimate; reports hypothetical eligibility and any near-boundary note separately, never assignments")
    p.add_argument("--expected-eligibility-config-sha256", help="Required with --experimental-eligibility-config")
    p.add_argument("--skani-binary", type=Path, help="Explicit local executable matching the package pin; required for a relocated package")
    a = p.parse_args(argv)
    try:
        result = run_assembly(query=a.query, reference_package=a.reference_package,
            expected_package_sha256=a.expected_package_sha256, output=a.output, threads=a.threads,
            development_guides=json.loads(a.development_guides.read_text()) if a.development_guides else None,
            quality_evidence_path=a.quality_evidence, expected_quality_sha256=a.expected_quality_sha256,
            experimental_rank_rule=json.loads(a.experimental_rank_rule.read_text()) if a.experimental_rank_rule else None,
            checkm2_config=a.checkm2_config, expected_checkm2_config_sha256=a.expected_checkm2_config_sha256,
            experimental_genomovars=a.experimental_genomovars,
            mlst_config=a.mlst_config, expected_mlst_config_sha256=a.expected_mlst_config_sha256,
            experimental_eligibility_config=a.experimental_eligibility_config,
            expected_eligibility_config_sha256=a.expected_eligibility_config_sha256, skani_binary=a.skani_binary)
        print(canonical_json(result), end="")
    except (ValueError, OSError, KeyError, TypeError, subprocess.SubprocessError) as exc:
        p.exit(2, f"development workflow error: {exc}\n")


if __name__ == "__main__":
    main()
