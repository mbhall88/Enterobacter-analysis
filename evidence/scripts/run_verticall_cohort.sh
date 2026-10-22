#!/usr/bin/env bash
#
# Verticall / FastME on any frozen cohort, one stage per invocation.
#
# Task 08's wrapper hardcoded the Task 06 pilot: its manifest, its derivative map and a single
# in-process pairwise run. Task 11 needs the same steps over a different cohort, split across many
# independently re-runnable Slurm jobs, so this script takes --manifest, --derivative-map,
# --run-dir, --threads and --part, and the word "pilot" appears nowhere in it.
#
# Why the stages are separate invocations
# ---------------------------------------
# Verticall 0.4.2 has **no resume**. Partitioning is the substitute: one part that dies is one part
# to re-run. That only works if the parts share a prepared input directory and never race to write
# it, so staging and indexing are their own stages and run to completion first:
#
#   stage     copy one derivative per tip into a private directory, verifying every copy against
#             the hash the freeze recorded. Real copies, never symlinks: verticall resolves the
#             assembly path, so a symlinked stage would point the run back at protected input.
#   index     'verticall pairwise --index_only' once. Without this each part would rebuild a
#             minimap2 index for every genome -- the single thing that makes 100 parts cheap.
#   pairwise  'verticall pairwise --part n/m --skip_check' for one part, into its own TSV.
#   merge     concatenate the parts with an explicit header check, and count every row.
#   downstream  QC gate, the four matrices, and FastME on each.
#
# The merge is the dangerous step, so it is checked rather than trusted. In this build only part 1
# writes a header (verticall/pairwise.py:38 'only include the header in the first part'), and
# 'verticall matrix' reads line 0 as the header (matrix.py:62). A merge that lost part 1 would
# therefore start with a data row; a merge that duplicated a header would carry one mid-file.
# Both abort 'verticall matrix' loudly, but this script refuses them earlier and says which.
#
# Usage:
#   run_verticall_cohort.sh --stage stage     --manifest IDS --derivative-map MAP --run-dir DIR
#   run_verticall_cohort.sh --stage index     --run-dir DIR [--threads N]
#   run_verticall_cohort.sh --stage pairwise  --run-dir DIR --part i/m [--threads N]
#   run_verticall_cohort.sh --stage merge     --run-dir DIR --parts M
#   run_verticall_cohort.sh --stage downstream --manifest IDS --run-dir DIR [--threads N]
#
set -euo pipefail

ROOT=${PROJECT_ROOT}
VERTICALL_ENV="$ROOT/.work/envs/verticall"

STAGE_NAME=""
MANIFEST=""
DERIVATIVE_MAP=""
RUN_DIR=""
THREADS=8
PART="1/1"
PARTS=1
LABEL="cohort"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stage)           STAGE_NAME="$2"; shift 2 ;;
    --manifest)        MANIFEST="$2"; shift 2 ;;
    --derivative-map)  DERIVATIVE_MAP="$2"; shift 2 ;;
    --run-dir)         RUN_DIR="$2"; shift 2 ;;
    --threads)         THREADS="$2"; shift 2 ;;
    --part)            PART="$2"; shift 2 ;;
    --parts)           PARTS="$2"; shift 2 ;;
    --label)           LABEL="$2"; shift 2 ;;
    *) echo "Error: unknown argument $1" >&2; exit 2 ;;
  esac
done

[[ -n "$RUN_DIR" ]] || { echo "Error: --run-dir is required" >&2; exit 2; }
case "$STAGE_NAME" in
  stage|index|pairwise|merge|downstream) ;;
  *) echo "Error: --stage must be stage, index, pairwise, merge or downstream" >&2; exit 2 ;;
esac

export PATH="$VERTICALL_ENV/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MPLBACKEND=Agg
cd "$ROOT"

STAGE_DIR="$RUN_DIR/staged/$LABEL"
PARTS_DIR="$RUN_DIR/parts"
LOGS="$RUN_DIR/logs"
TABLES="$RUN_DIR/tables"
LEDGER="$TABLES/stage-timings.tsv"
MERGED="$RUN_DIR/${LABEL}-pairwise.tsv"
mkdir -p "$LOGS" "$TABLES" "$PARTS_DIR" "$RUN_DIR/matrices" "$RUN_DIR/trees" "$RUN_DIR/qc"

if [[ ! -f "$LEDGER" ]]; then
  printf 'stage\tlabel\tpart\tstarted_utc\tfinished_utc\telapsed_seconds\texit_status\tcommand\tstdout_log\tstderr_log\n' > "$LEDGER"
fi

# run_stage NAME -- COMMAND...
# Runs under /usr/bin/time -v, captures stdout/stderr, appends one ledger row with the exact
# command, timing and exit status, and stops on a real failure rather than masking it.
run_stage() {
  local name="$1"; shift
  if [[ "${1:-}" == "--" ]]; then shift; fi
  local out="$LOGS/${LABEL}-${name}.stdout"
  local err="$LOGS/${LABEL}-${name}.stderr"
  local timing="$LOGS/${LABEL}-${name}.time"
  local started finished status start_s end_s
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; start_s="$(date +%s)"
  set +e
  /usr/bin/time -v -o "$timing" "$@" > "$out" 2> "$err"
  status=$?
  set -e
  end_s="$(date +%s)"; finished="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$name" "$LABEL" "$PART" "$started" "$finished" "$((end_s - start_s))" "$status" "$*" \
    "$out" "$err" >> "$LEDGER"
  if [[ $status -ne 0 ]]; then
    echo "Error: stage '$name' exited $status; see $err" >&2
    exit "$status"
  fi
}

# --------------------------------------------------------------------------------------- stage
do_stage() {
  [[ -n "$MANIFEST" ]]       || { echo "Error: --manifest is required for --stage stage" >&2; exit 2; }
  [[ -n "$DERIVATIVE_MAP" ]] || { echo "Error: --derivative-map is required" >&2; exit 2; }
  [[ -f "$MANIFEST" ]]       || { echo "Error: missing $MANIFEST" >&2; exit 2; }
  [[ -f "$DERIVATIVE_MAP" ]] || { echo "Error: missing $DERIVATIVE_MAP" >&2; exit 2; }
  mkdir -p "$STAGE_DIR"

  local checksums="$TABLES/${LABEL}-staged.sha256"
  local manifest_out="$TABLES/${LABEL}-staged-manifest.tsv"
  local started start_s copied=0 reused=0
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; start_s="$(date +%s)"
  : > "$checksums"
  printf 'tip_id\tsource_derivative_path\tsource_derivative_sha256\tstaged_path\tstaged_sha256\taction\n' > "$manifest_out"

  # Column positions are read from the map's own header, never assumed.
  local col_tip col_path col_sha
  col_tip=$(head -1 "$DERIVATIVE_MAP" | tr '\t' '\n' | grep -n '^tip_id$' | cut -d: -f1)
  col_path=$(head -1 "$DERIVATIVE_MAP" | tr '\t' '\n' | grep -n '^derivative_fasta_path$' | cut -d: -f1)
  col_sha=$(head -1 "$DERIVATIVE_MAP" | tr '\t' '\n' | grep -n '^derivative_sha256$' | cut -d: -f1)
  [[ -n "$col_tip" && -n "$col_path" && -n "$col_sha" ]] || {
    echo "Error: $DERIVATIVE_MAP lacks tip_id/derivative_fasta_path/derivative_sha256" >&2; exit 2; }

  while IFS= read -r tip; do
    [[ -n "$tip" ]] || continue
    local row src sha dest staged_sha action
    row="$(awk -F'\t' -v t="$tip" -v c="$col_tip" '$c == t {print; exit}' "$DERIVATIVE_MAP")"
    [[ -n "$row" ]] || { echo "Error: $tip is not in $DERIVATIVE_MAP" >&2; exit 2; }
    src="$(printf '%s' "$row" | cut -f"$col_path")"
    sha="$(printf '%s' "$row" | cut -f"$col_sha")"
    [[ "$src" = /* ]] || src="$ROOT/$src"
    dest="$STAGE_DIR/${tip}.fna"

    if [[ -f "$dest" ]] && [[ "$(sha256sum "$dest" | cut -d' ' -f1)" == "$sha" ]]; then
      staged_sha="$sha"; action=reused; reused=$((reused + 1))
    else
      # A stale copy would silently pair with a stale .mmi of the same sample name.
      cp -f --no-preserve=mode "$src" "$dest"
      rm -f "$STAGE_DIR/${tip}.mmi"
      staged_sha="$(sha256sum "$dest" | cut -d' ' -f1)"
      action=copied; copied=$((copied + 1))
    fi
    [[ "$staged_sha" == "$sha" ]] || {
      echo "Error: staged $tip hashes $staged_sha, the freeze recorded $sha" >&2; exit 1; }
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$tip" "$src" "$sha" "$dest" "$staged_sha" "$action" >> "$manifest_out"
    printf '%s  %s\n' "$sha" "$dest" >> "$checksums"
  done < "$MANIFEST"

  sha256sum -c "$checksums" > "$LOGS/${LABEL}-stage-verify.stdout" 2> "$LOGS/${LABEL}-stage-verify.stderr"
  local expected actual
  expected="$(grep -c . "$MANIFEST")"
  actual="$(find "$STAGE_DIR" -maxdepth 1 -name '*.fna' -type f | wc -l)"
  [[ "$expected" == "$actual" ]] || {
    echo "Error: staged $actual FASTAs, expected $expected" >&2; exit 1; }
  [[ -z "$(find "$STAGE_DIR" -maxdepth 1 -type l)" ]] || {
    echo "Error: $STAGE_DIR contains a symlink" >&2; exit 1; }

  printf '%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\t%s\t%s\n' \
    "stage" "$LABEL" "-" "$started" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(( $(date +%s) - start_s ))" \
    "stage $MANIFEST -> $STAGE_DIR (copied=$copied reused=$reused)" \
    "$LOGS/${LABEL}-stage-verify.stdout" "$LOGS/${LABEL}-stage-verify.stderr" >> "$LEDGER"
  echo "staged $actual FASTAs into $STAGE_DIR (copied=$copied reused=$reused)" >&2
}

# --------------------------------------------------------------------------------------- index
do_index() {
  [[ -d "$STAGE_DIR" ]] || { echo "Error: $STAGE_DIR does not exist; run --stage stage" >&2; exit 2; }
  run_stage index -- verticall pairwise -i "$STAGE_DIR" -o "$RUN_DIR/index-only.tsv" \
    -t "$THREADS" --index_only
  find "$STAGE_DIR" -maxdepth 1 -type f -name '*.mmi' -print0 | sort -z | xargs -0 sha256sum \
    > "$TABLES/${LABEL}-indexes.sha256"
  local built expected
  built="$(grep -c . "$TABLES/${LABEL}-indexes.sha256")"
  expected="$(find "$STAGE_DIR" -maxdepth 1 -name '*.fna' -type f | wc -l)"
  [[ "$built" == "$expected" ]] || {
    echo "Error: $built indexes for $expected assemblies" >&2; exit 1; }
  echo "built $built minimap2 indexes in $STAGE_DIR" >&2
}

# ------------------------------------------------------------------------------------ pairwise
do_pairwise() {
  [[ -d "$STAGE_DIR" ]] || { echo "Error: $STAGE_DIR does not exist" >&2; exit 2; }
  local n m out
  n="${PART%%/*}"; m="${PART##*/}"
  out="$PARTS_DIR/$(printf 'part_%s_of_%s.tsv' "$n" "$m")"
  # Never resume into a half-written part: this build has no resume, so a partial file is a wrong
  # answer, not a starting point.
  rm -f "$out"
  run_stage "pairwise-part-${n}-of-${m}" -- verticall pairwise \
    -i "$STAGE_DIR" -o "$out" -t "$THREADS" --part "$PART" --skip_check
  printf '%s\t%s\t%s\t%s\n' "$n" "$m" "$(wc -l < "$out")" "$(sha256sum "$out" | cut -d' ' -f1)" \
    >> "$TABLES/${LABEL}-part-rows.tsv"
  echo "part $PART wrote $(wc -l < "$out") lines to $out" >&2
}

# --------------------------------------------------------------------------------------- merge
do_merge() {
  local started start_s expected_header seen_header=0 data_rows=0 missing=()
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"; start_s="$(date +%s)"
  local ledger="$TABLES/${LABEL}-merge-ledger.tsv"
  printf 'part\tfile\tlines\theader_lines\tdata_rows\tsha256\n' > "$ledger"

  for ((i = 1; i <= PARTS; i++)); do
    local f="$PARTS_DIR/$(printf 'part_%s_of_%s.tsv' "$i" "$PARTS")"
    [[ -f "$f" ]] || { missing+=("$i"); continue; }
  done
  if [[ ${#missing[@]} -gt 0 ]]; then
    echo "Error: ${#missing[@]} of $PARTS parts are missing: ${missing[*]}" >&2
    echo "Refusing to merge an incomplete set. Re-run those parts." >&2
    exit 1
  fi

  : > "$MERGED"
  for ((i = 1; i <= PARTS; i++)); do
    local f="$PARTS_DIR/$(printf 'part_%s_of_%s.tsv' "$i" "$PARTS")"
    local lines head_lines rows
    lines="$(wc -l < "$f")"
    head_lines="$(awk -F'\t' '$1 == "assembly_a" && $2 == "assembly_b"' "$f" | wc -l)"
    rows=$((lines - head_lines))
    if [[ $i -eq 1 ]]; then
      [[ "$head_lines" == "1" ]] || {
        echo "Error: part 1 carries $head_lines header lines, expected exactly 1" >&2; exit 1; }
      [[ "$(head -1 "$f" | cut -f1,2)" == "$(printf 'assembly_a\tassembly_b')" ]] || {
        echo "Error: part 1 does not start with the header line" >&2; exit 1; }
      expected_header="$(head -1 "$f")"
      seen_header=1
    else
      [[ "$head_lines" == "0" ]] || {
        echo "Error: part $i carries $head_lines header lines; only part 1 may carry one" >&2
        exit 1; }
    fi
    cat "$f" >> "$MERGED"
    data_rows=$((data_rows + rows))
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$i" "$f" "$lines" "$head_lines" "$rows" \
      "$(sha256sum "$f" | cut -d' ' -f1)" >> "$ledger"
  done

  [[ $seen_header -eq 1 ]] || { echo "Error: no header was found in any part" >&2; exit 1; }
  [[ "$(head -1 "$MERGED")" == "$expected_header" ]] || {
    echo "Error: the merged file does not begin with the header" >&2; exit 1; }
  local merged_headers
  merged_headers="$(awk -F'\t' '$1 == "assembly_a" && $2 == "assembly_b"' "$MERGED" | wc -l)"
  [[ "$merged_headers" == "1" ]] || {
    echo "Error: the merged file carries $merged_headers header lines, expected 1" >&2; exit 1; }
  local merged_rows
  merged_rows=$(( $(wc -l < "$MERGED") - 1 ))
  [[ "$merged_rows" == "$data_rows" ]] || {
    echo "Error: merged $merged_rows data rows, parts hold $data_rows" >&2; exit 1; }

  sha256sum "$MERGED" > "$TABLES/${LABEL}-pairwise.sha256"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t0\t%s\t%s\t%s\n' \
    "merge" "$LABEL" "1..$PARTS" "$started" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    "$(( $(date +%s) - start_s ))" \
    "merge $PARTS parts -> $MERGED ($merged_rows data rows, 1 header)" "$ledger" "$ledger" >> "$LEDGER"
  echo "merged $PARTS parts into $MERGED: $merged_rows data rows, exactly 1 header line" >&2
}

# ---------------------------------------------------------------------------------- downstream
do_downstream() {
  [[ -n "$MANIFEST" ]] || { echo "Error: --manifest is required for --stage downstream" >&2; exit 2; }
  [[ -f "$MERGED" ]]   || { echo "Error: missing $MERGED; run --stage merge" >&2; exit 2; }

  for spec in "median-first:--multi first" "median-low:--multi low" "median-high:--multi high" \
              "asymmetrical:--multi first --asymmetrical"; do
    local label="${spec%%:*}" opts="${spec#*:}"
    # shellcheck disable=SC2086
    run_stage "matrix-${label}" -- verticall matrix -i "$MERGED" \
      -o "$RUN_DIR/matrices/${LABEL}-${label}.phylip" \
      --distance_type median_vertical_window $opts
  done

  # QC gates the trees: it exits non-zero on an incomplete pairwise run or an incomplete matrix,
  # and set -e stops here rather than building a tree from a matrix with an empty cell in it.
  run_stage qc -- python3 "$ROOT/scripts/qc_verticall.py" \
    --pairwise "$MERGED" --expected-ids "$MANIFEST" --output "$RUN_DIR/qc/$LABEL" \
    --matrix "primary=$RUN_DIR/matrices/${LABEL}-median-first.phylip" \
    --matrix "low=$RUN_DIR/matrices/${LABEL}-median-low.phylip" \
    --matrix "high=$RUN_DIR/matrices/${LABEL}-median-high.phylip" \
    --matrix "asymmetrical=$RUN_DIR/matrices/${LABEL}-asymmetrical.phylip" \
    --expect-asymmetric asymmetrical

  for label in median-first median-low median-high asymmetrical; do
    run_stage "fastme-${label}" -- fastme \
      -i "$RUN_DIR/matrices/${LABEL}-${label}.phylip" \
      -o "$RUN_DIR/trees/${LABEL}-${label}.nwk" \
      -I "$RUN_DIR/trees/${LABEL}-${label}.info" \
      -m B -n B -s -T "$THREADS"
  done

  find "$RUN_DIR/matrices" "$RUN_DIR/trees" -type f -print0 | sort -z | xargs -0 sha256sum \
    > "$RUN_DIR/${LABEL}-outputs.sha256"
}

run_stage versions -- "$VERTICALL_ENV/bin/verticall" --version
case "$STAGE_NAME" in
  stage)      do_stage ;;
  index)      do_index ;;
  pairwise)   do_pairwise ;;
  merge)      do_merge ;;
  downstream) do_downstream ;;
esac
echo "done: $STAGE_NAME ($LABEL${PART:+ part $PART})" >&2
