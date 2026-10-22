#!/usr/bin/env bash
# Tested commands; the batch pipeline is already submitted. Do not resubmit it.
python -m unittest discover -s tests -p test_task14_standard_kraken.py -v
bash -n results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/run-analysis.sh
csvtk nrow -t results/benchmark-design/2026-10-03/task14-astra-whole-assembly-design/candidate-whole-inputs.tsv
csvtk version
ssubmit --version
${USER_HOME}/sw/miniforge3/envs/kraken/bin/kraken2 --version
${USER_HOME}/sw/miniforge3/bin/python --version
# Read-only GitHub retrieval required network access; first sandbox attempt failed.
python scripts/fetch_github_decisions.py --id DEC-008 --id DEC-004 --output results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/github-start-retry
# Exact Slurm plan, validated test and submission argv are in logs/submission*.json.
# Pending post-completion check (not yet executed):
# python scripts/validate_task14_standard_kraken.py --run results/standard-kraken/2026-10-03/task14-sol61-whole-fragments --output results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/independent-validation.json
