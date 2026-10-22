#!/usr/bin/env bash
# These postprocessing commands completed with exit0. Existing outputs retained.
/usr/bin/time -v -o results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/logs/independent-validation.time python results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/source-snapshot/validate_task14_standard_kraken.py --run results/standard-kraken/2026-10-03/task14-sol61-whole-fragments --output results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/independent-validation.json
/usr/bin/time -v -o results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/logs/profile-details.time python results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/source-snapshot/task14_kraken_profile_details.py --run results/standard-kraken/2026-10-03/task14-sol61-whole-fragments
/usr/bin/time -v -o results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/logs/provenance-strata.time python results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/summarize-profile-strata.py
csvtk nrow -t results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/summaries/assembly-accounting.tsv
csvtk nrow -t results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/profile-details/native-contig-calls.tsv
csvtk nrow -t results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/database/assembly-overlap.tsv
python scripts/fetch_github_decisions.py --id DEC-008 --id DEC-004 --output results/standard-kraken/2026-10-03/task14-sol61-whole-fragments/github-handover
sacct -j 29315580 --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList --parsable2
sha256sum ${USER_HOME}/sw/miniforge3/envs/kraken/bin/kraken2 ${USER_HOME}/sw/miniforge3/envs/kraken/share/kraken2-2.1.6/libexec/classify ${USER_HOME}/sw/miniforge3/envs/kraken/share/kraken2-2.1.6/libexec/kraken2lib.pm
