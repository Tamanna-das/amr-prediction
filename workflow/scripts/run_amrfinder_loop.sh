#!/usr/bin/env bash
# Re-run Snakemake until every genome in the list has AMRFinderPlus output.
export PATH="/opt/miniconda3/envs/amr-amrfinder/bin:$PATH"
for round in $(seq 1 30); do
  done_n=$(ls results/amrfinder/*.tsv 2>/dev/null | wc -l | tr -d ' ')
  [ "$done_n" -ge 500 ] && { echo "all 500 done"; exit 0; }
  echo "$(date +%H:%M) round $round: $done_n done, starting snakemake"
  snakemake -s workflow/Snakefile --cores 6 --keep-going --rerun-triggers mtime --quiet
  sleep 60
done
