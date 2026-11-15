#!/usr/bin/env bash
# Download BV-BRC assemblies as FASTA via the HTTPS API (FTPS stalled, Nov 2026).
# Checks contig count + length against selected_genomes.tsv; resumable.
set -uo pipefail
TABLE=${1:?usage: download_genomes.sh SELECTED_GENOMES_TSV OUT_DIR}
OUT=${2:?usage: download_genomes.sh SELECTED_GENOMES_TSV OUT_DIR}
API="https://www.bv-brc.org/api/genome_sequence/"
mkdir -p "$OUT"
FREE_GB=$(df -Pk "$OUT" | awk 'NR==2 {print int($4 / 1024 / 1024)}')
[ "$FREE_GB" -lt 10 ] && { echo "Only ${FREE_GB} GB free; need 10. Stopping."; exit 1; }
echo "free disk: ${FREE_GB} GB"
JOBS=$(awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) c[$i]=i; next}
  {print $c["genome_id"] "\t" $c["contigs"] "\t" $c["genome_length"]}' "$TABLE")
total=$(printf '%s\n' "$JOBS" | grep -c .)
n=0; ok=0; skip=0; fail=0
: > "$OUT/failed.txt"
while IFS=$'\t' read -r id wc wl; do
  [ -z "$id" ] && continue
  n=$((n + 1)); dest="$OUT/$id.fna"; part="$dest.part"
  if [ -s "$dest" ]; then skip=$((skip + 1))
  else
    reason=""
    url="${API}?eq(genome_id,${id})&limit(25000)&http_accept=application/dna+fasta"
    if ! curl --fail --silent --show-error --globoff --retry 3 --retry-delay 5 \
         --connect-timeout 30 --max-time 300 -o "$part.raw" "$url"; then
      reason="download error"
    elif [ "$(head -c 1 "$part.raw")" != ">" ]; then
      reason="not FASTA"
    else
      awk '/^>/ {k = split($1, a, "|"); print ">" a[k]; next} {print toupper($0)}' "$part.raw" > "$part"
      gc=$(grep -c '^>' "$part")
      gl=$(awk '!/^>/ {s += length($0)} END {print s + 0}' "$part")
      awk -v gc="$gc" -v wc="$wc" -v gl="$gl" -v wl="$wl" \
        'BEGIN {exit !(gc+0 == wc+0 && gl+0 == wl+0)}' \
        || reason="mismatch: contigs $gc vs $wc, length $gl vs $wl"
    fi
    if [ -z "$reason" ]; then mv "$part" "$dest"; ok=$((ok + 1))
    else printf '%s\t%s\n' "$id" "$reason" >> "$OUT/failed.txt"; fail=$((fail + 1)); fi
    rm -f "$part" "$part.raw"; sleep 0.5
  fi
  [ $((n % 25)) -eq 0 ] && echo "[$n/$total] downloaded=$ok already_present=$skip failed=$fail"
done <<< "$JOBS"
echo "done: $n genomes, downloaded=$ok, already present=$skip, failed=$fail"
