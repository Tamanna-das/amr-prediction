# Week 5 decisions: data scope, labels and genome selection

All counts below come from files downloaded on 2026-10-31. Re-running the
pipeline on a later BV-BRC release will change them.

## Data sources

| File | Source | Notes |
|---|---|---|
| `PATRIC_genomes_AMR.txt` | `ftp://ftp.bv-brc.org/RELEASE_NOTES/` (FTPS) | 1,207,585 rows. The `ftp.bvbrc.org` hostname in the 2021 paper and plain `ftp://` no longer work. |
| `genome_metadata` | same directory | File dated 2024-09-13; still covered all 4,531 labelled genomes. |
| Assemblies | `ftp://ftp.bv-brc.org/genomes/<id>/<id>.fna` (FTPS) | |
| AMRFinderPlus | v4.2.7, database 2026-08-07.1 | Organism option `Klebsiella_pneumoniae`. |

## Species

Taxon IDs `573` (*K. pneumoniae*) and `72407` (subsp. *pneumoniae*), selected by
ID rather than name. 324 genomes named *K. pneumoniae* but filed under
strain-level taxon IDs were excluded.

## Labels

- **Laboratory phenotypes only.** `antimicrobial_resistance_evidence` is
  "AMR Panel" for every selected genome; no computationally predicted labels.
- **Provided calls only.** 38,179 *K. pneumoniae* rows have an MIC but no
  Resistant/Susceptible call. Re-calling them would require breakpoints matched
  to each testing standard and year, so this is left as future work.
- **Intermediate dropped**, along with "nonsusceptible". Only about 5 % of
  genomes per drug; dropping keeps the task cleanly binary.
- **Conflicting calls:** a genome-drug pair with more than one distinct call is
  dropped (58 pairs). The genome's other drugs are kept.
- **Known label noise:** calls come from mixed standards (CLSI, EUCAST, BSAC)
  and years.

## Drugs

| Drug | R | S | Reason |
|---|---|---|---|
| meropenem | 1,442 | 2,654 | key clinical drug; expected acquired-gene mechanism |
| ciprofloxacin | 2,613 | 1,204 | expected chromosomal point mutations (tests mutation calls) |
| gentamicin | 1,721 | 2,606 | different drug class, well balanced |
| trimethoprim/sulfamethoxazole | 2,416 | 1,153 | combination drug, two determinants |

Counts are usable labels after the filters above. Expected mechanisms are to be
verified against CARD and the NCBI Reference Gene Catalog in Week 7.
Not chosen: colistin (338 R; gene-disruption mechanisms), imipenem (redundant
with meropenem), "extended spectrum beta lactamase" (a phenotype test, not a drug).

## Assembly QC

Thresholds come from the distribution of the 3,045 genomes labelled for all four
drugs (1st to 99th percentiles: length 5.11 to 6.11 Mb, GC 56.3 to 57.7 %,
contigs up to 689).

| Rule | Keep |
|---|---|
| genome_status | not "Deprecated" |
| genome_length | 5.0 to 6.3 Mb |
| gc_content | 56.0 to 58.0 % |
| contigs | 500 or fewer |

Genomes failing each rule: status 1, length 39, GC 27, contigs 50 (72 fail any rule; 2,973 pass).

## Leakage controls at this stage

- **Same isolate:** 8 BioSamples had 2 genomes each; the assembly with fewer
  contigs is kept. Near-identical genomes without a BioSample are handled with
  Mash in Week 6.
- **Study concentration:** BioProject PRJNA376414 holds 1,370 of the 3,045
  candidate genomes (45 %). The sample is capped at 360 genomes per BioProject.
  What this project is: [look it up on NCBI BioProject and describe it here].

## Sample

500 genomes, random with seed 42, at most 150 per BioProject, from the 2,966 genomes passing QC and dedup (reduced from a planned 1,200 because of AMRFinderPlus run time).

| Drug | R | S |
|---|---|---|
| meropenem | 155 | 345 |
| ciprofloxacin | 319 | 181 |
| gentamicin | 182 | 318 |
| trimethoprim/sulfamethoxazole | 301 | 199 |
