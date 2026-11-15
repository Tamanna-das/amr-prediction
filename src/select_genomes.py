"""Load, summarise, label and QC genomes from the BV-BRC AMR collection.

Subcommands
-----------
summary   : per-species exploration report from the AMR table
labels    : apply the Week 5 labelling rules and write a long label table
qc-report : join labelled genomes to genome_metadata and report QC distributions
select    : apply QC, BioSample dedup and capped sampling; write the final genome list
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd

# Header as downloaded from ftp.bv-brc.org/RELEASE_NOTES on 2026-10-31.
AMR_COLUMNS = [
    "genome_id", "genome_name", "taxon_id", "antibiotic", "resistant_phenotype",
    "measurement", "measurement_sign", "measurement_value", "measurement_unit",
    "laboratory_typing_method", "laboratory_typing_method_version",
    "laboratory_typing_platform", "vendor", "testing_standard",
    "testing_standard_year", "source",
]

# Subset of genome_metadata columns (header checked 2026-10-31; file dated 2024-09-13).
META_COLUMNS = [
    "genome_id", "genome_name", "taxon_id", "genome_status",
    "biosample_accession", "bioproject_accession", "assembly_accession",
    "sequencing_platform", "contigs", "genome_length", "gc_content", "patric_cds",
    "collection_date", "isolation_country", "antimicrobial_resistance_evidence",
]
NUMERIC_META = ["contigs", "genome_length", "gc_content", "patric_cds"]

# Week 5 decisions (see README / docs for justification).
TAXON_IDS = ["573", "72407"]  # K. pneumoniae; K. pneumoniae subsp. pneumoniae
DRUGS = [
    "meropenem",
    "ciprofloxacin",
    "gentamicin",
    "trimethoprim/sulfamethoxazole",
]
LABEL_MAP = {"resistant": 1, "susceptible": 0}  # intermediate etc. are dropped

# QC thresholds, chosen from the observed distribution of the 3,045 complete-label
# genomes (1st-99th percentiles: length 5.11-6.11 Mb, GC 56.3-57.7 %, contigs <= 689).
QC = {
    "min_length": 5_000_000,
    "max_length": 6_300_000,
    "min_gc": 56.0,
    "max_gc": 58.0,
    "max_contigs": 500,
    "bad_status": ("Deprecated",),
}
SAMPLE_SIZE = 1200
MAX_PER_BIOPROJECT = 360  # 30 % of the sample; one BioProject holds 45 % of the pool
SEED = 42


# ---------------------------------------------------------------- loading
def load_amr_table(path: str | Path) -> pd.DataFrame:
    """Read the AMR table with every column as text and light normalisation.

    - dtype=str keeps genome IDs like '573.12340' intact (no float conversion).
    - QUOTE_NONE stops pandas guessing at the odd quoting in genome_name;
      we strip the quotes ourselves.
    - Empty fields become missing values (NaN).
    """
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        quoting=csv.QUOTE_NONE,
        keep_default_na=False,
        na_values=[""],
    )
    missing = [c for c in AMR_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"AMR table is missing columns: {missing}")
    df["genome_name"] = df["genome_name"].str.strip('"').str.strip()
    for col in ("antibiotic", "resistant_phenotype"):
        df[col] = df[col].str.strip().str.lower()
    return df


def load_genome_metadata(
    path: str | Path, genome_ids: set[str], chunksize: int = 200_000
) -> pd.DataFrame:
    """Read only META_COLUMNS for the given genome IDs, in chunks (low memory).

    Default CSV quoting is used here (unlike the AMR table) because
    genome_metadata has free-text columns that may contain quoted tabs.
    Malformed lines are skipped with a warning rather than crashing.
    """
    parts = []
    reader = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        usecols=META_COLUMNS,
        keep_default_na=False,
        na_values=[""],
        chunksize=chunksize,
        on_bad_lines="warn",
    )
    for chunk in reader:
        parts.append(chunk.loc[chunk["genome_id"].isin(genome_ids)])
    meta = pd.concat(parts, ignore_index=True)
    meta["genome_name"] = meta["genome_name"].str.strip('"').str.strip()
    for col in NUMERIC_META:
        meta[col] = pd.to_numeric(meta[col], errors="coerce")
    return meta


# ---------------------------------------------------------------- filtering
def filter_species(df: pd.DataFrame, species: str) -> pd.DataFrame:
    """Rows whose genome_name starts with the species name (exploration only)."""
    mask = df["genome_name"].str.startswith(species, na=False)
    return df.loc[mask].copy()


def filter_taxa(df: pd.DataFrame, taxon_ids: list[str]) -> pd.DataFrame:
    """Rows whose taxon_id is in an explicit allow-list (used for real selection)."""
    return df.loc[df["taxon_id"].isin(taxon_ids)].copy()


# ---------------------------------------------------------------- summaries
def phenotype_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Unique genomes per antibiotic x phenotype (duplicate rows counted once)."""
    labelled = df.dropna(subset=["resistant_phenotype"])
    per_genome = labelled.drop_duplicates(
        ["genome_id", "antibiotic", "resistant_phenotype"]
    )
    return pd.crosstab(per_genome["antibiotic"], per_genome["resistant_phenotype"])


def conflicting_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """Genome-antibiotic pairs that carry more than one distinct phenotype."""
    labelled = df.dropna(subset=["resistant_phenotype"])
    n = labelled.groupby(["genome_id", "antibiotic"])["resistant_phenotype"].nunique()
    return n[n > 1].reset_index()[["genome_id", "antibiotic"]]


def print_summary(df: pd.DataFrame, species: str, top: int = 25) -> None:
    sub = filter_species(df, species)
    print(f"species filter: {species!r}")
    print(f"rows: {len(sub):,}   unique genomes: {sub['genome_id'].nunique():,}")

    no_call = sub["resistant_phenotype"].isna() & sub["measurement"].notna()
    print(f"rows with a measurement but no phenotype call: {no_call.sum():,}")

    for col in ("taxon_id", "resistant_phenotype", "laboratory_typing_method",
                "testing_standard", "source"):
        print(f"\n# {col}")
        print(sub[col].value_counts(dropna=False).head(15).to_string())

    counts = phenotype_counts(sub)
    if {"resistant", "susceptible"} <= set(counts.columns):
        counts["min_R_S"] = counts[["resistant", "susceptible"]].min(axis=1)
        counts = counts.sort_values("min_R_S", ascending=False)
    print(f"\n# unique genomes per antibiotic x phenotype (top {top})")
    print(counts.head(top).to_string())

    conflicts = conflicting_pairs(sub)
    print(f"\n# genome-antibiotic pairs with conflicting phenotypes: {len(conflicts):,}")


# ---------------------------------------------------------------- labelling
def build_labels(df: pd.DataFrame, drugs: list[str]) -> pd.DataFrame:
    """Long table (genome_id, antibiotic, label) with 1 = resistant, 0 = susceptible.

    Rules, applied per genome-drug pair:
      1. keep only the requested drugs and rows that carry a phenotype call;
      2. drop the pair if it has more than one distinct call (conflict),
         counting intermediate as a distinct call;
      3. keep only resistant / susceptible (intermediate etc. dropped).
    """
    called = df.loc[df["antibiotic"].isin(drugs)].dropna(subset=["resistant_phenotype"])
    called = called.drop_duplicates(["genome_id", "antibiotic", "resistant_phenotype"])
    n_calls = called.groupby(["genome_id", "antibiotic"])["resistant_phenotype"].transform(
        "nunique"
    )
    consistent = called.loc[n_calls == 1]
    kept = consistent.loc[consistent["resistant_phenotype"].isin(list(LABEL_MAP))]
    labels = kept[["genome_id", "antibiotic"]].copy()
    labels["label"] = kept["resistant_phenotype"].map(LABEL_MAP).astype("int8")
    return labels.sort_values(["antibiotic", "genome_id"]).reset_index(drop=True)


def labels_wide(labels: pd.DataFrame) -> pd.DataFrame:
    """Genome x drug matrix of labels; missing = no usable label for that drug."""
    return labels.pivot(index="genome_id", columns="antibiotic", values="label").astype(
        "Int8"
    )


def load_labels(path: str | Path) -> pd.DataFrame:
    """Read the label table written by the `labels` subcommand."""
    return pd.read_csv(path, sep="\t", dtype={"genome_id": str, "antibiotic": str,
                                              "label": "int8"})


def print_label_report(df: pd.DataFrame, labels: pd.DataFrame, drugs: list[str]) -> None:
    kp = filter_species(df, "Klebsiella pneumoniae")
    excluded = kp.loc[~kp["taxon_id"].isin(TAXON_IDS)].drop_duplicates("genome_id")
    print(f"genomes named K. pneumoniae but outside taxa {TAXON_IDS}: {len(excluded):,}")
    print(excluded[["taxon_id", "genome_name"]].head(15).to_string(index=False))

    counts = labels.groupby("antibiotic")["label"].agg(n="size", resistant="sum")
    counts["susceptible"] = counts["n"] - counts["resistant"]
    print("\n# usable labels per drug (after taxa, conflict and intermediate filters)")
    print(counts.to_string())

    wide = labels_wide(labels)
    coverage = wide.notna().sum(axis=1).value_counts().sort_index()
    print(f"\n# genomes by number of drugs labelled (out of {len(drugs)})")
    print(coverage.to_string())

    complete = wide.dropna()
    print(f"\ngenomes labelled for all {len(drugs)} drugs: {len(complete):,}")
    if len(complete):
        print("R / S among those genomes:")
        for drug in drugs:
            if drug in complete.columns:
                r = int(complete[drug].sum())
                print(f"  {drug}: {r} R / {len(complete) - r} S")


# ---------------------------------------------------------------- QC
def biosample_duplicates(meta: pd.DataFrame) -> pd.DataFrame:
    """BioSamples linked to more than one genome_id (same isolate, multiple assemblies)."""
    with_bs = meta.dropna(subset=["biosample_accession"])
    n = with_bs.groupby("biosample_accession")["genome_id"].nunique()
    dup_bs = n[n > 1].index
    return (with_bs.loc[with_bs["biosample_accession"].isin(dup_bs),
                        ["biosample_accession", "genome_id"]]
            .sort_values(["biosample_accession", "genome_id"])
            .reset_index(drop=True))


def print_qc_report(meta: pd.DataFrame, labels: pd.DataFrame, drugs: list[str]) -> None:
    wide = labels_wide(labels)
    complete_ids = set(wide.dropna().index)
    all_ids = set(wide.index)

    found = set(meta["genome_id"])
    print(f"labelled genomes: {len(all_ids):,}   found in genome_metadata: "
          f"{len(found & all_ids):,}")
    print(f"labelled for all {len(drugs)} drugs: {len(complete_ids):,}   found: "
          f"{len(found & complete_ids):,}")
    print(f"duplicate genome_id rows in metadata: {meta['genome_id'].duplicated().sum():,}")

    sub = meta.loc[meta["genome_id"].isin(complete_ids)]
    print(f"\n=== QC on the {len(sub):,} complete-label genomes found in metadata ===")

    for col in ("genome_status", "antimicrobial_resistance_evidence",
                "sequencing_platform", "taxon_id"):
        print(f"\n# {col}")
        print(sub[col].value_counts(dropna=False).head(10).to_string())

    print("\n# numeric QC fields")
    print(sub[NUMERIC_META].describe(percentiles=[0.01, 0.05, 0.5, 0.95, 0.99])
          .to_string())

    print(f"\n# missing biosample_accession: {sub['biosample_accession'].isna().sum():,}")
    dups = biosample_duplicates(sub)
    print(f"# biosamples with >1 genome: {dups['biosample_accession'].nunique():,} "
          f"({len(dups):,} genomes involved)")

    print("\n# top bioprojects (study / outbreak concentration)")
    print(sub["bioproject_accession"].value_counts(dropna=False).head(10).to_string())


# ---------------------------------------------------------------- selection
def read_kpn_metadata(path: str | Path) -> pd.DataFrame:
    """Read the small per-genome metadata table written by `qc-report`."""
    meta = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    for col in NUMERIC_META:
        meta[col] = pd.to_numeric(meta[col], errors="coerce")
    return meta


def qc_flags(meta: pd.DataFrame, qc: dict = QC) -> pd.DataFrame:
    """One boolean column per QC rule; True means the genome FAILS that rule.

    Missing values fail the numeric rules (we cannot vouch for an assembly
    whose basic statistics are unknown).
    """
    flags = pd.DataFrame(index=meta.index)
    flags["status"] = meta["genome_status"].isin(qc["bad_status"])
    flags["length"] = ~meta["genome_length"].between(qc["min_length"], qc["max_length"])
    flags["gc"] = ~meta["gc_content"].between(qc["min_gc"], qc["max_gc"])
    flags["contigs"] = meta["contigs"].isna() | (meta["contigs"] > qc["max_contigs"])
    return flags


def dedup_biosample(meta: pd.DataFrame) -> pd.DataFrame:
    """Keep one genome per BioSample (fewest contigs, then lowest genome_id).

    Genomes with no BioSample are all kept; near-identical genomes among them
    are handled later with Mash distances.
    """
    has_bs = meta["biosample_accession"].notna()
    ranked = meta.loc[has_bs].sort_values(["contigs", "genome_id"])
    kept = ranked.drop_duplicates("biosample_accession", keep="first")
    out = pd.concat([kept, meta.loc[~has_bs]])
    return out.sort_values("genome_id").reset_index(drop=True)


def sample_genomes(
    meta: pd.DataFrame, n: int, max_per_project: int, seed: int
) -> pd.DataFrame:
    """Random sample of n genomes with at most max_per_project per BioProject.

    Genomes without a BioProject each count as their own project.
    Raises if the cap leaves fewer than n genomes.
    """
    shuffled = meta.sample(frac=1, random_state=seed)
    project = shuffled["bioproject_accession"].fillna("none:" + shuffled["genome_id"])
    capped = shuffled.loc[shuffled.groupby(project).cumcount() < max_per_project]
    if len(capped) < n:
        raise ValueError(f"only {len(capped)} genomes available after capping; "
                         f"lower n or raise max_per_project")
    return capped.head(n).sort_values("genome_id").reset_index(drop=True)


def print_selection_report(
    pool: pd.DataFrame, flags: pd.DataFrame, deduped: pd.DataFrame,
    selected: pd.DataFrame, drugs: list[str],
) -> None:
    print(f"complete-label genomes in pool: {len(pool):,}")
    print("\n# QC failures per rule (a genome can fail several)")
    print(flags.sum().to_string())
    print(f"failing any rule: {flags.any(axis=1).sum():,}   passing: "
          f"{(~flags.any(axis=1)).sum():,}")
    print(f"after BioSample dedup: {len(deduped):,}")
    print(f"selected: {len(selected):,}")

    print("\n# top bioprojects in the selection")
    print(selected["bioproject_accession"].value_counts(dropna=False).head(10).to_string())

    print("\n# class balance per drug in the selection")
    for drug in drugs:
        r = int(selected[drug].sum())
        print(f"  {drug}: {r} R / {len(selected) - r} S")


# ---------------------------------------------------------------- CLI
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_sum = sub.add_parser("summary", help="exploration report for one species")
    p_sum.add_argument("--amr", default="data/metadata/PATRIC_genomes_AMR.txt")
    p_sum.add_argument("--species", default="Klebsiella pneumoniae")
    p_sum.add_argument("--top", type=int, default=25)

    p_lab = sub.add_parser("labels", help="build the label table for the chosen drugs")
    p_lab.add_argument("--amr", default="data/metadata/PATRIC_genomes_AMR.txt")
    p_lab.add_argument("--out", default="data/metadata/kpn_labels.tsv")

    p_qc = sub.add_parser("qc-report", help="QC distributions for labelled genomes")
    p_qc.add_argument("--labels", default="data/metadata/kpn_labels.tsv")
    p_qc.add_argument("--meta", default="data/metadata/genome_metadata")
    p_qc.add_argument("--out", default="data/metadata/kpn_genome_metadata.tsv")

    p_sel = sub.add_parser("select", help="QC, dedup and sample the final genome set")
    p_sel.add_argument("--labels", default="data/metadata/kpn_labels.tsv")
    p_sel.add_argument("--meta", default="data/metadata/kpn_genome_metadata.tsv")
    p_sel.add_argument("--out", default="data/metadata/selected_genomes.tsv")
    p_sel.add_argument("--ids-out", default="data/metadata/genome_ids.txt")
    p_sel.add_argument("--n", type=int, default=SAMPLE_SIZE)
    p_sel.add_argument("--max-per-project", type=int, default=MAX_PER_BIOPROJECT)
    p_sel.add_argument("--seed", type=int, default=SEED)

    args = parser.parse_args()

    if args.command == "select":
        wide = labels_wide(load_labels(args.labels)).dropna()
        meta = read_kpn_metadata(args.meta)
        pool = meta.loc[meta["genome_id"].isin(wide.index)].reset_index(drop=True)
        flags = qc_flags(pool)
        passed = pool.loc[~flags.any(axis=1)]
        deduped = dedup_biosample(passed)
        sampled = sample_genomes(deduped, args.n, args.max_per_project, args.seed)
        keep_cols = ["genome_id", "genome_name", "biosample_accession",
                     "bioproject_accession", "contigs", "genome_length", "gc_content"]
        selected = sampled[keep_cols].merge(
            wide.astype(int), left_on="genome_id", right_index=True, how="left"
        )
        selected.to_csv(args.out, sep="\t", index=False)
        Path(args.ids_out).write_text("\n".join(selected["genome_id"]) + "\n")
        print(f"wrote {len(selected):,} genomes to {args.out} and {args.ids_out}\n")
        print_selection_report(pool, flags, deduped, selected, DRUGS)
        return

    if args.command == "summary":
        print_summary(load_amr_table(args.amr), args.species, args.top)
    elif args.command == "labels":
        df = load_amr_table(args.amr)
        labels = build_labels(filter_taxa(df, TAXON_IDS), DRUGS)
        labels.to_csv(args.out, sep="\t", index=False)
        print(f"wrote {len(labels):,} labels to {args.out}\n")
        print_label_report(df, labels, DRUGS)
    elif args.command == "qc-report":
        labels = load_labels(args.labels)
        meta = load_genome_metadata(args.meta, set(labels["genome_id"]))
        meta.to_csv(args.out, sep="\t", index=False)
        print(f"wrote metadata for {len(meta):,} genomes to {args.out}\n")
        print_qc_report(meta, labels, DRUGS)


if __name__ == "__main__":
    main()
