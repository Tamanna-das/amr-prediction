"""Build a genome x AMR-element presence/absence matrix from AMRFinderPlus output.

One AMRFinderPlus TSV per genome (results/amrfinder/<genome_id>.tsv, run with
--name <genome_id>). A genome with no hits still has a header-only file and
becomes a row of zeros, so it is not silently dropped.

Features are AMRFinderPlus "Element symbol" values with Type == AMR: acquired
genes (e.g. blaKPC-2) and point mutations (e.g. gyrA_S83I). Virulence, stress,
metal and biocide hits are excluded because they are not expected to cause
antibiotic resistance.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

REQUIRED = ["Element symbol", "Type", "Subtype", "Class", "Subclass"]


def read_amrfinder(path: str | Path) -> pd.DataFrame:
    """Read one AMRFinderPlus report; raise with the real header if columns differ."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}; header is {list(df.columns)}")
    return df


def collect_hits(genome_ids: list[str], result_dir: str | Path) -> pd.DataFrame:
    """Long table of AMR hits (genome_id, element, subtype, class, subclass)."""
    parts = []
    for gid in genome_ids:
        df = read_amrfinder(Path(result_dir) / f"{gid}.tsv")
        df = df.loc[df["Type"] == "AMR", ["Element symbol", "Subtype", "Class", "Subclass"]]
        df.insert(0, "genome_id", gid)
        parts.append(df)
    hits = pd.concat(parts, ignore_index=True)
    return hits.rename(columns={"Element symbol": "element", "Subtype": "subtype",
                                "Class": "class", "Subclass": "subclass"})


def presence_matrix(hits: pd.DataFrame, genome_ids: list[str]) -> pd.DataFrame:
    """Genome x element 0/1 matrix; every genome in genome_ids gets a row."""
    pairs = hits[["genome_id", "element"]].drop_duplicates()
    pairs = pairs.assign(present=1)
    mat = pairs.pivot(index="genome_id", columns="element", values="present")
    mat = mat.reindex(index=genome_ids).fillna(0).astype("int8")
    mat.index.name = "genome_id"
    return mat.sort_index(axis=1)


def element_annotation(hits: pd.DataFrame) -> pd.DataFrame:
    """One row per element: AMRFinderPlus subtype/class/subclass and genome count."""
    ann = (hits.groupby("element")
           .agg(subtype=("subtype", "first"), amr_class=("class", "first"),
                amr_subclass=("subclass", "first"), n_genomes=("genome_id", "nunique"))
           .sort_values("n_genomes", ascending=False))
    return ann


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ids", default="data/metadata/genome_ids_500.txt")
    p.add_argument("--results", default="results/amrfinder")
    p.add_argument("--out", default="data/features/amr_matrix.tsv")
    p.add_argument("--annotation", default="data/features/element_annotation.tsv")
    p.add_argument("--allow-missing", action="store_true",
                   help="skip genomes whose AMRFinderPlus output is not there yet")
    args = p.parse_args()

    ids = [l.strip() for l in open(args.ids) if l.strip()]
    present = [i for i in ids if (Path(args.results) / f"{i}.tsv").exists()]
    if len(present) < len(ids) and not args.allow_missing:
        raise SystemExit(f"{len(ids) - len(present)} genomes have no AMRFinderPlus output "
                         f"yet; finish the run or pass --allow-missing")

    hits = collect_hits(present, args.results)
    mat = presence_matrix(hits, present)
    ann = element_annotation(hits)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    mat.to_csv(args.out, sep="\t")
    ann.to_csv(args.annotation, sep="\t")

    print(f"genomes: {mat.shape[0]}   AMR elements: {mat.shape[1]}")
    print(f"genomes with zero AMR hits: {int((mat.sum(axis=1) == 0).sum())}")
    print(f"AMR elements per genome: median {int(mat.sum(axis=1).median())}, "
          f"max {int(mat.sum(axis=1).max())}")
    print("\n# subtypes")
    print(ann["subtype"].value_counts().to_string())
    print("\n# 15 most common elements")
    print(ann.head(15).to_string())


if __name__ == "__main__":
    main()
