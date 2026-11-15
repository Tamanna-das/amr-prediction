"""Group genomes into lineages from Mash distances (single-linkage clustering).

Input: the output of `mash dist sketch.msh sketch.msh` (query, reference,
distance, p-value, shared-hashes), where sequence names are genome FASTA paths.
Two genomes end up in the same group if a chain of genomes links them, each
step closer than the threshold. Single linkage is deliberately conservative for
leakage: anything close to anything in a group joins it.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

THRESHOLDS = [0.0005, 0.001, 0.0025, 0.005, 0.01, 0.02]


def genome_id_from_path(path: str) -> str:
    """data/genomes/573.12771.fna -> 573.12771"""
    name = Path(path).name
    return name[:-4] if name.endswith(".fna") else name


def read_mash_dist(path: str | Path) -> pd.DataFrame:
    """Square genome x genome distance matrix from `mash dist` output."""
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["a", "b", "dist", "p", "shared"], dtype={"dist": float})
    df["a"] = df["a"].map(genome_id_from_path)
    df["b"] = df["b"].map(genome_id_from_path)
    mat = df.pivot(index="a", columns="b", values="dist")
    ids = sorted(set(mat.index) | set(mat.columns))
    mat = mat.reindex(index=ids, columns=ids)
    if mat.isna().any().any():
        raise ValueError("distance matrix has gaps: run mash dist on all-vs-all")
    arr = ((mat + mat.T) / 2).to_numpy(copy=True)   # Mash is symmetric up to rounding
    np.fill_diagonal(arr, 0.0)
    return pd.DataFrame(arr, index=ids, columns=ids)


def cluster(dist: pd.DataFrame, threshold: float) -> pd.Series:
    """Single-linkage groups at a distance threshold; labels are 1..k."""
    if len(dist) == 1:
        return pd.Series([1], index=dist.index, name="lineage")
    z = linkage(squareform(dist.values, checks=False), method="single")
    labels = fcluster(z, t=threshold, criterion="distance")
    return pd.Series(labels, index=dist.index, name="lineage")


def threshold_table(dist: pd.DataFrame, thresholds=THRESHOLDS) -> pd.DataFrame:
    rows = []
    for t in thresholds:
        sizes = cluster(dist, t).value_counts()
        rows.append({"threshold": t, "n_groups": len(sizes),
                     "largest_group": int(sizes.max()),
                     "largest_share": round(sizes.max() / sizes.sum(), 3),
                     "singletons": int((sizes == 1).sum())})
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dist", default="results/mash/dist.tsv")
    p.add_argument("--threshold", type=float, default=0.005)
    p.add_argument("--out", default="data/features/lineage_groups.tsv")
    args = p.parse_args()

    dist = read_mash_dist(args.dist)
    print(f"genomes: {len(dist)}")
    off = dist.values[np.triu_indices(len(dist), k=1)]
    print(f"pairwise distance: median {np.median(off):.4f}, "
          f"pairs < 0.0001 (near-identical): {int((off < 0.0001).sum())}")
    print("\n# groups at each threshold")
    print(threshold_table(dist).to_string(index=False))

    groups = cluster(dist, args.threshold)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    groups.rename_axis("genome_id").to_csv(args.out, sep="\t")
    print(f"\nwrote {groups.nunique()} lineage groups at threshold {args.threshold} "
          f"to {args.out}")


if __name__ == "__main__":
    main()
