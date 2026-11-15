import pytest

import pandas as pd

from src.select_genomes import (
    AMR_COLUMNS,
    META_COLUMNS,
    biosample_duplicates,
    build_labels,
    conflicting_pairs,
    dedup_biosample,
    filter_taxa,
    load_amr_table,
    load_genome_metadata,
    phenotype_counts,
    qc_flags,
    sample_genomes,
)


def _row(genome_id, drug, phenotype, name='"""Klebsiella pneumoniae """', taxon="573"):
    row = [""] * len(AMR_COLUMNS)
    row[0], row[1], row[2], row[3], row[4] = genome_id, name, taxon, drug, phenotype
    return row


def _write(tmp_path, rows, columns=AMR_COLUMNS, name="amr.tsv"):
    path = tmp_path / name
    lines = ["\t".join(columns)] + ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------- loading
def test_genome_id_stays_text(tmp_path):
    df = load_amr_table(_write(tmp_path, [_row("573.12340", "meropenem", "Resistant")]))
    assert df.loc[0, "genome_id"] == "573.12340"


def test_genome_name_quotes_and_spaces_removed(tmp_path):
    df = load_amr_table(_write(tmp_path, [_row("573.1", "meropenem", "Resistant")]))
    assert df.loc[0, "genome_name"] == "Klebsiella pneumoniae"


def test_empty_phenotype_is_missing_and_values_lowercased(tmp_path):
    df = load_amr_table(_write(tmp_path, [
        _row("573.1", "Meropenem ", ""),
        _row("573.2", "meropenem", "Resistant"),
    ]))
    assert df["resistant_phenotype"].isna().sum() == 1
    assert set(df["antibiotic"]) == {"meropenem"}
    assert df.loc[1, "resistant_phenotype"] == "resistant"


def test_missing_column_raises(tmp_path):
    path = _write(tmp_path, [["573.1"]], columns=["genome_id"])
    with pytest.raises(ValueError, match="missing columns"):
        load_amr_table(path)


# ---------------------------------------------------------------- summaries
def test_phenotype_counts_counts_each_genome_once(tmp_path):
    df = load_amr_table(_write(tmp_path, [
        _row("573.1", "meropenem", "Resistant"),
        _row("573.1", "meropenem", "Resistant"),
        _row("573.2", "meropenem", "Susceptible"),
    ]))
    counts = phenotype_counts(df)
    assert counts.loc["meropenem", "resistant"] == 1
    assert counts.loc["meropenem", "susceptible"] == 1


def test_conflicting_pairs_found(tmp_path):
    df = load_amr_table(_write(tmp_path, [
        _row("573.1", "meropenem", "Resistant"),
        _row("573.1", "meropenem", "Susceptible"),
        _row("573.2", "meropenem", "Resistant"),
    ]))
    conflicts = conflicting_pairs(df)
    assert list(conflicts["genome_id"]) == ["573.1"]


# ---------------------------------------------------------------- labelling
def test_filter_taxa_keeps_only_allowed_ids(tmp_path):
    df = load_amr_table(_write(tmp_path, [
        _row("573.1", "meropenem", "Resistant", taxon="573"),
        _row("573.2", "meropenem", "Resistant", taxon="72407"),
        _row("573.3", "meropenem", "Resistant", taxon="1128952"),
    ]))
    kept = filter_taxa(df, ["573", "72407"])
    assert sorted(kept["genome_id"]) == ["573.1", "573.2"]


def test_build_labels_maps_r_s_and_drops_intermediate(tmp_path):
    df = load_amr_table(_write(tmp_path, [
        _row("573.1", "meropenem", "Resistant"),
        _row("573.2", "meropenem", "Susceptible"),
        _row("573.3", "meropenem", "Intermediate"),
        _row("573.4", "meropenem", ""),
    ]))
    labels = build_labels(df, ["meropenem"])
    assert dict(zip(labels["genome_id"], labels["label"])) == {"573.1": 1, "573.2": 0}


def test_build_labels_drops_only_the_conflicting_pair(tmp_path):
    df = load_amr_table(_write(tmp_path, [
        _row("573.1", "meropenem", "Resistant"),
        _row("573.1", "meropenem", "Intermediate"),   # conflict -> drop pair
        _row("573.1", "ciprofloxacin", "Resistant"),  # same genome, other drug: keep
    ]))
    labels = build_labels(df, ["meropenem", "ciprofloxacin"])
    assert list(zip(labels["genome_id"], labels["antibiotic"])) == [
        ("573.1", "ciprofloxacin")
    ]


def test_build_labels_ignores_other_drugs(tmp_path):
    df = load_amr_table(_write(tmp_path, [
        _row("573.1", "meropenem", "Resistant"),
        _row("573.1", "colistin", "Resistant"),
    ]))
    labels = build_labels(df, ["meropenem"])
    assert set(labels["antibiotic"]) == {"meropenem"}


# ---------------------------------------------------------------- QC
def _meta_row(genome_id, biosample="", contigs="50", length="5500000"):
    values = {c: "" for c in META_COLUMNS}
    values.update(genome_id=genome_id, genome_name="Klebsiella pneumoniae x",
                  taxon_id="573", biosample_accession=biosample, contigs=contigs,
                  genome_length=length, gc_content="57.2", patric_cds="5300")
    return [values[c] for c in META_COLUMNS] + ["ignored"]


def _write_meta(tmp_path, rows):
    return _write(tmp_path, rows, columns=META_COLUMNS + ["extra_column"],
                  name="genome_metadata")


def test_load_genome_metadata_keeps_requested_ids_and_casts_numbers(tmp_path):
    path = _write_meta(tmp_path, [
        _meta_row("573.12340", contigs="120"),
        _meta_row("573.999"),
    ])
    meta = load_genome_metadata(path, {"573.12340"}, chunksize=1)
    assert list(meta["genome_id"]) == ["573.12340"]
    assert meta.loc[0, "contigs"] == 120
    assert "extra_column" not in meta.columns


def test_biosample_duplicates_finds_shared_biosamples(tmp_path):
    path = _write_meta(tmp_path, [
        _meta_row("573.1", biosample="SAMN1"),
        _meta_row("573.2", biosample="SAMN1"),
        _meta_row("573.3", biosample="SAMN2"),
        _meta_row("573.4", biosample=""),
    ])
    meta = load_genome_metadata(path, {"573.1", "573.2", "573.3", "573.4"})
    dups = biosample_duplicates(meta)
    assert sorted(dups["genome_id"]) == ["573.1", "573.2"]


# ---------------------------------------------------------------- selection
def _meta(rows):
    """rows: dicts with any of genome_id, status, length, gc, contigs, biosample, project."""
    return pd.DataFrame({
        "genome_id": [r["genome_id"] for r in rows],
        "genome_status": [r.get("status", "WGS") for r in rows],
        "genome_length": [r.get("length", 5_500_000) for r in rows],
        "gc_content": [r.get("gc", 57.2) for r in rows],
        "contigs": [r.get("contigs", 100) for r in rows],
        "biosample_accession": [r.get("biosample") for r in rows],
        "bioproject_accession": [r.get("project") for r in rows],
    })


def test_qc_flags_each_rule():
    meta = _meta([
        {"genome_id": "ok"},
        {"genome_id": "short", "length": 17_021},
        {"genome_id": "gc", "gc": 46.35},
        {"genome_id": "frag", "contigs": 5155},
        {"genome_id": "nocontigs", "contigs": float("nan")},
        {"genome_id": "dep", "status": "Deprecated"},
    ])
    failed = meta.loc[qc_flags(meta).any(axis=1), "genome_id"]
    assert sorted(failed) == ["dep", "frag", "gc", "nocontigs", "short"]


def test_dedup_biosample_keeps_fewest_contigs_and_keeps_missing():
    meta = _meta([
        {"genome_id": "573.1", "biosample": "SAMN1", "contigs": 300},
        {"genome_id": "573.2", "biosample": "SAMN1", "contigs": 80},
        {"genome_id": "573.3", "biosample": None},
        {"genome_id": "573.4", "biosample": None},
    ])
    assert list(dedup_biosample(meta)["genome_id"]) == ["573.2", "573.3", "573.4"]


def test_sample_genomes_respects_cap_and_is_reproducible():
    rows = [{"genome_id": f"573.{i}", "project": "BIG"} for i in range(50)]
    rows += [{"genome_id": f"573.{100 + i}", "project": f"P{i}"} for i in range(20)]
    meta = _meta(rows)
    a = sample_genomes(meta, n=30, max_per_project=10, seed=42)
    b = sample_genomes(meta, n=30, max_per_project=10, seed=42)
    assert len(a) == 30
    assert (a["bioproject_accession"] == "BIG").sum() <= 10
    assert list(a["genome_id"]) == list(b["genome_id"])


def test_sample_genomes_raises_when_cap_too_tight():
    meta = _meta([{"genome_id": f"573.{i}", "project": "BIG"} for i in range(50)])
    with pytest.raises(ValueError, match="after capping"):
        sample_genomes(meta, n=30, max_per_project=10, seed=42)
