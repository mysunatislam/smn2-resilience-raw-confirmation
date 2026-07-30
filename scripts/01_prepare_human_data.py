from __future__ import annotations

import re

import pandas as pd

from common import ROOT, ensure_output_dirs, parse_geo_soft, sha256, write_manifest


def prepare_gse93939() -> dict[str, object]:
    raw_path = ROOT / "data/raw/GSE93939/GSE93939_GEO_counts_human_LCM.xlsx"
    soft_path = ROOT / "data/metadata/GSE93939_family.soft.gz"
    case_path = (
        ROOT / "data/metadata/PMC6565614_supplementary/mmc3.xlsx"
    )

    counts = pd.read_excel(raw_path, sheet_name="counts", index_col=0)
    counts.index = counts.index.astype(str)
    if counts.index.duplicated().any():
        counts = counts.groupby(level=0, sort=False).sum()

    metadata = parse_geo_soft(soft_path)
    metadata["sample_id"] = metadata["title"].str.extract(r"\[([^]]+)\]")
    metadata["group"] = metadata["sample_id"].str.extract(r"^(OMN|SC|Onuf)")
    metadata["donor_id"] = metadata[
        "case number (refer to table in paper)"
    ].astype(str)
    metadata["platform"] = metadata["platform_id"].map(
        {"GPL11154": "HiSeq2000", "GPL16791": "HiSeq2500"}
    )

    case_metadata = pd.read_excel(
        case_path,
        sheet_name="Human samples",
        skiprows=1,
        nrows=22,
    )
    case_metadata.columns = [
        "donor_id",
        "sex",
        "age_at_death",
        "postmortem_delay",
        "tissue_source",
    ]
    case_metadata["donor_id"] = (
        case_metadata["donor_id"].astype(int).astype(str)
    )
    delay_parts = case_metadata["postmortem_delay"].astype(str).str.split(":")
    case_metadata["postmortem_delay_hours"] = (
        delay_parts.str[0].astype(float)
        + delay_parts.str[1].astype(float) / 60.0
    )
    metadata = metadata.merge(case_metadata, on="donor_id", how="left")
    if metadata["sex"].isna().any():
        raise ValueError("Missing GSE93939 donor covariates after case-table merge")
    metadata = metadata.set_index("sample_id", drop=False)

    missing_metadata = counts.columns.difference(metadata.index)
    missing_counts = metadata.index.difference(counts.columns)
    if len(missing_metadata) or len(missing_counts):
        raise ValueError(
            f"GSE93939 sample mismatch: metadata missing {missing_metadata.tolist()}, "
            f"counts missing {missing_counts.tolist()}"
        )

    metadata = metadata.loc[counts.columns]
    counts.to_csv(
        ROOT / "data/processed/GSE93939_counts.csv.gz", compression="gzip"
    )
    metadata.to_csv(ROOT / "data/processed/GSE93939_metadata.tsv", sep="\t")

    return {
        "accession": "GSE93939",
        "raw_file": str(raw_path.relative_to(ROOT)),
        "raw_sha256": sha256(raw_path),
        "genes": int(counts.shape[0]),
        "samples": int(counts.shape[1]),
        "groups": metadata["group"].value_counts().to_dict(),
        "unique_donors": int(metadata["donor_id"].nunique()),
    }


def prepare_hgnc_annotation() -> dict[str, object]:
    raw_path = ROOT / "data/metadata/hgnc_complete_set.txt"
    annotation = pd.read_csv(raw_path, sep="\t", low_memory=False)
    annotation = annotation[
        ["symbol", "name", "locus_group", "locus_type", "location"]
    ].drop_duplicates("symbol")
    annotation["chromosome"] = annotation["location"].str.extract(
        r"^(X|Y|MT|\d+)", expand=False
    )
    annotation = annotation.set_index("symbol", drop=True)
    annotation.to_csv(
        ROOT / "data/processed/HGNC_gene_annotation.tsv", sep="\t"
    )
    return {
        "resource": "HGNC complete set",
        "raw_file": str(raw_path.relative_to(ROOT)),
        "raw_sha256": sha256(raw_path),
        "genes": int(annotation.shape[0]),
    }


def prepare_gse290979() -> dict[str, object]:
    raw_path = ROOT / "data/raw/GSE290979/GSE290979_count_matrix.txt.gz"
    soft_path = ROOT / "data/metadata/GSE290979_family.soft.gz"

    counts = pd.read_csv(raw_path, sep=r"\s+", index_col=0)
    counts.index = counts.index.astype(str)
    if counts.index.duplicated().any():
        counts = counts.groupby(level=0, sort=False).sum()

    metadata = parse_geo_soft(soft_path)
    metadata["sample_id"] = metadata["description"]
    metadata["replicate"] = (
        metadata["title"].str.extract(r"replicate\s+(\d+)", flags=re.IGNORECASE)[0]
    )
    metadata["analysis_group"] = "CTRL_NT"
    metadata.loc[metadata["genotype"].eq("SMA"), "analysis_group"] = (
        "SMA_" + metadata.loc[metadata["genotype"].eq("SMA"), "treatment"]
    )
    metadata = metadata.set_index("sample_id", drop=False)

    missing_metadata = counts.columns.difference(metadata.index)
    missing_counts = metadata.index.difference(counts.columns)
    if len(missing_metadata) or len(missing_counts):
        raise ValueError(
            f"GSE290979 sample mismatch: metadata missing {missing_metadata.tolist()}, "
            f"counts missing {missing_counts.tolist()}"
        )

    metadata = metadata.loc[counts.columns]
    counts.to_csv(
        ROOT / "data/processed/GSE290979_counts.csv.gz", compression="gzip"
    )
    metadata.to_csv(ROOT / "data/processed/GSE290979_metadata.tsv", sep="\t")

    return {
        "accession": "GSE290979",
        "raw_file": str(raw_path.relative_to(ROOT)),
        "raw_sha256": sha256(raw_path),
        "genes": int(counts.shape[0]),
        "samples": int(counts.shape[1]),
        "groups": metadata["analysis_group"].value_counts().to_dict(),
        "donor_lines": int(metadata["cell line"].nunique()),
    }


def main() -> None:
    ensure_output_dirs()
    manifest = [
        prepare_gse93939(),
        prepare_gse290979(),
        prepare_hgnc_annotation(),
    ]
    write_manifest(manifest, ROOT / "data/processed/human_data_manifest.json")
    for entry in manifest:
        if "accession" in entry:
            print(
                f"{entry['accession']}: {entry['genes']} genes, "
                f"{entry['samples']} samples"
            )
        else:
            print(f"{entry['resource']}: {entry['genes']} annotated genes")


if __name__ == "__main__":
    main()
