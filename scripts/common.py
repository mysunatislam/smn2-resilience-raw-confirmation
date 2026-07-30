from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]


def ensure_output_dirs() -> None:
    for relative in [
        "data/processed",
        "results/qc",
        "results/differential_expression",
        "results/pathway_enrichment",
        "results/splicing",
        "results/candidate_ranking",
        "results/figures",
        "docs",
    ]:
        (ROOT / relative).mkdir(parents=True, exist_ok=True)


def parse_geo_soft(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    record: dict[str, str] | None = None

    with gzip.open(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if record:
                    records.append(record)
                record = {"gsm": line.split("=", 1)[1].strip()}
                continue

            if record is None or not line.startswith("!Sample_") or " = " not in line:
                continue

            key, value = line.split(" = ", 1)
            key = key.removeprefix("!Sample_").strip()
            value = value.strip()
            if key == "characteristics_ch1" and ":" in value:
                characteristic, characteristic_value = value.split(":", 1)
                record[characteristic.strip().lower()] = characteristic_value.strip()
            elif key in {
                "title",
                "source_name_ch1",
                "platform_id",
                "instrument_model",
                "description",
            }:
                record[key] = value

    if record:
        records.append(record)

    return pd.DataFrame.from_records(records)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest(entries: Iterable[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(list(entries), handle, indent=2)
        handle.write("\n")


def bh_fdr(p_values: np.ndarray | pd.Series) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return adjusted

    p = values[valid]
    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    restored = np.empty_like(q)
    restored[order] = q
    adjusted[valid] = restored
    return adjusted


def poscounts_size_factors(counts: pd.DataFrame) -> pd.Series:
    values = counts.to_numpy(dtype=float)
    positive_genes = (values > 0).any(axis=1)
    geometric_means = np.full(values.shape[0], np.nan, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        positive_values = values[positive_genes]
        log_values = np.where(positive_values > 0, np.log(positive_values), np.nan)
        geometric_means[positive_genes] = np.exp(np.nanmean(log_values, axis=1))

    valid_genes = np.isfinite(geometric_means) & (geometric_means > 0)
    ratios = values[valid_genes] / geometric_means[valid_genes, None]
    ratios[ratios <= 0] = np.nan
    size_factors = np.nanmedian(ratios, axis=0)

    if not np.all(np.isfinite(size_factors)) or np.any(size_factors <= 0):
        library_sizes = values.sum(axis=0)
        size_factors = library_sizes / np.exp(np.mean(np.log(library_sizes)))

    size_factors = size_factors / np.exp(np.mean(np.log(size_factors)))
    return pd.Series(size_factors, index=counts.columns, name="size_factor")


def cpm(counts: pd.DataFrame) -> pd.DataFrame:
    library_sizes = counts.sum(axis=0).replace(0, np.nan)
    return counts.divide(library_sizes, axis=1) * 1_000_000


def expression_filter(
    counts: pd.DataFrame, min_cpm: float = 1.0, min_samples: int = 7
) -> pd.Series:
    return (cpm(counts) >= min_cpm).sum(axis=1) >= min_samples


def log2_normalized_counts(
    counts: pd.DataFrame, size_factors: pd.Series, pseudocount: float = 0.5
) -> pd.DataFrame:
    normalized = counts.divide(size_factors, axis=1)
    return np.log2(normalized + pseudocount)


def cluster_robust_ols(
    expression: pd.DataFrame,
    design: pd.DataFrame,
    clusters: pd.Series,
    coefficient: str,
) -> pd.DataFrame:
    sample_order = design.index.tolist()
    y = expression.loc[:, sample_order].T.to_numpy(dtype=float)
    x = design.to_numpy(dtype=float)

    if np.linalg.matrix_rank(x) < x.shape[1]:
        raise ValueError("Design matrix is not full rank")

    xtx_inv = np.linalg.inv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residuals = y - x @ beta
    coefficient_index = design.columns.get_loc(coefficient)
    contrast_row = xtx_inv[coefficient_index, :]

    cluster_values = clusters.loc[sample_order].astype(str).to_numpy()
    unique_clusters = np.unique(cluster_values)
    variance = np.zeros(y.shape[1], dtype=float)
    for cluster in unique_clusters:
        mask = cluster_values == cluster
        leverage = x[mask] @ contrast_row
        cluster_score = leverage @ residuals[mask]
        variance += np.square(cluster_score)

    n_samples, n_parameters = x.shape
    n_clusters = len(unique_clusters)
    correction = (n_clusters / (n_clusters - 1)) * (
        (n_samples - 1) / (n_samples - n_parameters)
    )
    standard_error = np.sqrt(np.maximum(variance * correction, 0))
    effect = beta[coefficient_index]

    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = effect / standard_error
    p_value = 2 * stats.t.sf(np.abs(statistic), df=n_clusters - 1)

    result = pd.DataFrame(
        {
            "effect": effect,
            "standard_error": standard_error,
            "statistic": statistic,
            "p_value": p_value,
        },
        index=expression.index,
    )
    result["q_value"] = bh_fdr(result["p_value"].to_numpy())
    result["n_samples"] = n_samples
    result["n_clusters"] = n_clusters
    return result


def ordinary_ols(
    expression: pd.DataFrame,
    design: pd.DataFrame,
    coefficient: str,
) -> pd.DataFrame:
    sample_order = design.index.tolist()
    y = expression.loc[:, sample_order].T.to_numpy(dtype=float)
    x = design.to_numpy(dtype=float)

    if np.linalg.matrix_rank(x) < x.shape[1]:
        raise ValueError("Design matrix is not full rank")

    xtx_inv = np.linalg.inv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    residuals = y - x @ beta
    degrees_freedom = x.shape[0] - x.shape[1]
    residual_variance = np.square(residuals).sum(axis=0) / degrees_freedom
    coefficient_index = design.columns.get_loc(coefficient)
    standard_error = np.sqrt(
        residual_variance * xtx_inv[coefficient_index, coefficient_index]
    )
    effect = beta[coefficient_index]

    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = effect / standard_error
    p_value = 2 * stats.t.sf(np.abs(statistic), df=degrees_freedom)

    result = pd.DataFrame(
        {
            "effect": effect,
            "standard_error": standard_error,
            "statistic": statistic,
            "p_value": p_value,
        },
        index=expression.index,
    )
    result["q_value"] = bh_fdr(result["p_value"].to_numpy())
    result["n_samples"] = x.shape[0]
    result["residual_df"] = degrees_freedom
    return result


def ranked_percentile(values: pd.Series, ascending: bool = True) -> pd.Series:
    return values.rank(method="average", pct=True, ascending=ascending).fillna(0.0)
