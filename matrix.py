import scipy.io
import numpy as np
import pandas as pd
from pathlib import Path
from inmoose.pycombat import pycombat_norm


def load_kb_counts(counts_dir: Path) -> pd.DataFrame:
    mtx_path = counts_dir / "cells_x_genes.mtx"
    barcodes_path = counts_dir / "cells_x_genes.barcodes.txt"
    genes_path = counts_dir / "cells_x_genes.genes.names.txt"

    matrix = scipy.io.mmread(mtx_path).tocsr()
    barcodes = pd.read_csv(barcodes_path, header=None)[0].tolist()
    genes = pd.read_csv(genes_path, header=None)[0].tolist()

    counts = pd.DataFrame(matrix.T.toarray(), index=genes, columns=barcodes)
    return counts


def normalize_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """Counts per million, using library size from the raw counts matrix."""
    lib_sizes = counts.sum(axis=0)
    return counts.div(lib_sizes, axis=1) * 1e6


def log_transform(cpm: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    return np.log2(cpm + pseudocount)


def batch_correct(log_cpm: pd.DataFrame, batch_labels: dict) -> pd.DataFrame:
    """
    Optional step — only called if the user opts in.
    batch_labels: {sample_id: batch_name}, must cover every column in log_cpm.
    Requires at least 2 samples per batch to estimate anything meaningful.
    """
    missing = set(log_cpm.columns) - set(batch_labels.keys())
    if missing:
        raise ValueError(f"Missing batch label for samples: {missing}")

    batch = pd.Series(batch_labels)[log_cpm.columns]

    batch_counts = batch.value_counts()
    if (batch_counts < 2).any():
        raise ValueError(
            f"Each batch needs >=2 samples to correct. Got: {batch_counts.to_dict()}"
        )

    corrected = pycombat_norm(log_cpm, batch)
    return corrected


def export_matrix(df: pd.DataFrame, out_path: Path) -> None:
    df.to_csv(out_path)
    print(f"Exported to {out_path}")


if __name__ == "__main__":
    counts = load_kb_counts(Path("counts_out/counts_unfiltered"))

    print(f"Shape: {counts.shape[0]} genes x {counts.shape[1]} samples")
    print(f"\nTotal counts per sample:\n{counts.sum(axis=0)}")
    print(f"\nGenes with nonzero counts: {(counts.sum(axis=1) > 0).sum()} / {len(counts)}")

    cpm = normalize_cpm(counts)
    log_cpm = log_transform(cpm)

    print(f"\nCPM (nonzero genes):\n{cpm[cpm.sum(axis=1) > 0]}")
    print(f"\nLog2(CPM+1) (nonzero genes):\n{log_cpm[cpm.sum(axis=1) > 0]}")

    export_matrix(cpm, Path("cpm_matrix.csv"))
    export_matrix(log_cpm, Path("log_cpm_matrix.csv"))

