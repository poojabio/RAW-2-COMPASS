from ast import arg
import scipy.io
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
#from inmoose.pycombat import pycombat_norm
from kb_script import APP_DATA_DIR



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


def export_matrix(df: pd.DataFrame, out_path: Path) -> None:
    df.to_csv(out_path)
    print(f"Exported to {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts-dir", default= APP_DATA_DIR/"counts_out_workdir/counts_out/counts_unfiltered")
    args = parser.parse_args()
    counts = load_kb_counts(Path(args.counts_dir))

    print(f"Shape: {counts.shape[0]} genes x {counts.shape[1]} samples")
    print(f"\nTotal counts per sample:\n{counts.sum(axis=0)}")
    print(f"\nGenes with nonzero counts: {(counts.sum(axis=1) > 0).sum()} / {len(counts)}")

    cpm = normalize_cpm(counts)
    log_cpm = log_transform(cpm)

    print(f"\nCPM (nonzero genes):\n{cpm[cpm.sum(axis=1) > 0]}")
    print(f"\nLog2(CPM+1) (nonzero genes):\n{log_cpm[cpm.sum(axis=1) > 0]}")

    export_matrix(cpm, Path("cpm_matrix.csv"))
    export_matrix(log_cpm, Path("log_cpm_matrix.csv"))

