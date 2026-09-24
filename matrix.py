from pathlib import Path

import pandas as pd
import numpy as np

from kb_pilot import APP_DATA_DIR

def load_kb_counts(counts_dir: Path | None = None) -> pd.DataFrame:
    if counts_dir is None:
        counts_dir = APP_DATA_DIR / "counts_out_workdir" / "counts_out"

    counts_dir = Path(counts_dir)
    t2g_candidates = list((APP_DATA_DIR / "kb_work").glob("*_t2g.txt"))
    if not t2g_candidates:
        raise FileNotFoundError(f"No t2g file found in {APP_DATA_DIR / 'kb_work'}")

    t2g_file = t2g_candidates[0]

    merged = None

    for sample_dir in sorted(counts_dir.iterdir()):
        if not sample_dir.is_dir():
            continue

        abundance_file = sample_dir / "abundance.tsv"
        if not abundance_file.exists():
            continue

        df = pd.read_csv(abundance_file, sep="\t")
        df = df[["target_id", "est_counts"]].copy()
        df = df.rename(columns={"est_counts": sample_dir.name})

        if merged is None:
            merged = df
        else:
            merged = merged.merge(df, on="target_id", how="outer")

    if merged is None:
        raise FileNotFoundError(f"No abundance.tsv files found in {counts_dir}")

    gene_map = {}
    with open(t2g_file, "r") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                gene_map[parts[0]] = parts[2]

    merged["gene_name"] = merged["target_id"].map(gene_map)
    merged = merged.dropna(subset=["gene_name"])
    merged = merged.drop(columns=["target_id"])

    counts = merged.groupby("gene_name").sum()
    return counts

def normalize_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib_sizes = counts.sum(axis=0)
    return counts.div(lib_sizes, axis=1) * 1_000_000

def log_transform(cpm: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    return np.log2(cpm + pseudocount)

def export_matrix(df: pd.DataFrame, out_path: str | Path) -> None:
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path)