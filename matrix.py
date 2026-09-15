import numpy as np
import pandas as pd
from pathlib import Path
import os
import glob
#from inmoose.pycombat import pycombat_norm
from kb_script import APP_DATA_DIR

'''
1. go into the workdir count output folder
2. for each sample abundance.tsv, rename its est_counts column to the sample name (from directory)
3. outer-merge each sample's [target_id, est_counts_<sample>] onto the running accumulator, on target_id
4. repeat via a loop across all sample directories
5. strip Ensembl version suffixes, then convert ENSG ids to gene names via batched Ensembl API lookups
6. for all duplicate gene names, sum via .groupby("gene_name").sum()
7. return the final combined counts table: genes x samples, raw counts
'''


def load_kb_counts() -> pd.DataFrame:
    kallisto_out_dir = APP_DATA_DIR / "counts_out_workdir" / "counts_out"
    pattern = os.path.join(APP_DATA_DIR,"kb_work", "*.txt")

    matches = glob.glob(pattern)
    t2g_file = matches[0] 

    prev = pd.DataFrame(columns=["target_id"])

    for obj in kallisto_out_dir.iterdir():
        if obj.is_dir():
            file = obj / "abundance.tsv"
            df = pd.read_csv(file, sep="\t")
            df = df[["target_id", "est_counts"]]
            df = df.rename(columns={"est_counts": f"est_counts_{obj.name}"})
            prev = pd.merge(prev, df, on="target_id", how="outer")

    df = prev 

    gene_maps = {}

    with open(t2g_file, "r") as ref:
        file = ref.readlines()
        for item in file:
            splits = item.strip().split("\t")
            print (splits)
            if len(splits) > 2:
                key = splits[0] #ENST
                gene = splits[2] #Gene Name
                gene_maps[key] = gene

    df['gene_id'] = df["target_id"].map(gene_maps)
    df.drop("target_id", axis=1,inplace=True)
    df = df.groupby("gene_id").sum()

    #print(df.head())

    return df


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
    counts = load_kb_counts()

    # sanity checks before trusting the table downstream
    print(f"dtypes:\n{counts.dtypes}")
    print(f"index sample: {counts.index[:10].tolist()}")
    print(f"NaNs per column:\n{counts.isna().sum()}")

    print(f"Shape: {counts.shape[0]} genes x {counts.shape[1]} samples")
    print(f"\nTotal counts per sample:\n{counts.sum(axis=0)}")
    print(f"\nGenes with nonzero counts: {(counts.sum(axis=1) > 0).sum()} / {len(counts)}")

    cpm = normalize_cpm(counts)
    log_cpm = log_transform(cpm)

    print(f"\nCPM (nonzero genes):\n{cpm[cpm.sum(axis=1) > 0]}")
    print(f"\nLog2(CPM+1) (nonzero genes):\n{log_cpm[cpm.sum(axis=1) > 0]}")

    export_matrix(cpm, APP_DATA_DIR/"kb_work"/ ("cpm_matrix.csv"))
    export_matrix(log_cpm, APP_DATA_DIR/"kb_work"/("log_cpm_matrix.csv"))