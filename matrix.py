import numpy as np
import pandas as pd
from pathlib import Path
import requests
import time
import threading
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
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


def batch_lookup(ids: list[str]) -> dict:
    """POST up to ~1000 Ensembl IDs at once, return {id: lookup_data_or_None}."""
    url = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    response = requests.post(url, headers=headers, json={"ids": ids})
    response.raise_for_status()
    return response.json()


def load_kb_counts(gene_map_cache: Path = Path("gene_map_cache.json")) -> pd.DataFrame:
    kallisto_out_dir = APP_DATA_DIR / "counts_out_workdir" / "counts_out"

    # ========================================================
    # DEBUG SECTIONAL 1: read + accumulate each sample's counts
    # ========================================================
    print("[1/4] Reading sample abundance files...")
    prev = pd.DataFrame(columns=["target_id"])
    n_samples = 0

    for obj in kallisto_out_dir.iterdir():
        if obj.is_dir():
            file = obj / "abundance.tsv"
            df = pd.read_csv(file, sep="\t")
            df = df[["target_id", "est_counts"]]
            df = df.rename(columns={"est_counts": f"est_counts_{obj.name}"})
            prev = pd.merge(prev, df, on="target_id", how="outer")
            n_samples += 1
            print(f"    merged sample: {obj.name}  (running shape: {prev.shape})")

    df = prev
    print(f"[1/4] Done. {n_samples} samples merged. Combined shape: {df.shape}")

    # ========================================================
    # DEBUG SECTIONAL 2: clean up target_id
    # ========================================================
    print("[2/4] Cleaning target_id (stripping version suffixes)...")
    df["target_id"] = df["target_id"].str.replace(r'\.\d+$', '', regex=True)
    df["gene_name"] = "N/A"

    all_ids = df["target_id"].unique().tolist()
    print(f"[2/4] Done. {len(all_ids)} unique transcript IDs to resolve.")

    # ========================================================
    # DEBUG SECTIONAL 3: batched + threaded Ensembl lookups
    #   - sized for large transcript counts (~460k+)
    #   - caches gene_map to disk so a failed/interrupted run
    #     doesn't force redoing all network calls
    # ========================================================
    if gene_map_cache.exists():
        print(f"[3/4] Found cached gene map at {gene_map_cache}, loading instead of querying API.")
        with open(gene_map_cache) as f:
            gene_map = json.load(f)
        print(f"[3/4] Loaded {len(gene_map)} cached gene mappings.")
    else:
        chunk_size = 1000
        gene_map = {}
        lock = threading.Lock()  # only one thread writes to gene_map at a time

        chunks = [all_ids[i:i + chunk_size] for i in range(0, len(all_ids), chunk_size)]
        total_chunks = len(chunks)
        print(f"[3/4] Querying Ensembl: {total_chunks} chunks of up to {chunk_size} IDs each.")

        def worker(chunk, chunk_num):
            try:
                result = batch_lookup(chunk)
                with lock:
                    gene_map.update(result)
                    if chunk_num % 10 == 0 or chunk_num == total_chunks:
                        print(f"    [{chunk_num}/{total_chunks} chunks done] "
                              f"{len(gene_map)} genes mapped so far")
            except requests.exceptions.RequestException as e:
                print(f"    chunk {chunk_num}/{total_chunks} (size {len(chunk)}) FAILED: {e}")
            time.sleep(1)  # throttle, per-worker rather than per-sequential-step

        max_workers = min(5, total_chunks) or 1  # stay under Ensembl's rate limit; guard empty case

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(worker, chunk, i + 1) for i, chunk in enumerate(chunks)]
            for future in as_completed(futures):
                future.result()  # re-raises anything that slipped past the try/except

        print(f"[3/4] Done. {len(gene_map)}/{len(all_ids)} IDs resolved "
              f"({len(all_ids) - len(gene_map)} unresolved, will fall back to raw ID).")

        # cache immediately, before any further processing that could fail
        with open(gene_map_cache, "w") as f:
            json.dump(gene_map, f)
        print(f"[3/4] Cached gene map to {gene_map_cache}")

    # ========================================================
    # DEBUG SECTIONAL 4: map gene names, collapse duplicates
    # ========================================================
    print("[4/4] Mapping gene names and collapsing duplicates...")
    df["gene_name"] = df["target_id"].map(
        lambda tid: (gene_map.get(tid) or {}).get("display_name", tid)
        if gene_map.get(tid) is not None else tid
    )

    df = df.drop(columns=["target_id"])
    n_before = len(df)
    df = df.groupby("gene_name").sum()
    print(f"[4/4] Done. Collapsed {n_before} rows -> {len(df)} unique genes.")

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

    export_matrix(cpm, Path("cpm_matrix.csv"))
    export_matrix(log_cpm, Path("log_cpm_matrix.csv"))