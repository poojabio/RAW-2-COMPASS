import numpy as np
import pandas as pd
from pathlib import Path
import requests
import time
#from inmoose.pycombat import pycombat_norm
from kb_script import APP_DATA_DIR

'''
1. go into the wordir count output folder 
2. for each sample abundance file.tsv, set the first sample as the constructor foundation for your replace counts column name with sample name (from directory)
3. for each following column save the sample name as the count column, join that revised dataframe on the gene column name or join on that column drop all others 
4. repeat via a loop
5. convert the ENSG ensembl gene ids to their gene name mappings
6. for all duplicate names sum the .groupby(gene_name).sum()
7. return the final overall table

'''

def load_kb_counts() -> pd.DataFrame:
    kallisto_our_dir = APP_DATA_DIR/"counts_out_workdir"/"counts_out"

# step 1 to 4 
    prev = pd.DataFrame(columns = ["target_id"]) ##using the target id column to do join on 
    
    #print("prev's head base df blank")
    ##prev.head()

    for obj in kallisto_our_dir.iterdir():

        if obj.is_dir():
            file = obj/ "abundance.tsv"
            df = pd.read_csv(file, sep="\t")
            df = df[["target_id","est_counts"]]
            df = df.rename(columns={"est_counts": f"est_counts_{obj.name}"}) ## to kep structure w obj name
            ##print("column grouping")
            ##print(df.head())

            ##df = df.join(prev,on="target_id") ##ABORTED this version to stop index based merging, join will allow for column-based
            prev = pd.merge(prev, df, on="target_id", how="outer")

            #print("JOINING with prev")
            ##df.head()

    df = prev

# step 2 mapping gene names
    ## cleaning up the name
    df["target_id"] = df["target_id"].str.replace(r'\.\d+$', '', regex=True)
    df["gene_name"] = "N/A"

    # REST API Calls from ENSEMBL 

    for i, targ_id in df["target_id"].items():
        url = f"https://rest.ensembl.org/lookup/id/{targ_id}"
        headers = {"Content-Type" : "application/json"}

        try:
            response = requests.get(url, headers=headers)
        except requests.exceptions.RequestException as e:
            print(f"Request failed for {targ_id}: {e}")
            df.at[i, "gene_name"] = "N/A"
            time.sleep(0.1)
            continue

        # requests = get, put, patch, post or delete we use get since its an id lookup to fetch
        if response.status_code == 200: ## a successful get 
            data = response.json() ##taking response content into a python dict format json -> dict -> string
            #print(data)
            gene_id = data.get('display_name', "N/A")
            df.at[i, "gene_name"] = gene_id
        else:
            df.at[i, "gene_name"] = targ_id

        time.sleep(0.1)

    df = df.drop(columns=["target_id"])          # drop the versioned/raw ID, keep gene_name
    df = df.groupby("gene_name").sum()            # collapse duplicate gene names, sums their counts

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
    #parser = argparse.ArgumentParser()
    #parser.add_argument("--counts-dir", default= APP_DATA_DIR/"counts_out_workdir/counts_out/counts_unfiltered")
    #args = parser.parse_args()
    counts = load_kb_counts()

    print(f"Shape: {counts.shape[0]} genes x {counts.shape[1]} samples")
    print(f"\nTotal counts per sample:\n{counts.sum(axis=0)}")
    print(f"\nGenes with nonzero counts: {(counts.sum(axis=1) > 0).sum()} / {len(counts)}")

    cpm = normalize_cpm(counts)
    log_cpm = log_transform(cpm)

    print(f"\nCPM (nonzero genes):\n{cpm[cpm.sum(axis=1) > 0]}")
    print(f"\nLog2(CPM+1) (nonzero genes):\n{log_cpm[cpm.sum(axis=1) > 0]}")

    export_matrix(cpm, Path("cpm_matrix.csv"))
    export_matrix(log_cpm, Path("log_cpm_matrix.csv"))

