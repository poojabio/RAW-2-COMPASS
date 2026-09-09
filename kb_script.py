import sys
import subprocess
from pathlib import Path
import time
#import gget
import argparse
import shutil
import gzip
import os

def find_binary(name, base_dir):
    local_path = base_dir / "bin" / name

    # 1. Check local bin/ first
    if local_path.exists() and local_path.is_file(): ##do the scripts and binaries exist
        if not os.access(local_path, os.X_OK): ##is it an executable file or not
            raise RuntimeError(f"{local_path} exists but is not executable")
        return str(local_path)

    # 2. Fallback to system PATH
    path = shutil.which(name)
    if path:
        return path
    raise RuntimeError(f"{name} not found in bin/ or PATH")

BASE_DIR = Path(__file__).resolve().parent.parent ##absolute path 2 levels above the scripts (src > raw2compass)

KALLISTO = find_binary("kallisto", BASE_DIR)
BUSTOOLS = find_binary("bustools", BASE_DIR)

APP_DATA_DIR = BASE_DIR / "workdir"
KB_CMD = [sys.executable, "-m", "kb_python.main"]
APP_DATA_DIR.mkdir(parents=True, exist_ok=True)

print("Using:")
print("  KALLISTO:", KALLISTO)
print("  BUSTOOLS:", BUSTOOLS)
print("  WORKDIR :", APP_DATA_DIR)

def _build_t2g_from_cdna(cdna_fasta: Path, t2g_out: Path) -> None:
    """
    Extracts transcript_id / gene_id / gene_name from Ensembl cDNA FASTA
    headers directly, since --workflow=custom doesn't generate t2g.txt.
    """
    with gzip.open(cdna_fasta, "rt") as f, open(t2g_out, "w") as out:
        for line in f:
            if not line.startswith(">"):
                continue
            fields = line[1:].split()
            tx = fields[0]
            gene = ""
            symbol = ""
            for field in fields:
                if field.startswith("gene:"):
                    gene = field.split(":", 1)[1]
                elif field.startswith("gene_symbol:"):
                    symbol = field.split(":", 1)[1]
            out.write(f"{tx}\t{gene}\t{symbol}\n")

def ref_builder_cdna(species: str) -> tuple[str, str]: ##homo_sapiens, mus_musculus, perhaps integrate zebrafish later too
    """
    Fast path: builds index directly from a pre-made cDNA transcript FASTA
    (e.g. Ensembl's cdna.all.fa.gz), skipping genome+GTF splitting entirely.
    Uses --workflow=custom under kb ref. Currently human/mouse only.
    """
    print("Hello! Starting cDNA Reference Build")
    start = time.time()

    work_dir = APP_DATA_DIR / "kb_work"
    work_dir.mkdir(exist_ok=True) ## if it exists it is idempotent
    index_file = f"{species}_index.idx"
    t2g_file = f"{species}_t2g.txt"

    index_path = work_dir / index_file
    t2g_path = work_dir / t2g_file

    if index_path.exists() and t2g_path.exists() and t2g_path.stat().st_size > 0: ##existing files and checking size for nonemptiness
        print(f"Reference for {species} already built — skipping.")
        return str(index_path), str(t2g_path)

    result = subprocess.run(
                ["gget",
                 "ref",
                 "--ftp",
                 "-w",
                 "cdna",
                 species
                ],check=True,capture_output=True, text=True)

    cdna_url = result.stdout.strip()
    cdna_fasta = work_dir / Path(cdna_url).name
    if not cdna_fasta.exists():
        print(f"Downloading cDNA FASTA for {species}...")
        subprocess.run(["curl", 
                        "-L", 
                        "-o", 
                        str(cdna_fasta),
                        cdna_url],
                        check=True)
        

    subprocess.run([
        *KB_CMD, "ref",
        "--workflow=custom",
        "-i", index_file,
        "-g", t2g_file,
        "--kallisto", KALLISTO,
        "-t" , "8",
        str(cdna_fasta.resolve()),
    ], cwd=work_dir, check=True)

    if not t2g_path.exists() or t2g_path.stat().st_size == 0:
        print("Generating transcript-to-gene mapping from FASTA headers...")
        _build_t2g_from_cdna(cdna_fasta, t2g_path)

    elapsed = time.time() - start
    print("Done! cDNA Reference Build")
    print(f"Took {elapsed:.1f} seconds")

    return str(work_dir / index_file), str(work_dir / t2g_file)



### FILE PAIR IDENTIFICATION AND DOWNSTREAM COUNTS FORMATION
def find_file_pairs(folder: Path) -> tuple[dict, list]:
    """
    Scans `folder` for paired fastqs. Returns (paired, unpaired).
    paired: {sample_id: (fwd_path, rev_path)}
    """
    print("Parsing through for file pairs")
    start = time.time()

    folder = Path(folder).expanduser()
    if not folder.is_dir():
        raise ValueError(f"Not a folder: {folder}")

    paired = {}
    unpaired = []

    PAIR_PATTERNS = [
        ("_R1", "_R2"),
        ("read1", "read2"),
        ("forward", "reverse"),
        ("_1", "_2"),
    ]

    for file in folder.iterdir():
        if not file.is_file():
            continue

        matched = False
        for fwd_token, rev_token in PAIR_PATTERNS:
            if fwd_token in file.name:
                rev_name = file.name.replace(fwd_token, rev_token) ##CREATING REV TOKEN TO SEARCH WITH
                rev = file.parent / rev_name
                if rev.exists():
                    sample = file.name[:file.name.find(fwd_token)] ##Splicing filename 
                    paired[sample] = (file, rev)
                    matched = True
                break  # stop checking other patterns once one token matches

        if not matched:
            unpaired.append(file)

    for sample, (fwd, rev) in paired.items():
        print(f"{sample}: {fwd.name} + {rev.name}")

    print("Done! Pair finding")
    elapsed = time.time() - start
    print(f"Took {elapsed:.1f} seconds")
    return paired, unpaired

def pair_selected_files(file_paths: list[Path]) -> tuple[dict, list]:
    """
    Pairs files from an explicit selection (e.g. multi-select file dialog).
    Two-pass: find all valid pairs first, then anything not claimed is unpaired.
    """
    PAIR_PATTERNS = [
        ("_R1", "_R2"),
        ("read1", "read2"),
        ("forward", "reverse"),
        ("_1", "_2"),
    ]

    file_set = {Path(f) for f in file_paths}
    paired = {}
    claimed = set()

    for file in file_set:
        if file in claimed:
            continue
        for fwd_token, rev_token in PAIR_PATTERNS:
            if fwd_token in file.name:
                rev_name = file.name.replace(fwd_token, rev_token)
                rev = file.parent / rev_name
                if rev in file_set:
                    sample = file.name[:file.name.find(fwd_token)]
                    paired[sample] = (file, rev)
                    claimed.add(file)
                    claimed.add(rev)
                break

    unpaired = [f for f in file_set if f not in claimed]
    return paired, unpaired


def selected_files_as_singles(file_paths: list[Path]) -> dict:
    """Single-end mode: each selected file is its own sample."""
    singles = {}
    for file in file_paths:
        file = Path(file)
        sample = file.stem.replace(".fastq", "").replace(".fq", "")
        singles[sample] = file
    return singles

def find_single_files(folder: Path) -> dict:
    """
    Single-end mode: every fastq in the folder is its own sample.
    {sample_id: file_path}
    """
    folder = Path(folder).expanduser()
    if not folder.is_dir():
        raise ValueError(f"Not a folder: {folder}")

    singles = {}
    for file in folder.iterdir():
        if file.is_file() and file.suffix in (".gz", ".fastq", ".fq"):
            sample = file.stem.replace(".fastq", "").replace(".fq", "")
            singles[sample] = file
    return singles


def write_batch_file(samples: dict, out_path: Path, parity: str) -> Path:
    """
    samples: {sample_id: (fwd, rev)} for paired, or {sample_id: file} for single.
    """
    with open(out_path, "w") as f:
        if parity == "paired":
            for sample_id, (fwd, rev) in samples.items():
                f.write(f"{sample_id}\t{fwd}\t{rev}\n")
        else:
            for sample_id, fwd in samples.items():
                f.write(f"{sample_id}\t{fwd}\n")
    return out_path


def run_count(samples: dict, index: str, t2g: str, parity: str = "paired",
              work_dir: Path=None) -> Path:
    """
    parity: "paired" or "single" — samples dict shape must match (see write_batch_file).
    Returns the output directory path.
    """
    print("Hello! Beginning Run Count function")

    if work_dir is None:
        work_dir = APP_DATA_DIR / "counts_out_workdir"
        work_dir.mkdir(exist_ok=True)

    start = time.time()

    
    batch_file = write_batch_file(samples, work_dir / "batch.txt", parity)
    out_dir = work_dir / "counts_out"

    subprocess.run([
        *KB_CMD, "count",
        "-i", index,
        "-g", t2g,
        "-x", "BULK",
        "--parity", parity,
        "-o", str(out_dir),
        "--kallisto", KALLISTO,
        "--bustools", BUSTOOLS,
        "-t", "8",
        str(batch_file),
    ], check=True, cwd=work_dir)

    print("Done! Run Count function")
    elapsed = time.time() - start
    print(f"Took {elapsed:.1f} seconds")
    return out_dir


if __name__ == "__main__":
    # CLI test path — mirrors what gui.py will eventually call to handle any CLIs and specs
    parser = argparse.ArgumentParser()

    #species and parity 
    parser.add_argument("--species", default="homo_sapiens")
    parser.add_argument("--parity", default="paired", choices=["single", "paired"])

    #fastq and files 
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--fastq-dir")
    group.add_argument("--files", nargs="+") ##file paths as a whole one after another 

    #parser object gets created
    args = parser.parse_args()

    #calling ref builder
    index, t2g = ref_builder_cdna(species=args.species)

    ##paired vs unpaied files
    if args.parity == "paired":
        if args.fastq_dir is not None:
            paired,unpaired = find_file_pairs(Path(args.fastq_dir))
            run_count(paired, index, t2g, parity=args.parity)
        else:
            pairs,unpaired = pair_selected_files(args.files)
            run_count(pairs, index, t2g, parity=args.parity)

    if args.parity == "single":
        if args.fastq_dir is not None:
            singles = find_single_files(Path(args.fastq_dir))
            run_count(singles , index, t2g, parity=args.parity)
        else:
            singular = selected_files_as_singles(args.files)
            run_count(singular , index, t2g, parity=args.parity)


