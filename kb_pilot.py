import sys
import subprocess
from pathlib import Path
import time
import gget
import argparse
import shutil
import gzip
import os
import concurrent.futures

def _bundle_root() -> Path:
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS)
    return Path(__file__).resolve().parent

BASE_DIR = _bundle_root()


def find_binary(name, base_dir):
    local_path = base_dir / "bin" / name
    if local_path.exists() and local_path.is_file():
        if not os.access(local_path, os.X_OK):
            raise RuntimeError(f"{local_path} exists but is not executable")
        return str(local_path)
    path = shutil.which(name)
    if path:
        return path
    raise RuntimeError(f"{name} not found in bin/ or PATH")


KALLISTO = find_binary("kallisto", BASE_DIR)
###BUSTOOLS = find_binary("bustools", BASE_DIR)

if getattr(sys, "frozen", False):
    APP_DATA_DIR = Path.home() / ".RAW2Compass" / "workdir"
else:
    APP_DATA_DIR = BASE_DIR / "workdir"

APP_DATA_DIR.mkdir(parents=True, exist_ok=True)


def _build_t2g_from_cdna(cdna_fasta: Path, t2g_out: Path) -> None:
    """
    Extracts transcript_id / gene_id / gene_name from Ensembl cDNA FASTA
    headers directly, since --workflow=custom doesn't generate t2g.txt. this is used for generation and matching

    file format is > header followed by
    Transcript ID | cDNA | Location | Gene Name | Gene Symbol | may have descriptions
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

def ref_builder_cdna(species: str) -> tuple[str, str]:
    """
    Fast path: builds index directly from a pre-made cDNA transcript FASTA
    (e.g. Ensembl's cdna.all.fa.gz), skipping genome+GTF splitting entirely.
    """
    print("Hello! Starting cDNA Reference Build")
    start = time.time()

    work_dir = APP_DATA_DIR / "kb_work"
    work_dir.mkdir(exist_ok=True)

    index_file = f"{species}_index.idx"
    t2g_file = f"{species}_t2g.txt"

    index_path = work_dir / index_file
    t2g_path = work_dir / t2g_file

    if index_path.exists() and t2g_path.exists() and t2g_path.stat().st_size > 0:
        print(f"Reference for {species} already built — skipping.")
        return str(index_path), str(t2g_path)

    result = subprocess.run(
        ["gget", "ref", "--ftp", "-w", "cdna", species],
        check=True,
        capture_output=True,
        text=True,
    )

    cdna_url = result.stdout.strip()
    if not cdna_url:
        raise RuntimeError(f"gget did not return a cDNA FASTA URL for {species}")
    cdna_fasta = work_dir / Path(cdna_url).name

    if not cdna_fasta.is_file():
        print(f"Downloading cDNA FASTA for {species}...")
        subprocess.run(["curl", "-L", "-o", str(cdna_fasta), cdna_url], check=True)

    if not cdna_fasta.is_file():
        raise RuntimeError(f"Expected a cDNA FASTA file, got: {cdna_fasta}")

    if not t2g_path.exists() or t2g_path.stat().st_size == 0:
        print("Generating transcript-to-gene mapping from FASTA headers...")
        _build_t2g_from_cdna(cdna_fasta, t2g_path)

    print("kallisto index running...")
    subprocess.run(
        [KALLISTO, "index", "-i", str(index_path), str(cdna_fasta.resolve())],
        cwd=work_dir,
        check=True,
    )

    elapsed = time.time() - start
    print("Done! cDNA Reference Build")
    print(f"Took {elapsed:.1f} seconds")

    return str(index_path), str(t2g_path)

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
                rev_name = file.name.replace(fwd_token, rev_token)
                rev = file.parent / rev_name
                if rev.exists():
                    sample = file.name[:file.name.find(fwd_token)]
                    paired[sample] = (file, rev)
                    matched = True
                break

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

## KALLISTO RUN HELPER ###
def _run_kallisto_sample(sample_name: str, fwd: str, rev: str | None, index: str, out_dir: Path, parity: str, threads: int) -> str:
    sample_outdir = out_dir / sample_name
    sample_outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        KALLISTO,
        "quant",
        "-i", index,
        "-o", str(sample_outdir),
        "-t", str(threads),
    ]

    if parity == "single":
        cmd.extend(["--single", fwd])
    else:
        cmd.extend([fwd, rev])

    print(f"[{sample_name}] starting kallisto quant")
    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    ) as proc:
        if proc.stdout is None:
            raise RuntimeError(f"No stdout captured for {sample_name}")

        for line in iter(proc.stdout.readline, ""):
            if line:
                print(f"[{sample_name}] {line.rstrip()}")

        returncode = proc.wait()
        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, cmd)

    print(f"[{sample_name}] completed kallisto quant")
    return sample_name


#####

def run_count(samples: dict,
              index: str,
              t2g: str,
              parity: str = "paired",
              work_dir: Path = None,
              parallel: bool = True,
              max_workers: int = None,
              threads_per_sample: int = 8) -> Path:
    """
    parity: "paired" or "single"
    If parallel=True, each sample runs in a separate worker.
    """
    print("Hello! Beginning Run Count function")

    if work_dir is None:
        work_dir = APP_DATA_DIR / "counts_out_workdir"
    work_dir.mkdir(parents=True, exist_ok=True)

    start = time.time()

    batch_file = write_batch_file(samples, work_dir / "batch.txt", parity)
    out_dir = work_dir / "counts_out"
    out_dir.mkdir(parents=True, exist_ok=True)

    if max_workers is None:
        available = max(1, os.cpu_count() or 1)
        max_workers = min(len(samples), max(1, available // 2)) ## desinate maxima using cpu presence

    if parallel:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for sample_name, sample_data in samples.items():
                if parity == "paired":
                    fwd, rev = sample_data
                    futures.append(
                        executor.submit(
                            _run_kallisto_sample, ## calling the helper to execute subprocess cmd 
                            sample_name,
                            str(fwd),
                            str(rev),
                            index,
                            out_dir,
                            parity,
                            threads_per_sample,
                        )
                    )
                else:
                    fwd = sample_data
                    futures.append(
                        executor.submit(
                            _run_kallisto_sample,
                            sample_name,
                            str(fwd),
                            None,
                            index,
                            out_dir,
                            parity,
                            threads_per_sample,
                        )
                    )

            for future in concurrent.futures.as_completed(futures):
                sample_name = future.result()
                print(f"Finished {sample_name}")
    else:
        for sample_name, sample_data in samples.items():
            if parity == "paired":
                fwd, rev = sample_data
                _run_kallisto_sample(
                    sample_name,
                    str(fwd),
                    str(rev),
                    index,
                    out_dir,
                    parity,
                    threads_per_sample,
                )
            else:
                fwd = sample_data
                _run_kallisto_sample(
                    sample_name,
                    str(fwd),
                    None,
                    index,
                    out_dir,
                    parity,
                    threads_per_sample,
                )

    print("Done! Kallisto Quant function")
    elapsed = time.time() - start
    print(f"Took {elapsed:.1f} seconds")
    return out_dir

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--species", default="homo_sapiens")
    parser.add_argument("--parity", default="paired", choices=["single", "paired"])
    parser.add_argument("--fastq-dir")
    parser.add_argument("--files", nargs="+")
    parser.add_argument("--max-workers", type=int, default=None)
    parser.add_argument("--threads-per-sample", type=int, default=8)
    parser.add_argument("--parallel", dest="parallel", action="store_true", default=True)
    parser.add_argument("--no-parallel", dest="parallel", action="store_false")

    args = parser.parse_args()

    index, t2g = ref_builder_cdna(species=args.species)

    if args.parity == "paired":
        if args.fastq_dir is not None:
            paired, unpaired = find_file_pairs(Path(args.fastq_dir))
            run_count(
                paired,
                index,
                t2g,
                parity=args.parity,
                parallel=args.parallel,
                max_workers=args.max_workers,
                threads_per_sample=args.threads_per_sample,
            )
        else:
            pairs, unpaired = pair_selected_files(args.files)
            run_count(
                pairs,
                index,
                t2g,
                parity=args.parity,
                parallel=args.parallel,
                max_workers=args.max_workers,
                threads_per_sample=args.threads_per_sample,
            )

    if args.parity == "single":
        if args.fastq_dir is not None:
            singles = find_single_files(Path(args.fastq_dir))
            run_count(
                singles,
                index,
                t2g,
                parity=args.parity,
                parallel=args.parallel,
                max_workers=args.max_workers,
                threads_per_sample=args.threads_per_sample,
            )
        else:
            singular = selected_files_as_singles(args.files)
            run_count(
                singular,
                index,
                t2g,
                parity=args.parity,
                parallel=args.parallel,
                max_workers=args.max_workers,
                threads_per_sample=args.threads_per_sample,
            )