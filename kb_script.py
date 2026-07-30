import subprocess
from pathlib import Path
import time
import shutil


def find_binary(name):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(
            f"{name} not found. Install with:\n"
            f"  conda install -c conda-forge -c bioconda {name}"
        )
    return path


KALLISTO = find_binary("kallisto")
BUSTOOLS = find_binary("bustools")

ALIASES = {
    "homo_sapiens": "human",
    "mus_musculus": "mouse",
    "canis_lupus_familiaris": "dog",
    "macaca_mulatta": "monkey",
    "danio_rerio": "zebrafish",
}

SUPPORTED_SPECIES = {"human", "mouse", "dog", "monkey", "zebrafish"}


def get_reference(species: str, work_dir: Path = Path("kb_work")) -> tuple[str, str]:
    """
    Checks for a locally cached prebuilt index; downloads via kb-python's
    -d flag (pachterlab/kallisto-transcriptome-indices) if missing.
    Fast — no genome build, just a direct download.
    """
    species = species.strip().lower()
    species = ALIASES.get(species, species)  # accepts either style of name

    if species not in SUPPORTED_SPECIES:
        raise ValueError(
            f"'{species}' not in prebuilt list: {sorted(SUPPORTED_SPECIES)}. "
            f"Use ref_builder() to build a custom index instead."
        )

    work_dir.mkdir(exist_ok=True)
    index_file = work_dir / f"{species}_index.idx"
    t2g_file = work_dir / f"{species}_t2g.txt"

    if index_file.exists() and t2g_file.exists():
        print(f"Using cached prebuilt index for {species}.")
        return str(index_file), str(t2g_file)

    print(f"Downloading prebuilt index for {species}")
    subprocess.run([
        "kb", "ref",
        "-d", species,
        "-i", str(index_file),
        "-g", str(t2g_file),
        "--kallisto", KALLISTO,
        "--bustools", BUSTOOLS,
    ], check=True)

    return str(index_file), str(t2g_file)


def ref_builder(species: str) -> tuple[str, str]:
    """
    Fallback for species NOT in the prebuilt list — builds from scratch via
    gget + kb ref. Slow (~20+ min). Only call this if get_reference() raises.
    """
    print("Hello! Starting Reference Transcript Build")
    start = time.time()

    species = species.strip().lower()
    index_file = f"{species}_index.idx"
    t2g_file = f"{species}_t2g.txt"
    cdna_file = f"{species}_cdna.fasta"

    work_dir = Path("kb_work")
    work_dir.mkdir(exist_ok=True)

    if (work_dir / index_file).exists() and (work_dir / t2g_file).exists():
        print(f"Reference for {species} already built — skipping.")
        return str(work_dir / index_file), str(work_dir / t2g_file)

    transcript_gget = subprocess.run(
        ["gget", "ref", "--ftp", "-w", "dna,gtf", species],
        capture_output=True, text=True, check=True,
    )
    gget_out = transcript_gget.stdout.split()

    shutil.rmtree(work_dir / "tmp", ignore_errors=True)
    subprocess.run([
        "kb", "ref",
        "-i", index_file,
        "-g", t2g_file,
        "-f1", cdna_file,
        "--kallisto", KALLISTO,
        "--bustools", BUSTOOLS,
        "-t", "4",
        *gget_out
    ], cwd=work_dir, check=True)

    elapsed = time.time() - start
    print("Done! Reference Transcript Build")
    print(f"Took {elapsed:.1f} seconds")

    return str(work_dir / index_file), str(work_dir / t2g_file)


##trying with cdna to minimize RAM usage
CDNA_URLS = {
    "human": "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    "mouse": "https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
}

def ref_builder_cdna(species: str) -> tuple[str, str]:
    """
    Fast path: builds index directly from a pre-made cDNA transcript FASTA
    (e.g. Ensembl's cdna.all.fa.gz), skipping genome+GTF splitting entirely.
    Uses --workflow=custom under kb ref. Currently human/mouse only.
    """
    print("Hello! Starting cDNA Reference Build")
    start = time.time()

    species = species.strip().lower()
    species = ALIASES.get(species, species)

    if species not in CDNA_URLS:
        raise ValueError(
            f"'{species}' not supported by ref_builder_cdna yet "
            f"(only {sorted(CDNA_URLS)}). Use ref_builder() instead for other species."
        )

    work_dir = Path("kb_work")
    work_dir.mkdir(exist_ok=True)
    index_file = f"{species}_index.idx"
    t2g_file = f"{species}_t2g.txt"

    if (work_dir / index_file).exists() and (work_dir / t2g_file).exists():
        print(f"Reference for {species} already built — skipping.")
        return str(work_dir / index_file), str(work_dir / t2g_file)

    cdna_fasta = work_dir / Path(CDNA_URLS[species]).name
    if not cdna_fasta.exists():
        print(f"Downloading cDNA FASTA for {species}...")
        subprocess.run(["curl", "-L", "-o", str(cdna_fasta), CDNA_URLS[species]], check=True)

    subprocess.run([
        "kb", "ref",
        "--workflow=custom",
        "-i", index_file,
        "-g", t2g_file,
        "--kallisto", KALLISTO,
        "-t" , "8",
        str(cdna_fasta.resolve()),
    ], cwd=work_dir, check=True)

    elapsed = time.time() - start
    print("Done! cDNA Reference Build")
    print(f"Took {elapsed:.1f} seconds")

    return str(work_dir / index_file), str(work_dir / t2g_file)



### FILE PAIR IDENTIFICATION AND DOWNSTREAM 
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
              work_dir: Path = Path(".")) -> Path:
    """
    parity: "paired" or "single" — samples dict shape must match (see write_batch_file).
    Returns the output directory path.
    """
    print("Hello! Beginning Run Count function")
    start = time.time()

    
    batch_file = write_batch_file(samples, work_dir / "batch.txt", parity)
    out_dir = work_dir / "counts_out"

    subprocess.run([
        "kb", "count",
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
    # CLI test path — mirrors what gui.py will eventually call
    index, t2g = get_reference("human")
    paired, unpaired = find_file_pairs(Path("testfq"))
    run_count(paired, index, t2g, parity="paired")