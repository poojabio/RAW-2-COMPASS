import subprocess
from pathlib import Path
import time
import shutil
#import pandas as pd

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

''' GENERAL IDEA
1) kb-python for kallisto-based processing 
2) pandas for normalization 
3) PyInstaller 
4) PyFreeze 
'''
def ref_builder():
    print("Hello! Starting Reference Transcript Build")
    start = time.time()

    ### reference 
    species = input("Species (e.g. homo_sapiens, mouse, dog, monkey, or zebrafish etc. ): ").strip().lower()

    ALIASES = {
    "human":     "homo_sapiens",
    "mouse":     "mus_musculus",
    "dog":       "canis_lupus_familiaris",
    "monkey":    "macaca_mulatta",
    "zebrafish": "danio_rerio",}

    species = species.strip().lower()
    species = ALIASES.get(species, species)


    ### index (use subprocess to sub out the initial bsh approach) 
    """
    Running kb ref --workflow=standard will generate three files:
    index.idx: Contains the kallisto index used for pseudoalignment and quantification.
    t2g.txt: A transcript-to-gene mapping file, linking each transcript in the index to its corresponding gene.
    cdna.fasta: A FASTA file containing the transcript sequences extracted from the input genome FASTA and GTF. This file is not required in downstream steps, but is useful to keep as a reference.
    """

    index_file = f"{species}_index.idx"
    t2g_file = f"{species}_t2g.txt"
    cdna_file = f"{species}_cdna.fasta"

    transcript_gget = subprocess.run(["gget","ref","--ftp","-w","dna,gtf",species], 
                                     capture_output=True, 
                                     text=True,
                                     check=True)
    
    gget_out = transcript_gget.stdout.split()

    work_dir = Path("kb_work")
    work_dir.mkdir(exist_ok=True)

    if (work_dir / index_file).exists() and (work_dir / t2g_file).exists():
        print(f"Reference for {species} already built — skipping.")
        return str(work_dir / index_file), str(work_dir / t2g_file)

    shutil.rmtree(work_dir / "tmp", ignore_errors=True) ## clean out anythign else impeding 
    subprocess.run([
        "kb", "ref",
        "-i", index_file,
        "-g", t2g_file,
        "-f1", cdna_file,
        "--kallisto", KALLISTO,
        "--bustools", BUSTOOLS,
        "-t", "8",
        *gget_out
    ], cwd=work_dir, check=True)

    elapsed = time.time() - start
    print( "Done! Reference Transcript Build")
    print(f"Took {elapsed:.1f} seconds")

    return str(work_dir / index_file), str(work_dir / t2g_file)

###ref_builder()

#finding pairs of fwd and rev reads
def find_file_pairs():
    print("Parsing through for file pairs")
    start = time.time()


    paired = {}
    unpaired = []
    folder = Path(input("Path to FASTQ folder: ").strip()).expanduser()

# exit or re-prompt for folder path
    while True:
        if  folder.is_dir():
            break
        print(f"Not a folder: {folder}")
        folder = Path(input("Path to FASTQ folder: ").strip()).expanduser()

    PAIR_PATTERNS =[
    ("_R1", "_R2"),
    ("read1", "read2"),
    ("forward", "reverse"),
    ("_1", "_2")]

    for file in folder.iterdir():
    #if file.suffix == ".gz":
        ###sample = file.stem
        for fwd_token, rev_token in PAIR_PATTERNS:
            if fwd_token in file.name:
                fwd = file
                rev_name = file.name.replace(fwd_token, rev_token) ## TEMPORARILY THE FILE LOOKING TO BE FOUND
                rev = file.parent / rev_name ##reverse file as a whole

                if rev.exists():
                    sample = file.name[:file.name.find(fwd_token)] # this takes the sample name essentially like .stem however considering naming conventions this is more optimal
                    paired[sample] = (file, rev)

                else:
                    unpaired.append(file)

                    
        else:
            unpaired.append(file)
                    ###paired[fwd] = rev

    for sample, (fwd, rev) in paired.items():
        print(f"{sample}: {fwd.name} + {rev.name}")
        #print(f"Oops we had some unpaired {unpaired}")


    print( "Done! Pair finding")

    elapsed = time.time() - start
    print(f"Took {elapsed:.1f} seconds")
    return paired, unpaired

def write_batch_file(paired, out_path: Path) -> Path:
    """
    paired: {sample_id: (fwd_path, rev_path)}
    Writes a tab-separated batch file kb count expects for multi-sample BULK runs:
        sample_id    fastq_1    fastq_2
    """
    with open(out_path, "w") as f:
        for sample_id, (fwd, rev) in paired.items():
            f.write(f"{sample_id}\t{fwd}\t{rev}\n")
    return out_path


def run_count(paired, index, t2g,  work_dir=Path(".")):
    print("Hello! Beginning Run Count function")
    start = time.time()

    batch_file = write_batch_file(paired, work_dir / "batch.txt")

    subprocess.run([
    "kb", "count",
    "-i", index,
    "-g", t2g,
    "-x", "BULK",
    "--parity", "paired",
    "-o", "counts_out",
    "--kallisto", KALLISTO,
    "--bustools", BUSTOOLS,
    "-t", "8",
    str(batch_file),
], check=True, cwd=work_dir)

    print("Done! Run Count function")
    elapsed = time.time() - start
    print(f"Took {elapsed:.1f} seconds")


            
index, t2g = ref_builder()
paired, unpaired = find_file_pairs()
run_count(paired, index, t2g)


