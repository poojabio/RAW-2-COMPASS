# RAW-2-COMPASS
### RAW-2-COMPASS is a desktop application for processing FASTQ sequencing reads into gene-level count tables and normalized expression matrices. It builds a transcript reference with kallisto, quantifies selected samples, and exports CPM and log-CPM matrices.

## Ways to use the app

- **Prebuilt macOS app:** If you have a prebuilt app bundle, open it. A separate Python environment is not needed to run the bundled app.
- **Run from source:** Downloadable source code is available from [GitHub](https://github.com/poojabio/RAW2Compass). Download or clone it, then create a local Python environment as described below.

## Processing workflow

1. **Choose a reference species.** The GUI offers human (*Homo sapiens*), mouse (*Mus musculus*), dog (*Canis lupus familiaris*), rhesus macaque (*Macaca mulatta*), and zebrafish (*Danio rerio*).
2. **Select FASTQ files and read type.** The application supports paired-end and single-end reads. For paired-end files, it looks for matching names containing `_R1`/`_R2`, `read1`/`read2`, `forward`/`reverse`, or `_1`/`_2`. Unmatched files are skipped. In single-end mode, each selected file is treated as a separate sample.
3. **Review samples.** Detected samples appear in a checklist; all are selected by default, and individual samples can be unchecked.
4. **Build or reuse the transcript reference.** The app looks up a cDNA FASTA URL with `gget`, downloads the FASTA if it is not already present, and runs kallisto to create a transcript index. It also parses transcript, gene, and gene-symbol identifiers from FASTA headers to create a transcript-to-gene mapping. An existing index and nonempty mapping for the species are reused. The latest reference is pulled from the ENSEMBL database if a prior index is not provided/synthesized by the user from earlier runs
5. **Quantify samples.** kallisto quantifies each selected sample against the reference. The GUI runs samples concurrently, assigning four kallisto threads per sample; progress and tool output are shown in the app.
6. **Create expression matrices.** The app reads each sample's `abundance.tsv`, joins transcript estimates to gene symbols, sums transcript estimated counts by gene symbol, and calculates:
   - **CPM:** counts per million, calculated per sample.
   - **log-CPM:** `log2(CPM + 1)`.

## Outputs and working files

For a source run, intermediate reference files and quantification results are kept under the project directory's `workdir/`. When running the frozen app, working files are stored under `~/.RAW2Compass/workdir/`.

The quantification results are organized by sample and include kallisto's `abundance.tsv`. The GUI exports these two files alongside the `counts_out` directory:

- `cpm_matrix.csv`
- `log_cpm_matrix.csv`

Matrix rows are gene symbols and columns are samples. Reference FASTA files, indexes, and quantification outputs can be large; make sure there is sufficient free disk space and an internet connection for the initial reference download.

## Run from source in a local environment

Python 3.10 or newer is required. From the downloaded or cloned project directory, create and activate a virtual environment, install the listed dependencies, and launch the GUI:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python gui.py
```

On Windows, activate the environment with `.venv\Scripts\activate` instead of the `source` command. If you prefer Conda, create and activate an environment before installing the same requirements:

```bash
conda create --prefix .conda python=3.10
conda activate ./.conda
python -m pip install -r requirements.txt
python gui.py
```

The project includes a `bin/kallisto` executable used for indexing and quantification. The application also needs `curl` on `PATH` to download reference FASTA files. Alternatively, `kb_pilot.py` exposes the reference-building and quantification workflow as a command-line entry point:

```bash
python kb_pilot.py --species homo_sapiens --parity paired --fastq-dir /path/to/fastq_files
```

Use `--parity single` for single-end reads. Run `python kb_pilot.py --help` for additional options, including worker and thread settings.

## Main code modules

- **`gui.py`** — PySide6 desktop interface, FASTQ selection and sample checklist, background pipeline worker, progress display, and automatic matrix export.
- **`kb_pilot.py`** — reference download and kallisto index creation, FASTQ pairing, batch creation, per-sample quantification, parallel execution, and command-line entry point.
- **`matrix.py`** — loading and merging kallisto abundance files, mapping transcript estimates to gene symbols, CPM normalization, log transformation, and CSV export.

## Notes

- The first run for a species may take longer because it downloads and indexes the cDNA reference. Later runs reuse the generated reference files when available.
- The app reports processing progress and errors in its log area. Keep the selected FASTQ files accessible until quantification finishes.
