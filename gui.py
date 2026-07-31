import sys
from pathlib import Path
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QLabel, QLineEdit, QPushButton, QComboBox,
    QFileDialog, QTextEdit, QListWidget, QListWidgetItem,
    QTableWidget, QTableWidgetItem, QMessageBox, QCheckBox
)
from PySide6.QtCore import QThread, Signal, Qt

from kb_script import ref_builder_cdna, pair_selected_files, selected_files_as_singles, run_count, APP_DATA_DIR
from matrix import load_kb_counts, normalize_cpm, log_transform, export_matrix

NAVY = "#0F2B46"
ACCENT_BLUE = "#2E5FA3"
ACCENT_BLUE_HOVER = "#25507F"
BODY_GRAY = "#5B6B7B"
BORDER = "#DDE5EC"
TEAL = "#1F9D82"

STYLESHEET = f"""
QMainWindow {{ background-color: white; }}
QWidget {{ background-color: white; }}

QLabel {{ color: {BODY_GRAY}; font-size: 13px; }}
QLabel#header {{ color: {NAVY}; font-size: 26px; font-weight: 800; }}
QLabel#eyebrow {{ color: {ACCENT_BLUE}; font-size: 11px; font-weight: 700; letter-spacing: 1px; }}

QTabWidget::pane {{ border: none; background-color: white; }}
QTabBar::tab {{ background: transparent; color: {BODY_GRAY}; padding: 10px 18px; font-weight: 600; font-size: 13px; }}
QTabBar::tab:selected {{ color: {NAVY}; border-bottom: 2px solid {ACCENT_BLUE}; }}

QLineEdit, QComboBox {{
    background-color: white; border: 1px solid {BORDER}; border-radius: 8px;
    padding: 8px 10px; color: {NAVY}; font-size: 13px;
}}
QLineEdit:focus, QComboBox:focus {{ border: 1px solid {ACCENT_BLUE}; }}

QListWidget, QTableWidget {{
    background-color: white; border: 1px solid {BORDER}; border-radius: 8px;
    color: {NAVY}; font-size: 13px;
}}
QListWidget::item, QTableWidget::item {{ padding: 6px; color: {NAVY}; }}
QListWidget::item:selected, QTableWidget::item:selected {{ background-color: {ACCENT_BLUE}; color: white; }}

QHeaderView::section {{
    background-color: white; color: {BODY_GRAY}; border: none;
    border-bottom: 1px solid {BORDER}; padding: 6px; font-weight: 600;
}}

QPushButton {{
    background-color: {ACCENT_BLUE}; color: white; border: none; border-radius: 20px;
    padding: 10px 22px; font-weight: 700; font-size: 13px;
}}
QPushButton:hover {{ background-color: {ACCENT_BLUE_HOVER}; }}
QPushButton:disabled {{ background-color: {BORDER}; color: {BODY_GRAY}; }}
QPushButton#secondary {{ background-color: white; color: {ACCENT_BLUE}; border: 1px solid {ACCENT_BLUE}; }}
QPushButton#secondary:hover {{ background-color: #EAF1F9; }}

QTextEdit {{
    background-color: white; border: 1px solid {BORDER}; border-radius: 8px;
    color: {NAVY}; font-family: Menlo, Consolas, monospace; font-size: 12px; padding: 8px;
}}
"""


class PipelineWorker(QThread):
    progress = Signal(str)
    finished_ok = Signal(object)
    failed = Signal(str)

    def __init__(self, species, samples, parity):
        super().__init__()
        self.species = species
        self.samples = samples  # already-selected/paired dict
        self.parity = parity

    def run(self):
        try:
            self.progress.emit("Fetching reference index...")
            index, t2g = ref_builder_cdna(self.species)

            self.progress.emit("Running kb count (this may take a while)...")
            out_dir = run_count(self.samples, index, t2g, parity=self.parity)

            self.finished_ok.emit(out_dir)
        except Exception as e:
            self.failed.emit(str(e))


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RAW2Compass")
        self.resize(960, 660)

        self.detected_samples = {}

        tabs = QTabWidget()
        tabs.addTab(self._build_run_panel(), "Run Pipeline")
        ##tabs.addTab(self._build_batch_panel(), "Batch Correction")
        tabs.addTab(self._build_results_panel(), "Results")
        self.setCentralWidget(tabs)

    # ---------- Run Pipeline tab ----------

    def _build_run_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(30, 30, 30, 30)
        layout.setSpacing(14)

        eyebrow = QLabel("FASTQ TO COMPASS-READY MATRIX")
        eyebrow.setObjectName("eyebrow")
        layout.addWidget(eyebrow)

        header = QLabel("RAW2Compass")
        header.setObjectName("header")
        layout.addWidget(header)
        layout.addSpacing(10)

        layout.addWidget(QLabel("Species"))
        self.species_dropdown = QComboBox()
        self.species_dropdown.addItems(["human", "mouse", "dog", "monkey", "zebrafish"])
        layout.addWidget(self.species_dropdown)

        layout.addWidget(QLabel("FASTQ Files"))
        file_row = QHBoxLayout()
        self.folder_input = QLineEdit()
        self.folder_input.setReadOnly(True)
        self.folder_input.setPlaceholderText("No files selected yet")
        browse_btn = QPushButton("Select FASTQ Files")
        browse_btn.setObjectName("secondary")
        browse_btn.clicked.connect(self._pick_files)
        file_row.addWidget(self.folder_input)
        file_row.addWidget(browse_btn)
        layout.addLayout(file_row)

        layout.addWidget(QLabel("Read Type"))
        self.parity_dropdown = QComboBox()
        self.parity_dropdown.addItems(["paired", "single"])
        layout.addWidget(self.parity_dropdown)

        select_row = QHBoxLayout()
        select_row.addWidget(QLabel("Samples — uncheck any to exclude from this run"))
        select_all_btn = QPushButton("Select All")
        select_all_btn.setObjectName("secondary")
        select_all_btn.clicked.connect(lambda: self._set_all_checked(True))
        none_btn = QPushButton("None")
        none_btn.setObjectName("secondary")
        none_btn.clicked.connect(lambda: self._set_all_checked(False))
        select_row.addWidget(select_all_btn)
        select_row.addWidget(none_btn)
        layout.addLayout(select_row)

        self.sample_list = QListWidget()
        layout.addWidget(self.sample_list)

        layout.addSpacing(6)
        self.run_btn = QPushButton("Run Pipeline  →")
        self.run_btn.clicked.connect(self._start_pipeline)
        layout.addWidget(self.run_btn, alignment=Qt.AlignLeft)

        layout.addSpacing(10)
        layout.addWidget(QLabel("Status"))
        self.status_log = QTextEdit()
        self.status_log.setReadOnly(True)
        layout.addWidget(self.status_log)

        panel.setLayout(layout)
        return panel

    def _pick_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select FASTQ files", "",
            "FASTQ files (*.fastq *.fastq.gz *.fq *.fq.gz)"
        )
        if not files:
            return

        file_paths = [Path(f) for f in files]
        self.folder_input.setText(f"{len(file_paths)} file(s) selected")
        self._detect_samples_from_files(file_paths)

    def _detect_samples_from_files(self, file_paths: list[Path]):
        self.status_log.append(f"Processing {len(file_paths)} selected file(s)...")
        try:
            parity = self.parity_dropdown.currentText()
            if parity == "paired":
                samples, unmatched = pair_selected_files(file_paths)
                if unmatched:
                    self.status_log.append(
                        f"Note: {len(unmatched)} file(s) didn't pair up and are excluded: "
                        f"{[f.name for f in unmatched]}"
                    )
            else:
                samples = selected_files_as_singles(file_paths)

            self.detected_samples = samples
            self.sample_list.clear()
            for sample_id in samples:
                item = QListWidgetItem(sample_id)
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setCheckState(Qt.Checked)
                self.sample_list.addItem(item)

            self.status_log.append(f"Found {len(samples)} sample(s).")
        except Exception as e:
            self.status_log.append(f"Error processing files: {e}")

    def _set_all_checked(self, checked: bool):
        state = Qt.Checked if checked else Qt.Unchecked
        for i in range(self.sample_list.count()):
            self.sample_list.item(i).setCheckState(state)

    def _get_selected_samples(self) -> dict:
        selected_ids = set()
        for i in range(self.sample_list.count()):
            item = self.sample_list.item(i)
            if item.checkState() == Qt.Checked:
                selected_ids.add(item.text())
        return {sid: path for sid, path in self.detected_samples.items() if sid in selected_ids}

    def _start_pipeline(self):
        selected = self._get_selected_samples()
        if not selected:
            QMessageBox.warning(self, "No samples selected",
                                 "Select at least one sample before running.")
            return

        self.run_btn.setEnabled(False)
        self.status_log.append(f"Starting pipeline on {len(selected)} sample(s)...")

        self.worker = PipelineWorker(
            species=self.species_dropdown.currentText(),
            samples=selected,
            parity=self.parity_dropdown.currentText(),
        )
        self.worker.progress.connect(self.status_log.append)
        self.worker.finished_ok.connect(self._on_pipeline_done)
        self.worker.failed.connect(self._on_pipeline_error)
        self.worker.start()

    def _on_pipeline_done(self, out_dir):
        self.status_log.append(f"Done! Output: {out_dir}")
        self.run_btn.setEnabled(True)

    def _on_pipeline_error(self, error_msg):
        self.status_log.append(f"Error: {error_msg}")
        self.run_btn.setEnabled(True)

    def _build_results_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(30, 30, 30, 30)
        layout.setSpacing(14)

        eyebrow = QLabel("NORMALIZATION & EXPORT")
        eyebrow.setObjectName("eyebrow")
        layout.addWidget(eyebrow)

        header = QLabel("Results")
        header.setObjectName("header")
        layout.addWidget(header)

        self.load_results_btn = QPushButton("Load Counts from Last Run")
        self.load_results_btn.setObjectName("secondary")
        self.load_results_btn.clicked.connect(self._load_results)
        layout.addWidget(self.load_results_btn, alignment=Qt.AlignLeft)

        # --- NEW: preview table goes here ---
        layout.addWidget(QLabel("Preview (first 20 genes)"))
        self.preview_table = QTableWidget(0, 0)
        self.preview_table.setMaximumHeight(300)
        layout.addWidget(self.preview_table)
        # --- end new ---

        self.export_btn = QPushButton("Normalize + Export  →")
        self.export_btn.clicked.connect(self._run_normalization)
        layout.addWidget(self.export_btn, alignment=Qt.AlignLeft)

        self.results_log = QTextEdit()
        self.results_log.setReadOnly(True)
        layout.addWidget(self.results_log)

        panel.setLayout(layout)
        return panel

    def _toggle_batch_table(self, state):
        checked = state == Qt.Checked.value
        self.results_batch_table.setVisible(checked)
        if checked:
            self.results_batch_table.setRowCount(0)
            for sample_id in self.detected_samples:
                row = self.results_batch_table.rowCount()
                self.results_batch_table.insertRow(row)
                self.results_batch_table.setItem(row, 0, QTableWidgetItem(sample_id))
                self.results_batch_table.setItem(row, 1, QTableWidgetItem(""))

    def _load_results(self):
        try:
            counts_dir = APP_DATA_DIR / "counts_out_workdir" / "counts_out" / "counts_unfiltered"
            counts = load_kb_counts(counts_dir)
            self._loaded_counts = counts

            n_genes, n_samples = counts.shape
            n_nonzero = (counts.values != 0).sum()
            pct_nonzero = 100 * n_nonzero / counts.size if counts.size else 0
            top_genes = counts.sum(axis=1).sort_values(ascending=False).head(5)

            self.results_log.clear()
            self.results_log.append(f"Loaded {n_genes} genes x {n_samples} sample(s)")
            self.results_log.append(f"Non-zero entries: {pct_nonzero:.1f}%")
            self.results_log.append("Top expressed genes:")
            for gene, val in top_genes.items():
                self.results_log.append(f"  {gene}: {val:.1f}")

            self._populate_preview_table(counts.head(20))
        except Exception as e:
            self.results_log.append(f"Error loading counts: {e}")

    def _populate_preview_table(self, df):
        self.preview_table.setRowCount(df.shape[0])
        self.preview_table.setColumnCount(df.shape[1])
        self.preview_table.setHorizontalHeaderLabels([str(c) for c in df.columns])
        self.preview_table.setVerticalHeaderLabels([str(i) for i in df.index])
        for i, (_, row) in enumerate(df.iterrows()):
            for j, val in enumerate(row):
                self.preview_table.setItem(i, j, QTableWidgetItem(f"{val:.1f}"))

    def _run_normalization(self):
        if not hasattr(self, "_loaded_counts"):
            QMessageBox.warning(self, "No data", "Load counts before normalizing.")
            return

        try:
            cpm = normalize_cpm(self._loaded_counts)
            log_cpm = log_transform(cpm)

            export_matrix(cpm, APP_DATA_DIR / "cpm_matrix.csv")
            export_matrix(log_cpm, APP_DATA_DIR / "log_cpm_matrix.csv")
            self.results_log.append(f"Exported to {APP_DATA_DIR}")

        except Exception as e:
            self.results_log.append(f"Error: {e}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyleSheet(STYLESHEET)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())