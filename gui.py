from pathlib import Path

from PySide6.QtCore import QThread, Signal, Qt
from PySide6.QtGui import QTextCursor
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QComboBox,
    QLineEdit,
    QPushButton,
    QFileDialog,
    QListWidget,
    QListWidgetItem,
    QTextEdit,
    QMessageBox,
)

from kb_pilot import (
    ref_builder_cdna,
    pair_selected_files,
    selected_files_as_singles,
    run_count,
)
from matrix import load_kb_counts, normalize_cpm, log_transform, export_matrix

NAVY = "#0F2B46"
ACCENT = "#2E5FA3"
ACCENT_HOVER = "#25507F"
BORDER = "#DDE5EC"
MUTED = "#5B6B7B"
WHITE = "#FFFFFF"

STYLESHEET = f"""
QMainWindow {{ background-color: {WHITE}; color: {NAVY}; }}
QWidget {{ background-color: {WHITE}; color: {NAVY}; }}
QLabel {{ color: {MUTED}; }}
QLabel#header {{ color: {NAVY}; font-size: 26px; font-weight: 700; }}
QLabel#eyebrow {{ color: {ACCENT}; font-size: 11px; font-weight: 700; letter-spacing: 1px; }}
QPushButton {{
    background-color: {ACCENT};
    color: white;
    border: none;
    border-radius: 12px;
    padding: 10px 18px;
    font-weight: 700;
}}
QPushButton:hover {{ background-color: {ACCENT_HOVER}; }}
QPushButton:disabled {{
    background-color: {BORDER};
    color: {MUTED};
    border: none;
}}
QLineEdit, QComboBox, QTextEdit, QListWidget {{
    background-color: white;
    color: {NAVY};
    border: 1px solid {BORDER};
    border-radius: 8px;
    padding: 8px 10px;
}}
QListWidget::item:selected {{
    background-color: {ACCENT};
    color: white;
}}
"""

class PipelineWorker(QThread):
    progress = Signal(str)
    output = Signal(str)
    finished_ok = Signal(object)
    failed = Signal(str)

    def __init__(self, species, samples, parity, max_workers=None, threads_per_sample=8):
        super().__init__()
        self.species = species
        self.samples = samples
        self.parity = parity
        self.max_workers = max_workers
        self.threads_per_sample = threads_per_sample

    def run(self):
        try:
            self.progress.emit("Building reference...")
            index, t2g = ref_builder_cdna(self.species, output=self.output.emit)

            self.progress.emit("Running kallisto quant...")
            out_dir = run_count(
                self.samples,
                index,
                t2g,
                parity=self.parity,
                parallel=True,
                max_workers=self.max_workers,
                threads_per_sample=self.threads_per_sample,
                output=self.output.emit,
            )

            self.finished_ok.emit(out_dir)
        except Exception as e:
            self.failed.emit(str(e))

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RAW2Compass")
        self.resize(900, 650)
        self.detected_samples = {}

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(12)

        eyebrow = QLabel("FASTQ TO COUNTS / MATRIX")
        eyebrow.setObjectName("eyebrow")
        layout.addWidget(eyebrow)

        header = QLabel("RAW2Compass")
        header.setObjectName("header")
        layout.addWidget(header)

        layout.addWidget(QLabel("Species"))
        self.species_combo = QComboBox()
        self.species_combo.addItem("Homo sapiens", "homo_sapiens")
        self.species_combo.addItem("Mus musculus", "mus_musculus")
        self.species_combo.addItem("Canis lupus familiaris", "canis_lupus_familiaris")
        self.species_combo.addItem("Macaca mulatta", "macaca_mulatta")
        self.species_combo.addItem("Danio rerio", "danio_rerio")
        layout.addWidget(self.species_combo)

        layout.addWidget(QLabel("FASTQ files"))
        row = QHBoxLayout()
        self.file_box = QLineEdit()
        self.file_box.setReadOnly(True)
        self.file_box.setPlaceholderText("No files selected")
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.pick_files)
        row.addWidget(self.file_box)
        row.addWidget(browse_btn)
        layout.addLayout(row)

        layout.addWidget(QLabel("Read type"))
        self.parity_combo = QComboBox()
        self.parity_combo.addItems(["paired", "single"])
        layout.addWidget(self.parity_combo)

        layout.addWidget(QLabel("Samples"))
        self.sample_list = QListWidget()
        layout.addWidget(self.sample_list)

        self.run_btn = QPushButton("Run Pipeline")
        self.run_btn.clicked.connect(self.run_pipeline)
        layout.addWidget(self.run_btn)

        self.log = QTextEdit()
        self.log.setReadOnly(True)
        layout.addWidget(self.log)

        self.setCentralWidget(central)

    def pick_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Select FASTQ files",
            "",
            "FASTQ files (*.fastq *.fastq.gz *.fq *.fq.gz)",
        )
        if not files:
            return

        self.file_box.setText(f"{len(files)} file(s) selected")

        file_paths = [Path(p) for p in files]
        parity = self.parity_combo.currentText()

        if parity == "paired":
            samples, unmatched = pair_selected_files(file_paths)
            if unmatched:
                self.log.append(f"Skipped {len(unmatched)} unpaired files")
        else:
            samples = selected_files_as_singles(file_paths)

        self.detected_samples = samples
        self.sample_list.clear()

        for sample_name in samples:
            item = QListWidgetItem(sample_name)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked)
            self.sample_list.addItem(item)

        self.log.append(f"Detected {len(samples)} sample(s)")

    def run_pipeline(self):
        selected = {}
        for i in range(self.sample_list.count()):
            item = self.sample_list.item(i)
            if item.checkState() == Qt.Checked:
                selected[item.text()] = self.detected_samples[item.text()]

        if not selected:
            QMessageBox.warning(self, "No samples selected", "Select at least one sample before running.")
            return

        self.run_btn.setEnabled(False)
        self.log.append("Starting pipeline...")

        self.worker = PipelineWorker(
            species=self.species_combo.currentData(),
            samples=selected,
            parity=self.parity_combo.currentText(),
            max_workers=None,
            threads_per_sample=4,
        )
        self.worker.progress.connect(self.log.append)
        self.worker.output.connect(self.append_pipeline_output)
        self.worker.finished_ok.connect(self.on_pipeline_done)
        self.worker.failed.connect(self.on_pipeline_error)
        self.worker.start()

    def append_pipeline_output(self, text):
        self.log.moveCursor(QTextCursor.End)
        self.log.insertPlainText(f"{text}\n")
        self.log.ensureCursorVisible()

    def on_pipeline_done(self, out_dir):
        self.log.append(f"Pipeline complete. Output: {out_dir}")

        try:
            counts = load_kb_counts(out_dir)
            cpm = normalize_cpm(counts)
            log_cpm = log_transform(cpm)

            export_matrix(cpm, out_dir.parent / "cpm_matrix.csv")
            export_matrix(log_cpm, out_dir.parent / "log_cpm_matrix.csv")

            self.log.append(f"Exported CPM and log-CPM matrices to {out_dir.parent}")
        except Exception as e:
            self.log.append(f"Matrix export error: {e}")

        self.run_btn.setEnabled(True)

    def on_pipeline_error(self, message):
        self.log.append(f"Error: {message}")
        self.run_btn.setEnabled(True)

if __name__ == "__main__":
    app = QApplication([])
    app.setStyleSheet(STYLESHEET)
    window = MainWindow()
    window.show()
    app.exec()