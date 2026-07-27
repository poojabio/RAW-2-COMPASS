import sys
from pathlib import Path
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QLabel, QLineEdit, QPushButton, QComboBox,
    QFileDialog, QTextEdit
)
from PySide6.QtCore import QThread, Signal, Qt

from kb_script import get_reference, find_file_pairs, find_single_files, run_count
from matrix import load_kb_counts, normalize_cpm, log_transform, export_matrix


# Colors pulled from the COMPASS Prep site
NAVY = "#0F2B46"
ACCENT_BLUE = "#2E5FA3"
ACCENT_BLUE_HOVER = "#25507F"
BODY_GRAY = "#5B6B7B"
BG = "#F7FAFC"
BORDER = "#DDE5EC"
TEAL = "#1F9D82"

STYLESHEET = f"""
QMainWindow {{
    background-color: {BG};
}}

QLabel {{
    color: {BODY_GRAY};
    font-size: 13px;
}}

QLabel#header {{
    color: {NAVY};
    font-size: 26px;
    font-weight: 800;
}}

QLabel#eyebrow {{
    color: {ACCENT_BLUE};
    font-size: 11px;
    font-weight: 700;
    letter-spacing: 1px;
}}

QTabWidget::pane {{
    border: none;
    background-color: {BG};
}}

QTabBar::tab {{
    background: transparent;
    color: {BODY_GRAY};
    padding: 10px 18px;
    font-weight: 600;
    font-size: 13px;
}}

QTabBar::tab:selected {{
    color: {NAVY};
    border-bottom: 2px solid {ACCENT_BLUE};
}}

QLineEdit, QComboBox {{
    background-color: white;
    border: 1px solid {BORDER};
    border-radius: 8px;
    padding: 8px 10px;
    color: {NAVY};
    font-size: 13px;
}}

QLineEdit:focus, QComboBox:focus {{
    border: 1px solid {ACCENT_BLUE};
}}

QPushButton {{
    background-color: {ACCENT_BLUE};
    color: white;
    border: none;
    border-radius: 20px;
    padding: 10px 22px;
    font-weight: 700;
    font-size: 13px;
}}

QPushButton:hover {{
    background-color: {ACCENT_BLUE_HOVER};
}}

QPushButton:disabled {{
    background-color: {BORDER};
    color: {BODY_GRAY};
}}

QPushButton#secondary {{
    background-color: transparent;
    color: {ACCENT_BLUE};
    border: 1px solid {ACCENT_BLUE};
}}

QPushButton#secondary:hover {{
    background-color: #EAF1F9;
}}

QTextEdit {{
    background-color: white;
    border: 1px solid {BORDER};
    border-radius: 8px;
    color: {NAVY};
    font-family: Menlo, Consolas, monospace;
    font-size: 12px;
    padding: 8px;
}}
"""


class PipelineWorker(QThread):
    progress = Signal(str)
    finished_ok = Signal(object)
    failed = Signal(str)

    def __init__(self, species, fastq_folder, parity):
        super().__init__()
        self.species = species
        self.fastq_folder = fastq_folder
        self.parity = parity

    def run(self):
        try:
            self.progress.emit("Fetching reference index...")
            index, t2g = get_reference(self.species)

            self.progress.emit("Finding fastq files...")
            if self.parity == "paired":
                samples, unmatched = find_file_pairs(Path(self.fastq_folder))
            else:
                samples = find_single_files(Path(self.fastq_folder))

            self.progress.emit("Running kb count (this may take a while)...")
            out_dir = run_count(samples, index, t2g, parity=self.parity)

            self.finished_ok.emit(out_dir)
        except Exception as e:
            self.failed.emit(str(e))


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RAW2Compass")
        self.resize(920, 620)

        tabs = QTabWidget()
        tabs.addTab(self._build_run_panel(), "Run Pipeline")
        tabs.addTab(self._build_results_panel(), "Results")
        self.setCentralWidget(tabs)

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

        layout.addWidget(QLabel("FASTQ Folder"))
        folder_row = QHBoxLayout()
        self.folder_input = QLineEdit()
        browse_btn = QPushButton("Browse")
        browse_btn.setObjectName("secondary")
        browse_btn.clicked.connect(self._pick_folder)
        folder_row.addWidget(self.folder_input)
        folder_row.addWidget(browse_btn)
        layout.addLayout(folder_row)

        layout.addWidget(QLabel("Read Type"))
        self.parity_dropdown = QComboBox()
        self.parity_dropdown.addItems(["paired", "single"])
        layout.addWidget(self.parity_dropdown)

        layout.addSpacing(10)

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

    def _build_results_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(30, 30, 30, 30)
        layout.addWidget(QLabel("Normalized matrix preview + export controls go here"))
        panel.setLayout(layout)
        return panel

    def _pick_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select FASTQ folder")
        if folder:
            self.folder_input.setText(folder)

    def _pick_files(self):
        folder = QFileDialog.getExistingDirectory(self, "Select FASTQ files")
        if folder:
            self.folder_input.setText(folder)

    def _start_pipeline(self):
        self.run_btn.setEnabled(False)
        self.status_log.append("Starting pipeline...")

        self.worker = PipelineWorker(
            species=self.species_dropdown.currentText(),
            fastq_folder=self.folder_input.text(),
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


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyleSheet(STYLESHEET)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())