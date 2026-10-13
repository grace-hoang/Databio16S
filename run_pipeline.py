"""Select paired-end 16S inputs and run the supplied R pipeline.

Python 3.9+ with Tkinter. No pip packages are required.
The analysis runs in one Rscript process, with stdout/stderr saved to disk.
"""

import csv
import json
import os
from pathlib import Path
import queue
import re
import subprocess
import sys
import threading
from datetime import datetime
from functools import partial
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

from app_setup import SetupWizard, find_cutadapt, find_rscript

BASE = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
LOCAL_APP_DATA = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData/Local")
SETTINGS_DIR = LOCAL_APP_DATA / "Databio16S"
SETTINGS = SETTINGS_DIR / "settings.json"
LEGACY_SETTINGS = Path(__file__).resolve().parent / "last_selection.json"
PAIR_PATTERN = re.compile(r"(.+)(?:_R([12])|_([12]))(?:_001)?\.(?:fastq|fq)(?:\.gz)?$")
TEXTS = {
    "en": {
        "file_missing": "File does not exist or is empty: {path}",
        "fastq_name": "FASTQ names must end with _R1/_R2 or _1/_2: {name}",
        "duplicate_sample": "Duplicate sample {sample}, R{direction}. Select one pair per sample.",
        "select_fastq": "Select paired-end FASTQ files to analyze.",
        "missing_mates": "Missing R1/R2 for samples: {samples}",
        "pair_order": "FASTQ pair order differs between folders. Keep R1/R2 in the same folder.",
        "metadata_csv": "Metadata must be a comma-separated CSV with a header and data.",
        "duplicate_columns": "Metadata contains duplicate column names.",
        "metadata_row": "A metadata row has a different number of fields than the CSV header.",
        "sample_ids": "The first metadata column must contain unique, non-empty SampleIDs.",
        "window_title": "Databio 16S — 16S amplicon analysis",
        "heading": "Select paired-end 16S data and configure the analysis",
        "language": "Language",
        "setup_dependencies": "Set up dependencies",
        "setup_complete": "Dependency setup completed.",
        "settings_error": "Could not save settings: {error}",
        "fastq_input": "Input FASTQ",
        "no_files": "No files selected",
        "choose_files": "Choose files…",
        "choose_folder": "Choose folder…",
        "choose": "Browse…",
        "rscript": "Rscript executable",
        "cutadapt": "Cutadapt executable",
        "metadata_file": "Metadata CSV",
        "silva_ref_path": "SILVA training set (DADA2)",
        "primer_fwd": "Forward primer",
        "primer_rev": "Reverse primer",
        "rarefaction_label": "Rarefaction",
        "output_folder": "Output folder",
        "options": "Parameters — check primers and lengths for your experiment",
        "group_column": "Group column",
        "reference_group": "Reference group",
        "asv_min": "ASV min (bp)",
        "asv_max": "ASV max (bp)",
        "threads": "Threads",
        "quality_threshold": "Quality threshold Q",
        "custom_depth": "Depth for custom mode",
        "run_ancombc": "Run ANCOM-BC (step 16)",
        "metadata_hint": "Metadata: first column is SampleID matching FASTQ names; requires ≥2 groups, ≥2 samples/group.",
        "check": "Check environment / data",
        "analyze": "Run analysis",
        "stop": "Stop",
        "open_results": "Open results",
        "ready": "Ready. Check the environment before running.",
        "browse_output": "Select output folder",
        "browse_key": "Select {label}",
        "browse_fastq": "Select all R1 and R2 files (Ctrl/Shift)",
        "all_files": "All files",
        "dialog_error": "Could not open the file picker: {error}",
        "selected_files": "{count} FASTQ files selected",
        "metadata_title": "Metadata",
        "invalid_file": "Select a valid file for {key}.",
        "select_output": "Select an output folder.",
        "invalid_group": "The group column must exist and use a name such as Group/Treatment (letters, numbers, underscore).",
        "missing_metadata": "FASTQ samples missing from metadata: {samples}",
        "invalid_group_value": "Group is empty, NA, or has surrounding whitespace for sample {sample}.",
        "minimum_groups": "Diversity analysis requires at least 2 groups with at least 2 samples per group.",
        "invalid_reference": "The reference group must be present among the selected FASTQ samples.",
        "integer": "{key} must be an integer.",
        "positive": "{key} must be greater than zero.",
        "asv_order": "ASV min must not exceed ASV max.",
        "rarefaction": "Choose rarefaction auto_rarefy, min, or custom.",
        "primer": "{key} must be a DNA/IUPAC sequence with no spaces.",
        "scripts": "The pipeline folder beside this Python file must contain all 17 scripts from 00 to 16.",
        "cannot_run": "Cannot run",
        "checking": "Checking…",
        "analyzing": "Analyzing…",
        "run_folder": "Run folder: {path}",
        "launch_error": "Launch error: {error}",
        "stopped_status": "Stopped. Current results are incomplete.",
        "check_success": "Check completed successfully.",
        "analysis_success": "Analysis completed.",
        "complete": "Complete",
        "failed": "Failed (code {code}). See the log below.",
        "pipeline_stopped": "Pipeline stopped",
        "log_at": "See pipeline.log at:\n{path}",
        "stopping": "Stopping…",
        "results": "Results",
        "running": "Running",
        "stop_before_close": "Stop the pipeline before closing the window.",
    },
    "vi": {
        "file_missing": "File không tồn tại hoặc rỗng: {path}",
        "fastq_name": "Tên FASTQ cần kết thúc bằng _R1/_R2 hoặc _1/_2: {name}",
        "duplicate_sample": "Trùng tên mẫu {sample}, R{direction}. Chọn một cặp cho mỗi mẫu.",
        "select_fastq": "Hãy chọn FASTQ paired-end cần phân tích.",
        "missing_mates": "Thiếu R1/R2 cho mẫu: {samples}",
        "pair_order": "Thứ tự cặp FASTQ khác nhau giữa các thư mục. Đặt R1/R2 cùng thư mục.",
        "metadata_csv": "Metadata phải là CSV phân cách bằng dấu phẩy, có tiêu đề và dữ liệu.",
        "duplicate_columns": "Metadata có tên cột trùng nhau.",
        "metadata_row": "Số cột trong một dòng metadata không khớp với tiêu đề CSV.",
        "sample_ids": "Cột đầu metadata phải chứa SampleID duy nhất và không được trống.",
        "window_title": "Databio 16S — Phân tích dữ liệu 16S",
        "heading": "Chọn dữ liệu 16S paired-end và thiết lập phân tích",
        "language": "Ngôn ngữ",
        "setup_dependencies": "Thiết lập thành phần",
        "setup_complete": "Đã thiết lập xong các thành phần.",
        "settings_error": "Không thể lưu cài đặt: {error}",
        "fastq_input": "FASTQ đầu vào",
        "no_files": "Chưa chọn file",
        "choose_files": "Chọn file…",
        "choose_folder": "Chọn thư mục…",
        "choose": "Chọn…",
        "rscript": "Rscript.exe",
        "cutadapt": "Cutadapt executable",
        "metadata_file": "Metadata CSV",
        "silva_ref_path": "SILVA trainset (DADA2)",
        "primer_fwd": "Primer forward",
        "primer_rev": "Primer reverse",
        "rarefaction_label": "Rarefaction",
        "output_folder": "Thư mục lưu kết quả",
        "options": "Thông số — kiểm tra primer và chiều dài theo thí nghiệm",
        "group_column": "Cột nhóm",
        "reference_group": "Nhóm tham chiếu",
        "asv_min": "ASV min (bp)",
        "asv_max": "ASV max (bp)",
        "threads": "Số luồng",
        "quality_threshold": "Ngưỡng chất lượng Q",
        "custom_depth": "Depth khi chọn custom",
        "run_ancombc": "Chạy ANCOM-BC (bước 16)",
        "metadata_hint": "Metadata: cột đầu là SampleID khớp tên FASTQ; cần ≥2 nhóm, ≥2 mẫu/nhóm.",
        "check": "Kiểm tra môi trường / dữ liệu",
        "analyze": "Chạy phân tích",
        "stop": "Dừng",
        "open_results": "Mở kết quả",
        "ready": "Sẵn sàng. Nên kiểm tra môi trường trước khi chạy.",
        "browse_output": "Chọn thư mục lưu kết quả",
        "browse_key": "Chọn {label}",
        "browse_fastq": "Chọn tất cả file R1 và R2 (Ctrl/Shift)",
        "all_files": "Tất cả",
        "dialog_error": "Không thể mở hộp thoại chọn file: {error}",
        "selected_files": "Đã chọn {count} file FASTQ",
        "metadata_title": "Metadata",
        "invalid_file": "Hãy chọn file hợp lệ cho {key}.",
        "select_output": "Hãy chọn thư mục lưu kết quả.",
        "invalid_group": "Cột nhóm phải tồn tại và có tên như Group/Treatment (chữ, số, dấu _).",
        "missing_metadata": "Mẫu FASTQ không có trong metadata: {samples}",
        "invalid_group_value": "Nhóm trống, NA hoặc thừa khoảng trắng cho mẫu {sample}.",
        "minimum_groups": "Pipeline so sánh đa dạng cần ít nhất 2 nhóm, mỗi nhóm ít nhất 2 mẫu.",
        "invalid_reference": "Nhóm tham chiếu phải có trong các mẫu FASTQ đã chọn.",
        "integer": "{key} phải là số nguyên.",
        "positive": "{key} phải lớn hơn 0.",
        "asv_order": "ASV min không được lớn hơn ASV max.",
        "rarefaction": "Chọn rarefaction auto_rarefy, min hoặc custom.",
        "primer": "{key} phải là trình tự DNA/IUPAC, không có khoảng trắng.",
        "scripts": "Thư mục pipeline cạnh file Python cần đủ 17 script từ 00 tới 16.",
        "cannot_run": "Chưa thể chạy",
        "checking": "Đang kiểm tra…",
        "analyzing": "Đang phân tích…",
        "run_folder": "Thư mục: {path}",
        "launch_error": "Lỗi khởi chạy: {error}",
        "stopped_status": "Đã dừng. Kết quả hiện tại chưa hoàn chỉnh.",
        "check_success": "Kiểm tra thành công.",
        "analysis_success": "Phân tích hoàn tất.",
        "complete": "Hoàn tất",
        "failed": "Thất bại (mã {code}). Xem log bên dưới.",
        "pipeline_stopped": "Pipeline dừng",
        "log_at": "Xem pipeline.log tại:\n{path}",
        "stopping": "Đang dừng…",
        "results": "Kết quả",
        "running": "Đang chạy",
        "stop_before_close": "Nhấn Dừng trước khi đóng cửa sổ.",
    },
}


def translate(language, key, **values):
    return TEXTS[language][key].format(**values)


def save_settings(settings):
    SETTINGS_DIR.mkdir(parents=True, exist_ok=True)
    SETTINGS.write_text(json.dumps(settings, ensure_ascii=False, indent=2), encoding="utf-8")


def pair_files(files, language="vi"):
    """Reject missing mates, duplicate SampleIDs, and unsupported names."""
    pairs = {}
    for filename in sorted(set(files)):
        path = Path(filename).resolve()
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(translate(language, "file_missing", path=path))
        match = PAIR_PATTERN.fullmatch(path.name)
        if not match:
            raise ValueError(translate(language, "fastq_name", name=path.name))
        sample, direction = match.group(1), match.group(2) or match.group(3)
        if sample in pairs and direction in pairs[sample]:
            raise ValueError(translate(language, "duplicate_sample", sample=sample, direction=direction))
        pairs.setdefault(sample, {})[direction] = path.as_posix()
    if not pairs:
        raise ValueError(translate(language, "select_fastq"))
    incomplete = [sample for sample, mates in pairs.items() if set(mates) != {"1", "2"}]
    if incomplete:
        raise ValueError(translate(language, "missing_mates", samples=", ".join(incomplete[:15])))
    # Step 02 sorts full paths before comparing names. Keep that contract.
    forward = sorted(mates["1"] for mates in pairs.values())
    reverse = sorted(mates["2"] for mates in pairs.values())
    forward_ids = [PAIR_PATTERN.fullmatch(Path(p).name).group(1) for p in forward]
    reverse_ids = [PAIR_PATTERN.fullmatch(Path(p).name).group(1) for p in reverse]
    if forward_ids != reverse_ids:
        raise ValueError(translate(language, "pair_order"))
    return pairs


def read_metadata(path, language="vi"):
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames or []
        rows = list(reader)
    if len(columns) < 2 or not rows:
        raise ValueError(translate(language, "metadata_csv"))
    if len(set(columns)) != len(columns):
        raise ValueError(translate(language, "duplicate_columns"))
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(translate(language, "metadata_row"))
    sample_ids = [row[columns[0]] for row in rows]
    if any(not value.strip() for value in sample_ids) or len(set(sample_ids)) != len(sample_ids):
        raise ValueError(translate(language, "sample_ids"))
    return columns, rows


def _parse_bool(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


class PipelineApp:
    def __init__(self, root):
        self.root = root
        self.root.geometry("1050x850")
        self.root.minsize(850, 720)
        self.events = queue.Queue()
        self.process = None
        self.busy = False
        self.cancelled = threading.Event()
        self.fastq_files = []
        self.last_run = None
        self.vars = {}
        self.controls = []
        self.language = "vi"
        self.language_var = tk.StringVar(value="English" if self.language == "en" else "Tiếng Việt")
        self.status_key = "ready"
        self.status_args = {}
        self.status = tk.StringVar()
        self.saved = {}
        settings_path = SETTINGS if SETTINGS.exists() else LEGACY_SETTINGS
        try:
            self.saved = json.loads(settings_path.read_text(encoding="utf-8"))
            if not isinstance(self.saved, dict):
                self.saved = {}
        except (OSError, ValueError, TypeError):
            self.saved = {}
        saved_language = self.saved.get("language", "vi")
        self.language = saved_language if saved_language in TEXTS else "vi"
        self.language_var.set("English" if self.language == "en" else "Tiếng Việt")
        defaults = {
            "rscript": find_rscript(), "cutadapt": find_cutadapt(BASE),
            "metadata_file": "", "silva_ref_path": "",
            "output_dir": str(Path.home() / "Documents" / "Databio16S" / "runs"), "group_column": "Group",
            "reference_group": "", "primer_fwd": "GGATTAGATACCCTGGTA",
            "primer_rev": "CRRCACGAGCTGACGAC", "asv_min_len": "230",
            "asv_max_len": "280", "trunc_qmin": "30", "threads": "1",
            "rarefaction_mode": "auto_rarefy", "custom_rarefaction_depth": "10000",
        }
        defaults.update({k: str(v) for k, v in self.saved.items() if k in defaults})
        if not Path(defaults["rscript"]).is_file():
            defaults["rscript"] = find_rscript()
        if not Path(defaults["cutadapt"]).is_file():
            defaults["cutadapt"] = find_cutadapt(BASE)
        self.fastq_files = self.saved.get("fastq_files", [])
        self.vars = {k: tk.StringVar(value=v) for k, v in defaults.items()}
        self.run_ancombc = tk.BooleanVar(value=_parse_bool(self.saved.get("run_ancombc", True)))
        self.build_ui()
        self.refresh_fastq_label()
        if self.vars["metadata_file"].get():
            self.refresh_metadata(show_errors=False)
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.root.after(100, self.poll)

    def t(self, key, **values):
        return translate(self.language, key, **values)

    def set_status(self, key, **values):
        self.status_key = key
        self.status_args = values
        self.status.set(self.t(key, **values))

    def build_ui(self):
        self.controls = []
        self.root.title(self.t("window_title"))
        self.main_frame = ttk.Frame(self.root, padding=12)
        self.main_frame.pack(fill="both", expand=True)
        outer = self.main_frame
        outer.columnconfigure(1, weight=1)
        outer.columnconfigure(4, weight=0)
        ttk.Label(outer, text=self.t("heading"),
                  font=("Segoe UI", 13, "bold")).grid(row=0, column=0, columnspan=3,
                                                       sticky="w", pady=(0, 10))
        ttk.Label(outer, text=self.t("language")).grid(row=0, column=3, sticky="e", padx=(8, 4))
        language_combo = ttk.Combobox(outer, textvariable=self.language_var,
                                      values=("English", "Tiếng Việt"), state="readonly", width=12)
        language_combo.grid(row=0, column=4, sticky="e", pady=(0, 10))
        language_combo.bind("<<ComboboxSelected>>", self.change_language)
        self.controls.append(language_combo)
        ttk.Label(outer, text=self.t("fastq_input")).grid(row=1, column=0, sticky="w", padx=(0, 10))
        self.fastq_label = ttk.Label(outer, text=self.t("no_files"))
        self.fastq_label.grid(row=1, column=1, sticky="w")
        buttons = ttk.Frame(outer)
        buttons.grid(row=1, column=2, pady=4)
        self.button(buttons, self.t("choose_files"), self.select_fastq).pack(side="left", padx=2)
        self.button(buttons, self.t("choose_folder"), self.select_fastq_folder).pack(side="left")
        self.path_row(outer, 2, self.t("metadata_file"), "metadata_file", [("CSV", "*.csv"), (self.t("all_files"), "*")])
        self.path_row(outer, 3, self.t("silva_ref_path"), "silva_ref_path",
                      [("FASTA / gzip", "*.fa *.fasta *.fa.gz *.fasta.gz"), (self.t("all_files"), "*")])
        self.path_row(outer, 4, self.t("output_folder"), "output_dir", None)
        options = ttk.LabelFrame(outer, text=self.t("options"), padding=8)
        options.grid(row=5, column=0, columnspan=3, sticky="ew", pady=8)
        options.columnconfigure(1, weight=1)
        options.columnconfigure(3, weight=1)
        self.group_combo = self.option(options, 0, 0, self.t("group_column"), "group_column", [])
        self.ref_combo = self.option(options, 0, 2, self.t("reference_group"), "reference_group", [])
        self.group_combo.bind("<<ComboboxSelected>>", lambda event: self.refresh_groups())
        self.option(options, 1, 0, self.t("primer_fwd"), "primer_fwd")
        self.option(options, 1, 2, self.t("primer_rev"), "primer_rev")
        self.option(options, 2, 0, self.t("asv_min"), "asv_min_len")
        self.option(options, 2, 2, self.t("asv_max"), "asv_max_len")
        self.option(options, 3, 0, self.t("threads"), "threads")
        self.option(options, 3, 2, self.t("quality_threshold"), "trunc_qmin")
        self.option(options, 4, 0, self.t("rarefaction_label"), "rarefaction_mode", ["auto_rarefy", "min", "custom"])
        self.option(options, 4, 2, self.t("custom_depth"), "custom_rarefaction_depth")
        checkbox = ttk.Checkbutton(options, text=self.t("run_ancombc"), variable=self.run_ancombc)
        checkbox.grid(row=5, column=0, columnspan=4, sticky="w", pady=5)
        self.controls.append(checkbox)
        ttk.Label(outer, text=self.t("metadata_hint"))\
            .grid(row=6, column=0, columnspan=3, sticky="w", pady=3)
        actions = ttk.Frame(outer)
        actions.grid(row=7, column=0, columnspan=3, sticky="w", pady=8)
        self.button(actions, self.t("check"), lambda: self.start(check_only=True)).pack(side="left", padx=(0, 8))
        self.button(actions, self.t("analyze"), self.start).pack(side="left", padx=(0, 8))
        self.stop_button = ttk.Button(actions, text=self.t("stop"), command=self.stop, state="disabled")
        self.stop_button.pack(side="left", padx=(0, 8))
        ttk.Button(actions, text=self.t("open_results"), command=self.open_results).pack(side="left")
        ttk.Button(actions, text=self.t("setup_dependencies"), command=self.open_setup).pack(side="left", padx=(8, 0))
        self.status.set(self.t(self.status_key, **self.status_args))
        ttk.Label(outer, textvariable=self.status).grid(row=8, column=0, columnspan=3, sticky="w")
        self.log = scrolledtext.ScrolledText(outer, height=12, wrap="word", state="disabled")
        self.log.grid(row=9, column=0, columnspan=3, sticky="nsew", pady=(8, 0))
        outer.rowconfigure(9, weight=1)

    def button(self, parent, text, command):
        widget = ttk.Button(parent, text=text, command=command)
        self.controls.append(widget)
        return widget

    def change_language(self, _event=None):
        language = "en" if self.language_var.get() == "English" else "vi"
        self.set_language(language)

    def set_language(self, language):
        if language not in TEXTS or language == self.language:
            return
        self.language = language
        self.language_var.set("English" if language == "en" else "Tiếng Việt")
        log_text = self.log.get("1.0", "end-1c") if hasattr(self, "log") else ""
        try:
            self.saved["language"] = self.language
            save_settings(self.saved)
        except OSError as exc:
            messagebox.showerror(
                self.t("cannot_run"), self.t("settings_error", error=exc), parent=self.root
            )
        if hasattr(self, "main_frame"):
            self.main_frame.destroy()
        self.build_ui()
        if log_text:
            self.log.configure(state="normal")
            self.log.insert("end", log_text)
            self.log.see("end")
            self.log.configure(state="disabled")
        self.refresh_fastq_label()
        if self.vars["metadata_file"].get():
            self.refresh_metadata(show_errors=False)

    def open_setup(self, setup_completed=False):
        if self.busy:
            messagebox.showinfo(self.t("running"), self.t("stop_before_close"), parent=self.root)
            return None
        self.setup_wizard = SetupWizard(
            self.root,
            BASE,
            self.language,
            self.vars["rscript"].get(),
            self.vars["cutadapt"].get(),
            setup_completed,
            self.setup_completed,
            self.set_language,
        )
        return self.setup_wizard

    def setup_completed(self, rscript, cutadapt):
        self.vars["rscript"].set(rscript)
        self.vars["cutadapt"].set(cutadapt)
        self.saved.update({
            "rscript": rscript,
            "cutadapt": cutadapt,
            "setup_completed": True,
            "language": self.language,
        })
        try:
            save_settings(self.saved)
        except OSError as exc:
            messagebox.showerror(
                self.t("cannot_run"), self.t("settings_error", error=exc), parent=self.root
            )
        self.root.deiconify()

    def path_row(self, parent, row, label, key, types):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 10))
        entry = ttk.Entry(parent, textvariable=self.vars[key])
        entry.grid(row=row, column=1, sticky="ew", pady=4)
        self.controls.append(entry)
        self.button(parent, self.t("choose"), partial(self.pick_path, key, types)).grid(
            row=row, column=2, sticky="ew", padx=(8, 0)
        )

    def option(self, parent, row, col, label, key, choices=None):
        ttk.Label(parent, text=label).grid(row=row, column=col, sticky="w", padx=5, pady=3)
        widget = (ttk.Entry(parent, textvariable=self.vars[key]) if choices is None else
                  ttk.Combobox(parent, textvariable=self.vars[key], values=choices))
        widget.grid(row=row, column=col + 1, sticky="ew", padx=5, pady=3)
        self.controls.append(widget)
        return widget

    def pick_path(self, key, types):
        try:
            if types is None:
                path = filedialog.askdirectory(
                    parent=self.root, title=self.t("browse_output")
                )
            else:
                path = filedialog.askopenfilename(
                    parent=self.root,
                    title=self.t("browse_key", label=self.t(key)),
                    filetypes=types,
                )
        except (OSError, tk.TclError) as exc:
            messagebox.showerror(
                self.t("metadata_title"),
                self.t("dialog_error", error=exc),
                parent=self.root,
            )
            return
        if path:
            self.vars[key].set(path)
            if key == "metadata_file":
                self.refresh_metadata()

    def select_fastq(self):
        files = filedialog.askopenfilenames(parent=self.root, title=self.t("browse_fastq"),
                                            filetypes=[("FASTQ", "*.fastq *.fq *.fastq.gz *.fq.gz"), (self.t("all_files"), "*")])
        if files:
            self.fastq_files = list(files)
            self.refresh_fastq_label()

    def select_fastq_folder(self):
        folder = filedialog.askdirectory(parent=self.root, title=self.t("browse_fastq"))
        if not folder:
            return
        files = []
        excluded = {"result", "results", "runs", "filtered", "primer_trimmed", ".git"}
        output_text = self.vars["output_dir"].get().strip()
        output = Path(output_text).resolve() if output_text else None
        for directory, dirs, names in os.walk(folder):
            dirs[:] = [name for name in dirs if name not in excluded and
                       (output is None or (Path(directory) / name).resolve() != output)]
            files.extend(str(Path(directory) / name) for name in names
                         if name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")))
        self.fastq_files = sorted(files)
        self.refresh_fastq_label()
        metadata = Path(folder) / "metadata.csv"
        if metadata.is_file() and not self.vars["metadata_file"].get():
            self.vars["metadata_file"].set(str(metadata))
            self.refresh_metadata()

    def refresh_fastq_label(self):
        self.fastq_label.configure(text=self.t("selected_files", count=len(self.fastq_files)))

    def refresh_metadata(self, show_errors=True):
        try:
            columns, _ = read_metadata(self.vars["metadata_file"].get(), self.language)
            self.group_combo.configure(values=columns[1:])
            if self.vars["group_column"].get() not in columns[1:]:
                self.vars["group_column"].set(columns[1])
            self.refresh_groups()
        except (OSError, ValueError, UnicodeError) as exc:
            if show_errors:
                messagebox.showerror(self.t("metadata_title"), str(exc), parent=self.root)

    def refresh_groups(self):
        try:
            _, rows = read_metadata(self.vars["metadata_file"].get(), self.language)
            groups = sorted({row.get(self.vars["group_column"].get(), "") for row in rows} - {""})
            self.ref_combo.configure(values=groups)
            if self.vars["reference_group"].get() not in groups:
                self.vars["reference_group"].set(groups[0] if groups else "")
        except (OSError, ValueError, UnicodeError):
            pass

    def configuration(self):
        cfg = {key: var.get().strip() for key, var in self.vars.items()}
        for key in ("rscript", "cutadapt", "metadata_file", "silva_ref_path"):
            if not cfg[key] or not Path(cfg[key]).is_file():
                raise ValueError(self.t("invalid_file", key=key))
            cfg[key] = Path(cfg[key]).resolve().as_posix()
        if not cfg["output_dir"]:
            raise ValueError(self.t("select_output"))
        cfg["output_dir"] = Path(cfg["output_dir"]).resolve().as_posix()
        pairs = pair_files(self.fastq_files, self.language)
        columns, rows = read_metadata(cfg["metadata_file"], self.language)
        group = cfg["group_column"]
        if group not in columns[1:] or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", group):
            raise ValueError(self.t("invalid_group"))
        meta = {row[columns[0]]: row for row in rows}
        missing = set(pairs) - set(meta)
        if missing:
            raise ValueError(self.t("missing_metadata", samples=", ".join(sorted(missing)[:15])))
        counts = {}
        for sample in pairs:
            value = meta[sample][group]
            if not value or value != value.strip() or value == "NA":
                raise ValueError(self.t("invalid_group_value", sample=sample))
            counts[value] = counts.get(value, 0) + 1
        if len(counts) < 2 or min(counts.values()) < 2:
            raise ValueError(self.t("minimum_groups"))
        if cfg["reference_group"] not in counts:
            raise ValueError(self.t("invalid_reference"))
        for key in ("asv_min_len", "asv_max_len", "threads", "trunc_qmin", "custom_rarefaction_depth"):
            try:
                cfg[key] = int(cfg[key])
            except ValueError:
                raise ValueError(self.t("integer", key=key)) from None
            if cfg[key] <= 0:
                raise ValueError(self.t("positive", key=key))
        if cfg["asv_min_len"] > cfg["asv_max_len"]:
            raise ValueError(self.t("asv_order"))
        if cfg["rarefaction_mode"] not in {"auto_rarefy", "min", "custom"}:
            raise ValueError(self.t("rarefaction"))
        for key in ("primer_fwd", "primer_rev"):
            cfg[key] = cfg[key].upper()
            if not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", cfg[key]):
                raise ValueError(self.t("primer", key=key))
        cfg["fastq_files"] = [path for sample in sorted(pairs) for path in pairs[sample].values()]
        cfg["project_dir"] = Path(cfg["fastq_files"][0]).parent.as_posix()
        cfg["pipeline_dir"] = (BASE / "pipeline").as_posix()
        cfg["run_ancombc"] = self.run_ancombc.get()
        cfg["language"] = self.language
        scripts = sorted((BASE / "pipeline").glob("[0-9][0-9]_*.R"))
        if len(scripts) != 17 or [p.name[:2] for p in scripts] != [f"{i:02}" for i in range(17)]:
            raise ValueError(self.t("scripts"))
        return cfg

    def start(self, check_only=False):
        if self.busy:
            return
        try:
            cfg = self.configuration()
            prefix = "check" if check_only else "analysis"
            run_dir = Path(cfg["output_dir"]) / f"{prefix}_{datetime.now():%Y%m%d_%H%M%S_%f}"
            run_dir.mkdir(parents=True, exist_ok=False)
            cfg["run_dir"] = run_dir.as_posix()
            (run_dir / "config.json").write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8")
            with (run_dir / "input_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["SampleID", "R1", "R2"])
                for sample, mates in sorted(pair_files(cfg["fastq_files"], self.language).items()):
                    writer.writerow([sample, mates["1"], mates["2"]])
            try:
                save_settings({**self.saved, **cfg, "setup_completed": True})
            except OSError as exc:
                messagebox.showerror(
                    self.t("cannot_run"), self.t("settings_error", error=exc), parent=self.root
                )
        except (OSError, ValueError, UnicodeError) as exc:
            messagebox.showerror(self.t("cannot_run"), str(exc), parent=self.root)
            return
        self.last_run = run_dir
        self.busy = True
        self.cancelled.clear()
        for widget in self.controls:
            widget.configure(state="disabled")
        self.stop_button.configure(state="normal")
        self.set_status("checking" if check_only else "analyzing")
        self.append_log("\n" + self.t("run_folder", path=run_dir) + "\n")
        threading.Thread(target=self.worker, args=(cfg, check_only), daemon=True).start()

    def worker(self, cfg, check_only):
        run_dir = Path(cfg["run_dir"])
        command = [cfg["rscript"], "--no-save", "--no-restore", "--encoding=UTF-8",
                   str(BASE / "run_pipeline.R"), str(run_dir / "config.json")]
        if check_only:
            command.append("--check")
        try:
            with (run_dir / "pipeline.log").open("w", encoding="utf-8") as log:
                process = subprocess.Popen(command, cwd=str(BASE), stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT, text=True, encoding="utf-8",
                                           errors="replace", creationflags=(subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0))
                self.process = process
                if self.cancelled.is_set():
                    self.kill_process(process)
                for line in process.stdout:
                    log.write(line)
                    log.flush()
                    self.events.put(("log", line))
                code = process.wait()
            result = {"exit_code": code, "cancelled": self.cancelled.is_set(), "check_only": check_only}
            (run_dir / "status.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
            self.events.put(("done", result))
        except Exception as exc:
            error = "\n" + self.t("launch_error", error=exc) + "\n"
            result = {"exit_code": -1, "cancelled": self.cancelled.is_set(), "check_only": check_only}
            try:
                with (run_dir / "pipeline.log").open("a", encoding="utf-8") as log:
                    log.write(error)
                (run_dir / "status.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
            except OSError:
                pass
            self.events.put(("log", error))
            self.events.put(("done", result))
        finally:
            self.process = None

    def poll(self):
        try:
            while True:
                kind, data = self.events.get_nowait()
                if kind == "log":
                    self.append_log(data)
                    if "[STEP " in data:
                        self.status.set(data.strip())
                elif kind == "done":
                    self.busy = False
                    for widget in self.controls:
                        widget.configure(state="normal")
                    self.stop_button.configure(state="disabled")
                    if data["cancelled"]:
                        self.set_status("stopped_status")
                    elif data["exit_code"] == 0:
                        self.set_status("check_success" if data["check_only"] else "analysis_success")
                        messagebox.showinfo(self.t("complete"), self.status.get() + f"\n{self.last_run}", parent=self.root)
                    else:
                        self.set_status("failed", code=data["exit_code"])
                        messagebox.showerror(self.t("pipeline_stopped"), self.t("log_at", path=self.last_run), parent=self.root)
        except queue.Empty:
            pass
        self.root.after(100, self.poll)

    def append_log(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text)
        # Keep the interface responsive during verbose, long analyses.
        if int(self.log.index("end-1c").split(".")[0]) > 4000:
            self.log.delete("1.0", "1000.0")
        self.log.see("end")
        self.log.configure(state="disabled")

    @staticmethod
    def kill_process(process):
        if process.poll() is not None:
            return
        if os.name == "nt":
            subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                           capture_output=True, creationflags=subprocess.CREATE_NO_WINDOW)
        else:
            process.terminate()

    def stop(self):
        self.cancelled.set()
        self.stop_button.configure(state="disabled")
        self.set_status("stopping")
        if self.process is not None:
            threading.Thread(target=self.kill_process, args=(self.process,), daemon=True).start()

    def open_results(self):
        path = self.last_run or Path(self.vars["output_dir"].get())
        if Path(path).is_dir() and os.name == "nt":
            os.startfile(str(path))
        else:
            messagebox.showinfo(self.t("results"), str(path), parent=self.root)

    def close(self):
        if self.busy:
            messagebox.showinfo(self.t("running"), self.t("stop_before_close"), parent=self.root)
            return
        self.root.destroy()


if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()
    app = PipelineApp(root)
    app.open_setup(setup_completed=_parse_bool(app.saved.get("setup_completed")))
    root.mainloop()
