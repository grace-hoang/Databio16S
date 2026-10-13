"""Runtime dependency discovery and the first-run setup wizard."""

import os
from pathlib import Path
import queue
import shutil
import subprocess
import threading
import tkinter as tk
from tkinter import filedialog, scrolledtext, ttk


WIZARD_TEXT = {
    "en": {
        "title": "Databio 16S — First-run setup",
        "heading": "Preparing Databio 16S",
        "intro": "The app needs R and Cutadapt. R packages can be installed automatically; R itself stays separate.",
        "language": "Language",
        "r_label": "Rscript executable",
        "cutadapt_label": "Cutadapt executable",
        "browse": "Browse…",
        "setup": "Check and install dependencies",
        "checking": "Checking dependencies…",
        "installing": "Installing missing R packages. Internet access is required…",
        "continue": "Continue",
        "download_r": "Download R",
        "close": "Close",
        "r_missing": "Rscript was not found. Install R, then browse to its Rscript.exe.",
        "cutadapt_missing": "Cutadapt was not found. Select cutadapt.exe or rebuild the app with the bundled Cutadapt package.",
        "choose_r": "Select Rscript.exe",
        "choose_cutadapt": "Select cutadapt executable",
        "invalid_r": "The Rscript path is invalid or did not start successfully.",
        "invalid_cutadapt": "The Cutadapt path is invalid or did not start successfully.",
        "package_failed": "R package setup failed. Review the output below, check your internet connection, then try again.",
        "setup_complete": "Dependencies are ready. Continue to the application.",
        "install_notice": "Missing R packages will be downloaded from CRAN/Bioconductor into your R user library.",
        "install_r": "R is not bundled. Install it from the official R for Windows download page.",
        "dependencies": "Dependencies",
        "error": "Setup error",
    },
    "vi": {
        "title": "Databio 16S — Thiết lập lần đầu",
        "heading": "Chuẩn bị Databio 16S",
        "intro": "Ứng dụng cần R và Cutadapt. Có thể tự cài các package R; R được cài riêng.",
        "language": "Ngôn ngữ",
        "r_label": "Tệp thực thi Rscript",
        "cutadapt_label": "Tệp thực thi Cutadapt",
        "browse": "Chọn…",
        "setup": "Kiểm tra và cài đặt thành phần",
        "checking": "Đang kiểm tra thành phần…",
        "installing": "Đang cài package R còn thiếu. Cần kết nối Internet…",
        "continue": "Tiếp tục",
        "download_r": "Tải R",
        "close": "Đóng",
        "r_missing": "Không tìm thấy Rscript. Hãy cài R rồi chọn tệp Rscript.exe.",
        "cutadapt_missing": "Không tìm thấy Cutadapt. Hãy chọn cutadapt.exe hoặc đóng gói lại ứng dụng kèm Cutadapt.",
        "choose_r": "Chọn Rscript.exe",
        "choose_cutadapt": "Chọn tệp thực thi Cutadapt",
        "invalid_r": "Đường dẫn Rscript không hợp lệ hoặc không thể khởi chạy.",
        "invalid_cutadapt": "Đường dẫn Cutadapt không hợp lệ hoặc không thể khởi chạy.",
        "package_failed": "Cài đặt package R thất bại. Xem nội dung bên dưới, kiểm tra Internet rồi thử lại.",
        "setup_complete": "Các thành phần đã sẵn sàng. Tiếp tục vào ứng dụng.",
        "install_notice": "Package R còn thiếu sẽ được tải từ CRAN/Bioconductor vào thư viện R của người dùng.",
        "install_r": "R không được đóng gói cùng ứng dụng. Hãy cài từ trang tải R chính thức cho Windows.",
        "dependencies": "Thành phần cần thiết",
        "error": "Lỗi thiết lập",
    },
}


def _default_r_roots():
    roots = [
        Path(os.environ.get("ProgramFiles", "C:/Program Files")) / "R",
        Path(os.environ.get("ProgramFiles(x86)", "C:/Program Files (x86)")) / "R",
        Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData/Local")) / "Programs/R",
        Path("C:/R"),
    ]
    return roots


def _registry_r_roots():
    if os.name != "nt":
        return []
    try:
        import winreg
    except ImportError:
        return []
    roots = []
    for hive in (winreg.HKEY_CURRENT_USER, winreg.HKEY_LOCAL_MACHINE):
        for registry_path in (r"SOFTWARE\R-core\R64", r"SOFTWARE\R-core\R"):
            try:
                with winreg.OpenKey(hive, registry_path) as key:
                    for index in range(winreg.QueryInfoKey(key)[0]):
                        version = winreg.EnumKey(key, index)
                        with winreg.OpenKey(key, version) as version_key:
                            for value_name in ("InstallPath", ""):
                                try:
                                    value, _ = winreg.QueryValueEx(version_key, value_name)
                                except OSError:
                                    continue
                                if value:
                                    roots.append(Path(value))
                                    break
            except OSError:
                continue
    return roots


def find_rscript(search_path=None, install_roots=None, registry_roots=None):
    """Find Rscript through PATH, common installation locations, or Windows registry."""
    found = shutil.which("Rscript", path=search_path)
    if found:
        return str(Path(found).resolve())

    roots = list(install_roots) if install_roots is not None else _default_r_roots()
    if registry_roots is None:
        roots.extend(_registry_r_roots())
    else:
        roots.extend(Path(root) for root in registry_roots)

    candidates = []
    for root in roots:
        root = Path(root)
        candidates.extend((
            root / "bin" / "Rscript.exe",
            root / "bin" / "Rscript",
            root / "bin" / "x64" / "Rscript.exe",
            root / "bin" / "i386" / "Rscript.exe",
        ))
        try:
            candidates.extend(root.glob("R-*/bin/Rscript.exe"))
            candidates.extend(root.glob("R-*/bin/Rscript"))
            candidates.extend(root.glob("R-*/bin/x64/Rscript.exe"))
            candidates.extend(root.glob("R-*/bin/i386/Rscript.exe"))
        except OSError:
            continue
    existing = []
    for candidate in candidates:
        if candidate.is_file():
            try:
                existing.append((candidate.stat().st_mtime, candidate))
            except OSError:
                continue
    return str(max(existing, key=lambda item: item[0])[1]) if existing else ""


def find_cutadapt(resource_root, preferred_path=""):
    """Find a user-selected or packaged Cutadapt executable."""
    candidates = []
    if preferred_path:
        candidates.append(Path(preferred_path))
    candidates.extend((
        Path(resource_root) / "tools" / "cutadapt.exe",
        Path(resource_root) / "cutadapt.exe",
    ))
    found = shutil.which("cutadapt")
    if found:
        candidates.append(Path(found))
    candidates.append(Path.home() / "Downloads" / "cutadapt.exe")
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate.resolve())
    return ""


def verify_executable(executable, timeout=20):
    """Return version output, or raise a useful error if the command cannot run."""
    if not executable or not Path(executable).is_file():
        raise FileNotFoundError(executable)
    completed = subprocess.run(
        [str(executable), "--version"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
    )
    output = (completed.stdout + "\n" + completed.stderr).strip()
    if completed.returncode != 0 or not output:
        raise RuntimeError(output or f"Process exited with code {completed.returncode}.")
    return output


class SetupWizard:
    def __init__(
        self,
        parent,
        resource_root,
        language,
        rscript_path,
        cutadapt_path,
        setup_completed,
        on_success,
        on_language_change=None,
    ):
        self.parent = parent
        self.resource_root = Path(resource_root)
        self.language = language if language in WIZARD_TEXT else "vi"
        self.text = WIZARD_TEXT[self.language]
        self.on_success = on_success
        self.on_language_change = on_language_change
        self.events = queue.Queue()
        self.busy = False
        self.ready = False
        self.window = tk.Toplevel(parent)
        self.window.title(self.text["title"])
        self.window.geometry("680x510")
        self.window.minsize(580, 440)
        if parent.state() != "withdrawn":
            self.window.transient(parent)
        self.window.protocol("WM_DELETE_WINDOW", self.close)

        self.rscript_var = tk.StringVar(value=rscript_path or find_rscript())
        self.cutadapt_var = tk.StringVar(
            value=cutadapt_path or find_cutadapt(self.resource_root)
        )
        self.status_var = tk.StringVar(value=self.text["intro"])
        self.language_var = tk.StringVar(
            value="English" if self.language == "en" else "Tiếng Việt"
        )
        self._build()
        self.window.after(100, self.poll)
        self.window.grab_set()
        if setup_completed:
            self.begin(install=False, quiet=True)

    def _build(self):
        frame = ttk.Frame(self.window, padding=16)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)
        ttk.Label(frame, text=self.text["heading"], font=("Segoe UI", 14, "bold")).grid(
            row=0, column=0, columnspan=2, sticky="w", pady=(0, 8)
        )
        ttk.Label(frame, text=self.text["language"]).grid(row=0, column=2, sticky="e", padx=4)
        self.language_combo = ttk.Combobox(
            frame,
            textvariable=self.language_var,
            values=("English", "Tiếng Việt"),
            state="readonly",
            width=12,
        )
        self.language_combo.grid(row=0, column=3, sticky="e")
        self.language_combo.bind("<<ComboboxSelected>>", self.change_language)
        ttk.Label(frame, text=self.text["intro"], wraplength=630).grid(
            row=1, column=0, columnspan=4, sticky="w", pady=(0, 12)
        )
        self._path_row(frame, 2, self.text["r_label"], self.rscript_var, self.browse_rscript)
        self._path_row(frame, 3, self.text["cutadapt_label"], self.cutadapt_var, self.browse_cutadapt)
        ttk.Label(frame, text=self.text["install_r"], wraplength=630).grid(
            row=4, column=0, columnspan=4, sticky="w", pady=(10, 4)
        )
        self.download_button = ttk.Button(frame, text=self.text["download_r"], command=self.open_r_download)
        self.download_button.grid(row=5, column=0, sticky="w", pady=(0, 10))
        ttk.Label(frame, text=self.text["install_notice"], wraplength=630).grid(
            row=6, column=0, columnspan=4, sticky="w", pady=(0, 8)
        )
        self.log = scrolledtext.ScrolledText(frame, height=10, wrap="word", state="disabled")
        self.log.grid(row=7, column=0, columnspan=4, sticky="nsew")
        frame.rowconfigure(7, weight=1)
        self.status_label = ttk.Label(frame, textvariable=self.status_var, wraplength=630)
        self.status_label.grid(row=8, column=0, columnspan=4, sticky="w", pady=8)
        self.setup_button = ttk.Button(frame, text=self.text["setup"], command=self.begin)
        self.setup_button.grid(row=9, column=0, sticky="w")
        self.continue_button = ttk.Button(
            frame, text=self.text["continue"], command=self.continue_to_app, state="disabled"
        )
        self.continue_button.grid(row=9, column=1, sticky="e")
        ttk.Button(frame, text=self.text["close"], command=self.close).grid(
            row=9, column=2, sticky="e", padx=(8, 0)
        )

    def _path_row(self, parent, row, label, variable, command):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 8), pady=4)
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, sticky="ew", pady=4)
        ttk.Button(parent, text=self.text["browse"], command=command).grid(
            row=row, column=2, sticky="e", padx=(8, 0), pady=4
        )

    def browse_rscript(self):
        path = filedialog.askopenfilename(
            parent=self.window,
            title=self.text["choose_r"],
            filetypes=[("Rscript", "Rscript.exe"), ("All files", "*")],
        )
        if path:
            self.rscript_var.set(path)

    def browse_cutadapt(self):
        path = filedialog.askopenfilename(
            parent=self.window,
            title=self.text["choose_cutadapt"],
            filetypes=[("Executable", "*.exe"), ("All files", "*")],
        )
        if path:
            self.cutadapt_var.set(path)

    def open_r_download(self):
        import webbrowser

        webbrowser.open("https://cran.r-project.org/bin/windows/base/")

    def change_language(self, _event=None):
        self.language = "en" if self.language_var.get() == "English" else "vi"
        self.text = WIZARD_TEXT[self.language]
        log_text = self.log.get("1.0", "end-1c")
        self.window.title(self.text["title"])
        for child in self.window.winfo_children():
            child.destroy()
        self._build()
        if log_text:
            self.log.configure(state="normal")
            self.log.insert("end", log_text)
            self.log.see("end")
            self.log.configure(state="disabled")
        if self.on_language_change is not None:
            self.on_language_change(self.language)
        if self.ready:
            self.status_var.set(self.text["setup_complete"])
            self.continue_button.configure(state="normal")
        elif not self.busy:
            self.status_var.set(self.text["intro"])

    def begin(self, install=True, quiet=False):
        if self.busy:
            return
        rscript = self.rscript_var.get().strip()
        cutadapt = self.cutadapt_var.get().strip()
        if not rscript or not Path(rscript).is_file():
            self.status_var.set(self.text["r_missing"])
            return
        if not cutadapt or not Path(cutadapt).is_file():
            self.status_var.set(self.text["cutadapt_missing"])
            return
        self.busy = True
        self.ready = False
        self.continue_button.configure(state="disabled")
        self.setup_button.configure(state="disabled")
        self.language_combo.configure(state="disabled")
        self.status_var.set(self.text["installing"] if install else self.text["checking"])
        threading.Thread(
            target=self._run_setup,
            args=(rscript, cutadapt, install, quiet),
            daemon=True,
        ).start()

    def _run_setup(self, rscript, cutadapt, install, quiet):
        try:
            try:
                r_version = verify_executable(rscript)
            except (OSError, RuntimeError, subprocess.TimeoutExpired):
                self.events.put(("failure", "invalid_r"))
                return
            self.events.put(("log", f"Rscript: {r_version}"))
            try:
                cutadapt_version = verify_executable(cutadapt)
            except (OSError, RuntimeError, subprocess.TimeoutExpired):
                self.events.put(("failure", "invalid_cutadapt"))
                return
            self.events.put(("log", f"Cutadapt: {cutadapt_version}"))

            args = [rscript, "--vanilla", str(self.resource_root / "install_packages.R")]
            if not install:
                args.append("--check")
            process = subprocess.Popen(
                args,
                cwd=str(self.resource_root),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
            if process.stdout is not None:
                for line in process.stdout:
                    self.events.put(("log", line.rstrip()))
            code = process.wait()
            if code != 0:
                self.events.put(("failure", "package_failed"))
                return
            self.events.put(("success", (rscript, cutadapt, quiet)))
        except (OSError, subprocess.SubprocessError) as exc:
            self.events.put(("failure_detail", str(exc)))

    def poll(self):
        try:
            while True:
                event, value = self.events.get_nowait()
                if event == "log":
                    self.log.configure(state="normal")
                    self.log.insert("end", value + "\n")
                    self.log.see("end")
                    self.log.configure(state="disabled")
                elif event == "failure":
                    self.busy = False
                    self.setup_button.configure(state="normal")
                    self.language_combo.configure(state="readonly")
                    self.status_var.set(self.text[value])
                elif event == "failure_detail":
                    self.busy = False
                    self.setup_button.configure(state="normal")
                    self.language_combo.configure(state="readonly")
                    self.status_var.set(self.text["package_failed"] + "\n" + value)
                elif event == "success":
                    rscript, cutadapt, quiet = value
                    self.busy = False
                    self.ready = True
                    self.status_var.set(self.text["setup_complete"])
                    self.setup_button.configure(state="normal")
                    self.language_combo.configure(state="readonly")
                    self.continue_button.configure(state="normal")
                    self.on_success(rscript, cutadapt)
                    if quiet:
                        self.continue_to_app()
        except queue.Empty:
            pass
        if self.window.winfo_exists():
            self.window.after(100, self.poll)

    def continue_to_app(self):
        if not self.ready:
            return
        self.window.grab_release()
        self.window.destroy()

    def close(self):
        if self.busy:
            return
        self.window.grab_release()
        self.parent.destroy()
