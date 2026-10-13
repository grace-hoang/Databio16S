"""Input validation and subprocess orchestration; no sequencing data needed."""
import importlib.util
import json
from pathlib import Path
import queue
import tempfile
import threading
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("launcher", ROOT / "run_pipeline.py")
launcher = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launcher)


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.folder = Path(self.temp.name) / "Dữ liệu có dấu cách"
        self.folder.mkdir()

    def tearDown(self):
        self.temp.cleanup()

    def make_files(self, names):
        paths = []
        for name in names:
            path = self.folder / name
            path.write_bytes(b"placeholder")
            paths.append(str(path))
        return paths

    def test_pairing_preserves_full_sample_id_and_zero_padding(self):
        files = self.make_files(["001_R1_001.fastq.gz", "001_R2_001.fastq.gz",
                                 "site_A_1.fq", "site_A_2.fq"])
        pairs = launcher.pair_files(files)
        self.assertEqual(set(pairs), {"001", "site_A"})
        self.assertEqual(set(pairs["001"]), {"1", "2"})

    def test_missing_reverse_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Thiếu R1/R2"):
            launcher.pair_files(self.make_files(["sample_R1.fastq.gz"]))

    def test_duplicate_sample_across_folders_is_rejected(self):
        files = self.make_files(["sample_R1.fastq.gz", "sample_R2.fastq.gz"])
        other = self.folder / "other"
        other.mkdir()
        duplicate = other / "sample_R1.fastq.gz"
        duplicate.write_bytes(b"placeholder")
        with self.assertRaisesRegex(ValueError, "Trùng tên mẫu"):
            launcher.pair_files(files + [str(duplicate)])

    def test_non_fastq_name_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Tên FASTQ"):
            launcher.pair_files(self.make_files(["metadata.csv"]))

    def test_empty_file_is_rejected(self):
        path = self.folder / "empty_R1.fastq"
        path.touch()
        with self.assertRaisesRegex(ValueError, "rỗng"):
            launcher.pair_files([str(path)])

    def test_metadata_bom_leading_zeros_and_unicode(self):
        path = self.folder / "metadata.csv"
        path.write_text("SampleID,Group\n001,Đối chứng\n002,Điều trị\n", encoding="utf-8-sig")
        columns, rows = launcher.read_metadata(path)
        self.assertEqual(columns, ["SampleID", "Group"])
        self.assertEqual(rows[0]["SampleID"], "001")
        self.assertEqual(rows[0]["Group"], "Đối chứng")

    def test_duplicate_metadata_sample_is_rejected(self):
        path = self.folder / "metadata.csv"
        path.write_text("SampleID,Group\nsample,A\nsample,B\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duy nhất"):
            launcher.read_metadata(path)

    def test_malformed_csv_is_rejected(self):
        path = self.folder / "metadata.csv"
        path.write_text("SampleID,Group\nsample,A,unexpected\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Số cột"):
            launcher.read_metadata(path)

    def test_translations_cover_english_and_vietnamese(self):
        self.assertEqual(launcher.translate("en", "selected_files", count=3), "3 FASTQ files selected")
        self.assertEqual(launcher.translate("vi", "selected_files", count=3), "Đã chọn 3 file FASTQ")
        with self.assertRaisesRegex(ValueError, "Select paired-end FASTQ"):
            launcher.pair_files([], language="en")

    def test_path_picker_routes_files_and_output_folder_to_dialogs(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.root = object()
        app.language = "en"
        metadata_var = Mock()
        output_var = Mock()
        app.vars = {"metadata_file": metadata_var, "output_dir": output_var}
        app.refresh_metadata = Mock()
        csv_types = [("CSV", "*.csv"), ("All files", "*")]
        with patch.object(launcher.filedialog, "askopenfilename", return_value="metadata.csv") as open_file:
            app.pick_path("metadata_file", csv_types)
        open_file.assert_called_once_with(
            parent=app.root,
            title="Select Metadata CSV",
            filetypes=csv_types,
        )
        metadata_var.set.assert_called_once_with("metadata.csv")
        app.refresh_metadata.assert_called_once_with()

        with patch.object(launcher.filedialog, "askdirectory", return_value="results") as open_folder:
            app.pick_path("output_dir", None)
        open_folder.assert_called_once_with(parent=app.root, title="Select output folder")
        output_var.set.assert_called_once_with("results")

    def test_path_picker_surfaces_tk_dialog_errors(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.root = object()
        app.language = "en"
        app.vars = {"metadata_file": Mock()}
        with patch.object(
            launcher.filedialog,
            "askopenfilename",
            side_effect=launcher.tk.TclError("dialog unavailable"),
        ), patch.object(launcher.messagebox, "showerror") as showerror:
            app.pick_path("metadata_file", [("CSV", "*.csv")])
        showerror.assert_called_once()
        self.assertIn("dialog unavailable", showerror.call_args.args[1])

    def test_metadata_browse_button_invokes_path_picker(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.root = object()
        app.language = "en"
        app.controls = []
        app.vars = {"metadata_file": Mock()}
        app.refresh_metadata = Mock()
        button = Mock()
        app.button = Mock(return_value=button)
        csv_types = [("CSV", "*.csv")]
        with patch.object(launcher.ttk, "Label", return_value=Mock()), patch.object(
            launcher.ttk, "Entry", return_value=Mock()
        ), patch.object(launcher.filedialog, "askopenfilename", return_value="metadata.csv") as open_file:
            app.path_row(object(), 0, "Metadata CSV", "metadata_file", csv_types)
            app.button.call_args.args[2]()
        open_file.assert_called_once()
        app.vars["metadata_file"].set.assert_called_once_with("metadata.csv")
        app.refresh_metadata.assert_called_once_with()

    def test_worker_uses_one_r_process_and_preserves_failure_status(self):
        app = launcher.PipelineApp.__new__(launcher.PipelineApp)
        app.events = queue.Queue()
        app.cancelled = threading.Event()
        cfg = {"rscript": "C:/R path/Rscript.exe", "run_dir": self.folder.as_posix()}
        fake_process = unittest.mock.Mock()
        fake_process.stdout = iter(["[STEP 03]\n", "Cutadapt failed\n"])
        fake_process.wait.return_value = 1
        with patch.object(launcher.subprocess, "Popen", return_value=fake_process) as popen:
            app.worker(cfg, check_only=False)
        popen.assert_called_once()
        command = popen.call_args.args[0]
        self.assertEqual(command[0], cfg["rscript"])
        self.assertEqual(command[-2], str(ROOT / "run_pipeline.R"))
        self.assertEqual(json.loads((self.folder / "status.json").read_text())["exit_code"], 1)
        self.assertIn("Cutadapt failed", (self.folder / "pipeline.log").read_text())


if __name__ == "__main__":
    unittest.main(verbosity=2)
