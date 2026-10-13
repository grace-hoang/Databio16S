"""Dependency discovery and setup helpers; no R installation is required."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import app_setup


class DependencyDiscoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def make_executable(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"test")
        return path

    def test_find_rscript_in_standard_install_root(self):
        executable = self.make_executable(self.root / "R-4.6.1" / "bin" / "Rscript.exe")
        found = app_setup.find_rscript(
            search_path="",
            install_roots=[self.root],
            registry_roots=[],
        )
        self.assertTrue(Path(found).samefile(executable))

    def test_find_rscript_in_registry_root(self):
        executable = self.make_executable(self.root / "bin" / "Rscript.exe")
        found = app_setup.find_rscript(
            search_path="",
            install_roots=[],
            registry_roots=[self.root],
        )
        self.assertEqual(Path(found), executable)

    def test_find_cutadapt_prefers_packaged_executable(self):
        bundled = self.make_executable(self.root / "tools" / "cutadapt.exe")
        with patch.object(app_setup.shutil, "which", return_value=""):
            found = app_setup.find_cutadapt(self.root)
        self.assertTrue(Path(found).samefile(bundled))

    def test_find_cutadapt_uses_saved_executable(self):
        preferred = self.make_executable(self.root / "custom" / "cutadapt.exe")
        with patch.object(app_setup.shutil, "which", return_value=""):
            found = app_setup.find_cutadapt(self.root, preferred)
        self.assertTrue(Path(found).samefile(preferred))

    def test_verify_executable_reports_version(self):
        executable = self.make_executable(self.root / "tool.exe")
        with patch.object(
            app_setup.subprocess,
            "run",
            return_value=app_setup.subprocess.CompletedProcess(
                ["tool", "--version"], 0, stdout="Tool 1.2\n", stderr=""
            ),
        ):
            self.assertEqual(app_setup.verify_executable(str(executable)), "Tool 1.2")

    def test_verify_executable_rejects_nonzero_exit(self):
        executable = self.make_executable(self.root / "tool.exe")
        with patch.object(
            app_setup.subprocess,
            "run",
            return_value=app_setup.subprocess.CompletedProcess(
                ["tool", "--version"], 1, stdout="", stderr="not runnable"
            ),
        ):
            with self.assertRaisesRegex(RuntimeError, "not runnable"):
                app_setup.verify_executable(str(executable))


if __name__ == "__main__":
    unittest.main(verbosity=2)
