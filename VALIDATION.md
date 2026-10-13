# Validation status

## Current source checks

- Python syntax checks passed for `run_pipeline.py`, `app_setup.py`, and
  `tests/test_app_setup.py`.
- All 19 Python unit tests pass, including Rscript discovery, Cutadapt lookup,
  executable failure handling, FASTQ pairing, metadata validation, and
  file/folder browse-dialog routing and error reporting.
- A Tkinter smoke test confirmed that the first-run wizard is visible while the
  main window is hidden, that it can change the UI language, and that it closes
  cleanly.
- The Windows packaging PowerShell script parses successfully.
- PyInstaller built the `dist/Databio16S` one-folder application using Python
  3.13.15, including the 16S R scripts, resources, and frozen Cutadapt.
- The frozen Cutadapt executable returned version `5.2`.
- The packaged GUI executable stayed running and displayed the first-run setup
  wizard.
- Inno Setup 7 produced `dist/Databio16S-Setup.exe`; the build also produced
  `dist/Databio16S-portable.zip`.
- The browse-dialog fix was included in a successful PyInstaller and Inno Setup
  rebuild. The installer and portable ZIP were regenerated.

## Not verified on this machine

- No R runtime is available here, so the first-run CRAN/Bioconductor installation
  and full 16S analysis were not exercised.
- The generated installer was compiled but not run through an install/uninstall
  cycle. The setup wizard still requires users to install R separately and
  internet access to install missing R packages.
- Review third-party license and notice files in each produced package before
  distribution.

The package/discovery unit tests use temporary files or mocked subprocess
responses. The Cutadapt version and packaged GUI checks exercised the actual
build output, but they do not establish that R/Bioconductor works on a user's
computer or validate scientific results from research data.
