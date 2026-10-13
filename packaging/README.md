# Windows packaging

## Recommended release format

Use **PyInstaller one-folder + Inno Setup** for the supported Windows release.
The installer bundles the Python runtime, Tkinter application, R scripts,
resources, and a separately frozen Cutadapt executable. Inno Setup then creates
a per-user `Databio16S-Setup.exe`, shortcuts, and an uninstaller. This is more
maintainable than a large PyInstaller one-file executable and avoids requiring
Python or a PATH edit from end users.

Also publish the portable ZIP created by the same script for users who cannot
install software. The ZIP still needs an existing R installation.

## Deliberately external dependency

Keep **R itself external**. It is a large, independently versioned scientific
runtime with platform-specific installers and a broad system footprint.
Bundling it would greatly increase release size and complicate security,
licensing, and update maintenance. RStudio is not required. The app discovers
`Rscript.exe`, guides the user to the official R download if needed, and installs
missing CRAN/Bioconductor packages into the user's R library.

Cutadapt is supplied in the packaged build as its own frozen executable, so it
does not need a Python installation or a PATH entry. Source-code users may still
select an existing Cutadapt executable.

## Build

On a Windows build machine, install Python 3.13 and Inno Setup, then run:

```powershell
.\packaging\build_windows.ps1
```

The script creates an isolated build environment, installs pinned build
dependencies, freezes Cutadapt and the application, assembles a portable ZIP,
and invokes Inno Setup when `ISCC.exe` is available. End users do not need any
of these build tools.

The build machine needs internet access for the pinned Python build
dependencies. End users need internet access on first setup to install missing
R packages. R must be installed separately from the official R for Windows
installer. The wizard uses `Rscript.exe` only; RStudio is neither detected nor
required.

Review `LICENSES.md` and the licenses of all included third-party components
before distributing each build. The repository does not contain a generated
installer or a Cutadapt binary; the build script creates both as appropriate.
