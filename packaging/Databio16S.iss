#define AppName "Databio 16S"
#define AppVersion "1.0.0"
#define AppPublisher "Databio 16S"
#define AppExeName "Databio16S.exe"
#define AppSource AddBackslash(SourcePath) + "..\dist\Databio16S"

[Setup]
AppId={{7E6BFAC8-4F68-44A9-9A2A-58154458F1BD}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher={#AppPublisher}
DefaultDirName={localappdata}\Programs\Databio16S
DefaultGroupName={#AppName}
PrivilegesRequired=lowest
ArchitecturesInstallIn64BitMode=x64
WizardStyle=modern
OutputDir=..\dist
OutputBaseFilename=Databio16S-Setup
Compression=lzma2
SolidCompression=yes
UninstallDisplayIcon={app}\{#AppExeName}

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional shortcuts:"

[Files]
Source: "{#AppSource}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\Databio 16S"; Filename: "{app}\{#AppExeName}"
Name: "{autodesktop}\Databio 16S"; Filename: "{app}\{#AppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExeName}"; Description: "Launch Databio 16S"; Flags: postinstall nowait skipifsilent
