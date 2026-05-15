#Requires -Version 5.1
<#
.SYNOPSIS
    P2MAT Windows Installer  –  v1.0.0
    Developer: Methun Kamruzzaman
    Requires : Windows 10 (64-bit, build 1903+) or Windows 11

.DESCRIPTION
    Standalone installer that:
      1. Verifies OS and architecture
      2. Installs Java 21 LTS (Microsoft OpenJDK) if needed
      3. Installs Miniconda3 if needed
      4. Creates the 'qsai' conda environment
      5. Installs all Python dependencies (including lightgbm, torch)
      6. Copies P2MAT to  %LOCALAPPDATA%\P2MAT
      7. Creates a desktop shortcut and Start Menu entry

.EXAMPLE
    Right-click Install-P2MAT.ps1 → "Run with PowerShell"
    (Or: powershell -ExecutionPolicy Bypass -File Install-P2MAT.ps1)
#>

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# ── Helpers ──────────────────────────────────────────────────
function Write-Banner([string]$msg) {
    Write-Host ""
    Write-Host ("=" * 50) -ForegroundColor Cyan
    Write-Host "  $msg" -ForegroundColor Cyan
    Write-Host ("=" * 50) -ForegroundColor Cyan
}
function Write-Step([string]$msg)  { Write-Host "`n$msg" -ForegroundColor White }
function Write-OK([string]$msg)    { Write-Host "  [OK] $msg" -ForegroundColor Green }
function Write-Warn([string]$msg)  { Write-Host "  [!]  $msg" -ForegroundColor Yellow }
function Write-Fail([string]$msg)  { Write-Host "`n[ERROR] $msg`n" -ForegroundColor Red; Read-Host "Press ENTER to exit"; exit 1 }

function Get-FileFromWeb([string]$Url, [string]$Dest) {
    Write-Host "  Downloading $(Split-Path $Url -Leaf) ..."
    $wc = New-Object System.Net.WebClient
    $wc.DownloadFile($Url, $Dest)
}

function Test-CommandExists([string]$cmd) {
    return [bool](Get-Command $cmd -ErrorAction SilentlyContinue)
}

# ── Configuration ────────────────────────────────────────────
$VERSION     = "1.0.0"
$ENV_NAME    = "qsai"
$SCRIPT_DIR  = Split-Path -Parent $MyInvocation.MyCommand.Definition
$SOURCE_DIR  = Join-Path $SCRIPT_DIR "P2MAT"
$INSTALL_DIR = Join-Path $env:LOCALAPPDATA "P2MAT"
$LOG_DIR     = Join-Path $env:LOCALAPPDATA "P2MAT\Logs"
$MINI_DIR    = Join-Path $env:USERPROFILE "miniconda3"
$CONDA_EXE   = Join-Path $MINI_DIR "Scripts\conda.exe"
$CONDA_HOOK  = Join-Path $MINI_DIR "shell\condabin\conda-hook.ps1"

# URLs
$MINICONDA_URL  = "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe"
$JAVA_URL       = "https://aka.ms/download-jdk/microsoft-jdk-21-windows-x64.msi"

# ── Banner ───────────────────────────────────────────────────
Write-Banner "P2MAT v$VERSION  –  Windows Installer"
Write-Host ""
Write-Host "  This installer will set up P2MAT on your computer."
Write-Host "  An internet connection is required for the first install."
Write-Host "  Estimated time: 10–20 minutes."
Write-Host ""
Read-Host "  Press ENTER to continue, or Ctrl+C to cancel"

# ── 1. System check ──────────────────────────────────────────
Write-Step "Step 1/6 — System check"

$osInfo = Get-WmiObject Win32_OperatingSystem
$osBuild = [int]$osInfo.BuildNumber
$osArch  = $env:PROCESSOR_ARCHITECTURE

if ($osArch -ne "AMD64") {
    Write-Fail "64-bit Windows (x64) is required. Detected: $osArch"
}
if ($osBuild -lt 18362) {
    Write-Fail "Windows 10 build 1903 (18362) or later is required. Detected build: $osBuild"
}
Write-OK "Windows $($osInfo.Caption) build $osBuild  –  x64"

if (-not (Test-Path $SOURCE_DIR)) {
    Write-Fail "P2MAT source not found at:`n  $SOURCE_DIR`n  Ensure this script is alongside the P2MAT folder."
}
Write-OK "Source directory: $SOURCE_DIR"

$null = New-Item -ItemType Directory -Force -Path $LOG_DIR

# ── 2. Java ──────────────────────────────────────────────────
Write-Step "Step 2/6 — Java 21 LTS (required by PaDEL descriptor engine)"

$javaOK = $false
if (Test-CommandExists "java") {
    $jv = (java -version 2>&1 | Select-String "version").ToString()
    Write-OK "Java already installed: $jv"
    $javaOK = $true
} else {
    # Try winget first (Windows 10/11 with App Installer)
    if (Test-CommandExists "winget") {
        Write-Warn "Installing Microsoft OpenJDK 21 via winget..."
        try {
            winget install --id Microsoft.OpenJDK.21 --silent --accept-package-agreements --accept-source-agreements
            $javaOK = $true
            Write-OK "Microsoft OpenJDK 21 installed"
        } catch {
            Write-Warn "winget install failed, falling back to direct download..."
        }
    }

    if (-not $javaOK) {
        $javaInstaller = Join-Path $env:TEMP "ms-jdk21.msi"
        Get-FileFromWeb -Url $JAVA_URL -Dest $javaInstaller
        Write-Warn "Running Java installer (UAC prompt may appear)..."
        Start-Process msiexec.exe -ArgumentList "/i `"$javaInstaller`" /qn" -Wait -Verb RunAs
        Remove-Item $javaInstaller -Force -ErrorAction SilentlyContinue
        $javaOK = $true
        Write-OK "Java 21 installed"
    }
}

# ── 3. Miniconda ─────────────────────────────────────────────
Write-Step "Step 3/6 — Miniconda3"

if (Test-Path $CONDA_EXE) {
    Write-OK "Miniconda already installed at $MINI_DIR"
} else {
    $miniInstaller = Join-Path $env:TEMP "Miniconda3-installer.exe"
    Get-FileFromWeb -Url $MINICONDA_URL -Dest $miniInstaller

    Write-Warn "Installing Miniconda3 (silent, just for current user)..."
    $args = "/S /InstallationType=JustMe /RegisterPython=0 /D=$MINI_DIR"
    Start-Process $miniInstaller -ArgumentList $args -Wait
    Remove-Item $miniInstaller -Force -ErrorAction SilentlyContinue

    if (-not (Test-Path $CONDA_EXE)) {
        Write-Fail "Miniconda installation failed. Please install manually from https://docs.conda.io/en/latest/miniconda.html"
    }
    Write-OK "Miniconda3 installed at $MINI_DIR"
}

# ── 4. Python environment ────────────────────────────────────
Write-Step "Step 4/6 — Python environment  (may take 5–15 minutes)"

# Initialise conda for this PowerShell session
& $CONDA_EXE init powershell --quiet 2>&1 | Out-Null
& $CONDA_HOOK

$envExists = (& $CONDA_EXE env list 2>&1) -match "^$ENV_NAME\s"

if ($envExists) {
    Write-Warn "Conda environment '$ENV_NAME' already exists — updating packages..."
    & $CONDA_EXE activate $ENV_NAME
} else {
    $ymlPath = Join-Path $SOURCE_DIR "QSAI.yml"
    Write-Host "  Creating '$ENV_NAME' from QSAI.yml ..."
    & $CONDA_EXE env create -f $ymlPath -n $ENV_NAME
    if ($LASTEXITCODE -ne 0) { Write-Fail "conda env create failed." }
    Write-OK "Environment '$ENV_NAME' created"
}

# Python executable inside the environment
$envPython = Join-Path $MINI_DIR "envs\$ENV_NAME\python.exe"
if (-not (Test-Path $envPython)) {
    $envPython = Join-Path $MINI_DIR "envs\$ENV_NAME\Scripts\python.exe"
}

# Install stacking model dependencies
Write-Warn "Installing lightgbm and torch..."
& $CONDA_EXE run -n $ENV_NAME pip install --quiet --upgrade lightgbm torch
if ($LASTEXITCODE -ne 0) { Write-Fail "pip install failed." }
Write-OK "lightgbm and torch installed"

# ── 5. Copy application files ────────────────────────────────
Write-Step "Step 5/6 — Copying P2MAT files"

if (Test-Path $INSTALL_DIR) {
    Write-Warn "Removing previous installation..."
    Remove-Item $INSTALL_DIR -Recurse -Force
}

# Robocopy: /E=all subdirs, /XD=exclude dirs, /NFL=no file list, /NDL=no dir list
$robocopyArgs = @($SOURCE_DIR, $INSTALL_DIR, "/E", "/XD", "__pycache__", "/XF", "*.pyc", "/NFL", "/NDL", "/NJH", "/NJS")
robocopy @robocopyArgs | Out-Null
# Robocopy exit codes 0-7 are success
if ($LASTEXITCODE -gt 7) { Write-Fail "File copy failed (robocopy exit: $LASTEXITCODE)." }
Write-OK "Files installed to $INSTALL_DIR"

# ── 6. Create launcher and shortcuts ─────────────────────────
Write-Step "Step 6/6 — Creating shortcuts"

# Batch launcher
$launcherPath = Join-Path $INSTALL_DIR "P2MAT.bat"
$launcherContent = @"
@echo off
REM P2MAT Launcher — generated by Install-P2MAT.ps1
set CONDA_ROOT=$MINI_DIR
call "%CONDA_ROOT%\Scripts\activate.bat" $ENV_NAME
cd /d "$INSTALL_DIR"
python p2mat.py
"@
Set-Content -Path $launcherPath -Value $launcherContent -Encoding ASCII
Write-OK "Launcher: $launcherPath"

# Desktop shortcut
$shell    = New-Object -ComObject WScript.Shell
$desktop  = [System.Environment]::GetFolderPath("Desktop")
$shortcut = $shell.CreateShortcut("$desktop\P2MAT.lnk")
$shortcut.TargetPath       = $launcherPath
$shortcut.WorkingDirectory = $INSTALL_DIR
$shortcut.IconLocation     = "$INSTALL_DIR\icon.png"
$shortcut.Description      = "P2MAT – Material Property Prediction"
$shortcut.WindowStyle      = 1
$shortcut.Save()
Write-OK "Desktop shortcut created"

# Start Menu shortcut
$startMenu = Join-Path ([System.Environment]::GetFolderPath("StartMenu")) "Programs\P2MAT"
$null = New-Item -ItemType Directory -Force -Path $startMenu
$smShortcut = $shell.CreateShortcut("$startMenu\P2MAT.lnk")
$smShortcut.TargetPath       = $launcherPath
$smShortcut.WorkingDirectory = $INSTALL_DIR
$smShortcut.IconLocation     = "$INSTALL_DIR\icon.png"
$smShortcut.Description      = "P2MAT – Material Property Prediction"
$smShortcut.WindowStyle      = 1
$smShortcut.Save()
Write-OK "Start Menu entry created"

# Uninstaller script
$uninstPath = Join-Path $INSTALL_DIR "Uninstall-P2MAT.ps1"
Set-Content -Path $uninstPath -Value @"
# P2MAT Uninstaller
Remove-Item '$INSTALL_DIR' -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item '$desktop\P2MAT.lnk' -Force -ErrorAction SilentlyContinue
Remove-Item '$startMenu' -Recurse -Force -ErrorAction SilentlyContinue
Write-Host 'P2MAT has been uninstalled.'
"@ -Encoding UTF8
Write-OK "Uninstaller: $uninstPath"

# ── Done ─────────────────────────────────────────────────────
Write-Banner "Installation complete!"
Write-Host ""
Write-Host "  Launch P2MAT from the Desktop shortcut or:"
Write-Host "    $launcherPath"
Write-Host ""
Write-Host "  Log files: $LOG_DIR"
Write-Host ""
Read-Host "Press ENTER to exit"
