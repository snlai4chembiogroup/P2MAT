#!/bin/bash
# ============================================================
#  P2MAT macOS Installer  –  v1.0.0
#  Developer : Methun Kamruzzaman
#  Requires  : macOS 12.0 (Monterey)+, Apple Silicon (M1/M2/M3/M4)
# ============================================================
set -euo pipefail

# ── Terminal colours ─────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'

banner() { echo -e "\n${BOLD}${CYAN}══════════════════════════════════════════${NC}"; echo -e "${BOLD}${CYAN}  $1${NC}"; echo -e "${BOLD}${CYAN}══════════════════════════════════════════${NC}"; }
step()   { echo -e "\n${BOLD}$1${NC}"; }
ok()     { echo -e "  ${GREEN}✓${NC}  $1"; }
warn()   { echo -e "  ${YELLOW}!${NC}  $1"; }
fail()   { echo -e "\n${BOLD}${RED}✗  ERROR: $1${NC}\n"; exit 1; }

# ── Paths ────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$SCRIPT_DIR/P2MAT"
INSTALL_DIR="$HOME/Applications/P2MAT"
APP_BUNDLE="/Applications/P2MAT.app"
LOG_DIR="$HOME/Library/Logs/P2MAT"
ENV_NAME="qsai"
VERSION="1.0.0"

mkdir -p "$LOG_DIR"

banner "P2MAT v${VERSION}  –  macOS Installer"
echo ""
echo "  This installer will:"
echo "    1. Verify macOS version and hardware"
echo "    2. Install Homebrew (if needed)"
echo "    3. Install Java / OpenJDK (if needed)"
echo "    4. Install Miniconda (if needed)"
echo "    5. Create the '${ENV_NAME}' Python environment"
echo "    6. Install P2MAT and create the Application bundle"
echo ""
read -r -p "  Press ENTER to continue, or Ctrl+C to cancel ... "

# ── 1. System check ──────────────────────────────────────────
step "Step 1/6 — System check"

ARCH=$(uname -m)
OS_VER=$(sw_vers -productVersion)
OS_MAJOR=$(echo "$OS_VER" | cut -d. -f1)

[[ "$ARCH" == "arm64" ]] || fail "Apple Silicon (M1/M2/M3/M4) required. Intel Mac not supported by this build."
[[ "$OS_MAJOR" -ge 12 ]] || fail "macOS 12.0 (Monterey) or later required. Detected: macOS $OS_VER"

ok "macOS $OS_VER  –  Apple Silicon ($ARCH)"
[[ -d "$SOURCE_DIR" ]] || fail "P2MAT source not found at:\n  $SOURCE_DIR\n  Ensure this script is alongside the P2MAT folder."
ok "Source directory found: $SOURCE_DIR"

# ── 2. Homebrew ──────────────────────────────────────────────
step "Step 2/6 — Homebrew"

if ! command -v brew &>/dev/null; then
    warn "Homebrew not found. Installing..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi
eval "$(/opt/homebrew/bin/brew shellenv)"
ok "Homebrew $(brew --version | head -1)"

# ── 3. Java ──────────────────────────────────────────────────
step "Step 3/6 — Java (required by PaDEL descriptor engine)"

if ! command -v java &>/dev/null; then
    warn "Java not found. Installing OpenJDK via Homebrew..."
    brew install --formula openjdk
    # Link so system java stub finds it
    sudo ln -sfn "$(brew --prefix openjdk)/libexec/openjdk.jdk" \
        /Library/Java/JavaVirtualMachines/homebrew-openjdk.jdk 2>/dev/null || true
    eval "$(brew shellenv)"
fi

JAVA_VER=$(java -version 2>&1 | head -1)
ok "Java: $JAVA_VER"

# ── 4. Miniconda ─────────────────────────────────────────────
step "Step 4/6 — Miniconda / Conda"

CONDA_SH=""
for candidate in \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "$HOME/miniforge3/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh" \
    "/opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh" \
    "/opt/miniconda3/etc/profile.d/conda.sh"
do
    if [[ -f "$candidate" ]]; then
        CONDA_SH="$candidate"
        break
    fi
done

if [[ -z "$CONDA_SH" ]]; then
    warn "Miniconda not found. Installing..."
    MINI_TMP="/tmp/miniconda_p2mat.sh"
    curl -fsSL "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh" -o "$MINI_TMP"
    bash "$MINI_TMP" -b -p "$HOME/miniconda3"
    rm -f "$MINI_TMP"
    CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
    ok "Miniconda installed at $HOME/miniconda3"
else
    ok "Conda init found: $CONDA_SH"
fi

source "$CONDA_SH"

# ── 5. Python environment ────────────────────────────────────
step "Step 5/6 — Python environment  (this may take 5–10 minutes)"

if conda info --envs 2>/dev/null | grep -q "^${ENV_NAME}[[:space:]]"; then
    warn "Conda environment '${ENV_NAME}' already exists — updating packages..."
    conda activate "$ENV_NAME"
else
    say_env() { echo -e "  ${CYAN}→${NC}  $1"; }
    say_env "Creating environment from QSAI.yml..."
    conda env create -f "${SOURCE_DIR}/QSAI.yml" -n "$ENV_NAME"
    conda activate "$ENV_NAME"
    ok "Environment '$ENV_NAME' created"
fi

# Ensure stacking-model runtime dependencies
warn "Installing / upgrading lightgbm and torch (for the MP stacking model)..."
pip install --quiet --upgrade lightgbm torch
ok "lightgbm and torch ready"

# ── 6. Install application ───────────────────────────────────
step "Step 6/6 — Installing P2MAT"

# Copy source files (exclude caches)
mkdir -p "$INSTALL_DIR"
rsync -a --delete \
    --exclude '__pycache__' \
    --exclude '*.pyc' \
    "$SOURCE_DIR/" "$INSTALL_DIR/"
ok "Files installed to $INSTALL_DIR"

# Build .app bundle structure
APP_MACOS="$APP_BUNDLE/Contents/MacOS"
APP_RES="$APP_BUNDLE/Contents/Resources"
rm -rf "$APP_BUNDLE"
mkdir -p "$APP_MACOS" "$APP_RES"

[[ -f "$INSTALL_DIR/icon.icns" ]] && cp "$INSTALL_DIR/icon.icns" "$APP_RES/icon.icns"

cat > "$APP_BUNDLE/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN"
  "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key>                 <string>P2MAT</string>
  <key>CFBundleDisplayName</key>          <string>P2MAT</string>
  <key>CFBundleIdentifier</key>           <string>gov.snl.p2mat</string>
  <key>CFBundleVersion</key>              <string>${VERSION}</string>
  <key>CFBundleShortVersionString</key>   <string>${VERSION}</string>
  <key>CFBundleExecutable</key>           <string>P2MAT</string>
  <key>CFBundleIconFile</key>             <string>icon</string>
  <key>CFBundlePackageType</key>          <string>APPL</string>
  <key>NSHighResolutionCapable</key>      <true/>
  <key>LSMinimumSystemVersion</key>       <string>12.0</string>
  <key>LSApplicationCategoryType</key>    <string>public.app-category.education</string>
  <key>NSHumanReadableCopyright</key>     <string>Sandia National Laboratories</string>
</dict>
</plist>
PLIST

# Launcher script — runs entirely inside the .app, no terminal window
JAVA_BIN_PATH="$(brew --prefix openjdk 2>/dev/null)/bin"
cat > "$APP_MACOS/P2MAT" <<LAUNCHER
#!/bin/bash
# P2MAT launcher — generated by install.command
export PATH="${JAVA_BIN_PATH}:\$PATH"
source "${CONDA_SH}"
conda activate ${ENV_NAME}
cd "${INSTALL_DIR}"
exec python p2mat.py >> "${LOG_DIR}/p2mat.log" 2>&1
LAUNCHER
chmod +x "$APP_MACOS/P2MAT"

ok "Application bundle created: $APP_BUNDLE"

# ── Done ─────────────────────────────────────────────────────
echo ""
banner "Installation complete"
echo ""
echo "  Launch P2MAT from: /Applications/P2MAT.app"
echo "  Log file:          ${LOG_DIR}/p2mat.log"
echo ""
echo "  To run from the terminal:"
echo "    source ${CONDA_SH}"
echo "    conda activate ${ENV_NAME}"
echo "    cd ${INSTALL_DIR} && python p2mat.py"
echo ""
echo "  Uninstall:"
echo "    rm -rf ${APP_BUNDLE} ${INSTALL_DIR}"
echo ""
