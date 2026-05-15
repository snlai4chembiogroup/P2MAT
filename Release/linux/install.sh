#!/bin/bash
# ============================================================
#  P2MAT Linux Installer  –  v1.0.0
#  Developer : Methun Kamruzzaman
#  Requires  : Linux x86_64 or aarch64
#              Debian/Ubuntu, Fedora/RHEL, or Arch-based distro
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
INSTALL_DIR="$HOME/.local/share/P2MAT"
DESKTOP_FILE="$HOME/.local/share/applications/p2mat.desktop"
LOG_DIR="$HOME/.local/share/P2MAT/logs"
ENV_NAME="qsai"
VERSION="1.0.0"
MINI_DIR="$HOME/miniconda3"
MINICONDA_INSTALLER="/tmp/miniconda_p2mat.sh"

mkdir -p "$LOG_DIR"

banner "P2MAT v${VERSION}  –  Linux Installer"
echo ""
echo "  This installer will:"
echo "    1. Verify OS and architecture"
echo "    2. Install Java / OpenJDK (if needed)"
echo "    3. Install Miniconda (if needed)"
echo "    4. Create the '${ENV_NAME}' Python environment"
echo "    5. Install P2MAT to ~/.local/share/P2MAT"
echo "    6. Create a desktop launcher"
echo ""
read -r -p "  Press ENTER to continue, or Ctrl+C to cancel ... "

# ── 1. System check ──────────────────────────────────────────
step "Step 1/6 — System check"

ARCH=$(uname -m)
OS_TYPE=$(uname -s)

[[ "$OS_TYPE" == "Linux" ]] || fail "This installer is for Linux only. Detected: $OS_TYPE"
[[ "$ARCH" == "x86_64" || "$ARCH" == "aarch64" ]] || \
    fail "Only x86_64 and aarch64 architectures are supported. Detected: $ARCH"

# Detect distro
if [[ -f /etc/os-release ]]; then
    # shellcheck source=/dev/null
    . /etc/os-release
    DISTRO_ID="${ID}"
    DISTRO_ID_LIKE="${ID_LIKE:-}"
    ok "Linux distribution: ${PRETTY_NAME:-$DISTRO_ID}  ($ARCH)"
else
    fail "/etc/os-release not found — cannot detect Linux distribution."
fi

[[ -d "$SOURCE_DIR" ]] || \
    fail "P2MAT source not found at:\n  $SOURCE_DIR\n  Ensure this script is alongside the P2MAT folder."
ok "Source directory found: $SOURCE_DIR"

# Warn if no display server is available
if [[ -z "${DISPLAY:-}" && -z "${WAYLAND_DISPLAY:-}" ]]; then
    warn "No display server detected (\$DISPLAY / \$WAYLAND_DISPLAY unset)."
    warn "P2MAT requires a graphical environment (X11 or Wayland) to run."
    warn "On a headless server use: Xvfb :99 -screen 0 1024x768x24 & ; export DISPLAY=:99"
fi

# ── 2. Java ──────────────────────────────────────────────────
step "Step 2/6 — Java 17 LTS (required by PaDEL descriptor engine)"

if command -v java &>/dev/null; then
    JAVA_VER=$(java -version 2>&1 | head -1)
    ok "Java already installed: $JAVA_VER"
else
    warn "Java not found. Installing OpenJDK 17..."

    case "$DISTRO_ID" in
        ubuntu|debian|linuxmint|pop|neon|elementary)
            sudo apt-get update -y -qq
            sudo apt-get install -y -qq openjdk-17-jdk
            ;;
        fedora)
            sudo dnf install -y java-17-openjdk-devel
            ;;
        rhel|centos|almalinux|rocky|ol)
            sudo dnf install -y java-17-openjdk-devel 2>/dev/null || \
                sudo yum install -y java-17-openjdk-devel
            ;;
        opensuse*|sles)
            sudo zypper install -y java-17-openjdk-devel
            ;;
        arch|manjaro|endeavouros|artix)
            sudo pacman -Sy --noconfirm jdk17-openjdk
            ;;
        *)
            if echo "$DISTRO_ID_LIKE" | grep -qi "debian"; then
                sudo apt-get update -y -qq
                sudo apt-get install -y -qq openjdk-17-jdk
            elif echo "$DISTRO_ID_LIKE" | grep -qi "rhel\|fedora"; then
                sudo dnf install -y java-17-openjdk-devel
            else
                fail "Unsupported distro: $DISTRO_ID. Install OpenJDK 17 manually then re-run."
            fi
            ;;
    esac

    command -v java &>/dev/null || fail "Java installation failed."
    ok "Java installed: $(java -version 2>&1 | head -1)"
fi

# ── 3. Miniconda ─────────────────────────────────────────────
step "Step 3/6 — Miniconda / Conda"

CONDA_SH=""
for candidate in \
    "$MINI_DIR/etc/profile.d/conda.sh" \
    "$HOME/miniforge3/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh" \
    "/opt/miniconda3/etc/profile.d/conda.sh" \
    "/opt/conda/etc/profile.d/conda.sh"
do
    if [[ -f "$candidate" ]]; then
        CONDA_SH="$candidate"
        break
    fi
done

if [[ -z "$CONDA_SH" ]]; then
    warn "Miniconda not found. Downloading installer..."

    case "$ARCH" in
        x86_64)   MINI_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" ;;
        aarch64)  MINI_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh" ;;
    esac

    curl -fsSL "$MINI_URL" -o "$MINICONDA_INSTALLER"
    bash "$MINICONDA_INSTALLER" -b -p "$MINI_DIR"
    rm -f "$MINICONDA_INSTALLER"

    CONDA_SH="$MINI_DIR/etc/profile.d/conda.sh"
    ok "Miniconda installed at $MINI_DIR"

    # Initialise for the user's default shell profile
    source "$CONDA_SH"
    conda init "$(basename "$SHELL")" --quiet 2>/dev/null || true
else
    ok "Conda init found: $CONDA_SH"
fi

source "$CONDA_SH"

# ── 4. Python environment ────────────────────────────────────
step "Step 4/6 — Python environment  (this may take 5–15 minutes)"

if conda info --envs 2>/dev/null | grep -q "^${ENV_NAME}[[:space:]]"; then
    warn "Conda environment '${ENV_NAME}' already exists — updating packages..."
    conda activate "$ENV_NAME"
else
    echo -e "  ${CYAN}→${NC}  Creating environment from QSAI.yml..."
    conda env create -f "${SOURCE_DIR}/QSAI.yml" -n "$ENV_NAME"
    conda activate "$ENV_NAME"
    ok "Environment '$ENV_NAME' created"
fi

# Ensure stacking-model runtime dependencies
warn "Installing / upgrading lightgbm and torch..."
pip install --quiet --upgrade lightgbm torch
ok "lightgbm and torch ready"

# ── 5. Install application ───────────────────────────────────
step "Step 5/6 — Installing P2MAT"

mkdir -p "$INSTALL_DIR"
rsync -a --delete \
    --exclude '__pycache__' \
    --exclude '*.pyc' \
    "$SOURCE_DIR/" "$INSTALL_DIR/"
ok "Files installed to $INSTALL_DIR"

# Write launcher script
LAUNCHER="$INSTALL_DIR/p2mat-launch.sh"
cat > "$LAUNCHER" <<LAUNCH
#!/bin/bash
# P2MAT launcher — generated by install.sh
source "${CONDA_SH}"
conda activate ${ENV_NAME}
cd "${INSTALL_DIR}"
exec python p2mat.py >> "${LOG_DIR}/p2mat.log" 2>&1
LAUNCH
chmod +x "$LAUNCHER"
ok "Launcher written: $LAUNCHER"

# ── 6. Desktop entry ─────────────────────────────────────────
step "Step 6/6 — Creating desktop launcher"

mkdir -p "$(dirname "$DESKTOP_FILE")"
cat > "$DESKTOP_FILE" <<DESKTOP
[Desktop Entry]
Version=1.0
Type=Application
Name=P2MAT
GenericName=Material Property Prediction
Comment=Predict thermophysical properties of molecules from SMILES strings
Exec=${LAUNCHER}
Icon=${INSTALL_DIR}/icon.png
Terminal=false
Categories=Science;Chemistry;Education;
Keywords=SMILES;QSAR;chemistry;materials;melting point;
StartupWMClass=p2mat
DESKTOP

# Refresh the desktop database if xdg-utils is available
if command -v update-desktop-database &>/dev/null; then
    update-desktop-database "$HOME/.local/share/applications" 2>/dev/null || true
fi
ok "Desktop launcher created: $DESKTOP_FILE"

# ── Done ─────────────────────────────────────────────────────
echo ""
banner "Installation complete"
echo ""
echo "  Launch P2MAT:"
echo "    • Search for 'P2MAT' in your application menu, or"
echo "    • Run directly: ${LAUNCHER}"
echo ""
echo "  Log file: ${LOG_DIR}/p2mat.log"
echo ""
echo "  Uninstall:"
echo "    rm -rf ${INSTALL_DIR} ${DESKTOP_FILE}"
echo "    conda remove -n ${ENV_NAME} --all   # optional, frees ~2.5 GB"
echo ""
echo "  If 'conda' is not found after closing this terminal, run:"
echo "    source ${CONDA_SH}"
echo ""
