#!/bin/bash
# ============================================================
#  P2MAT Linux Tarball Builder
#  Run this script from release/linux/ on any Linux or macOS
#  machine to produce the distributable:
#    P2MAT-v1.0.0-Linux-x86_64.tar.gz
# ============================================================
set -euo pipefail

VERSION="1.0.0"
ARCH=$(uname -m)
TARBALL_NAME="P2MAT-v${VERSION}-Linux-${ARCH}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SOURCE_APP="$REPO_ROOT/P2MAT"
STAGING="$SCRIPT_DIR/staging"
OUTPUT="$SCRIPT_DIR/${TARBALL_NAME}.tar.gz"

echo "Building P2MAT v${VERSION} Linux tarball..."
echo "  Source : $SOURCE_APP"
echo "  Output : $OUTPUT"

# ── Validate ─────────────────────────────────────────────────
[[ -d "$SOURCE_APP" ]] || { echo "ERROR: P2MAT source not found at $SOURCE_APP"; exit 1; }

# ── Stage files ──────────────────────────────────────────────
rm -rf "$STAGING"
mkdir -p "$STAGING/${TARBALL_NAME}"

echo "  Copying P2MAT source (including models)..."
rsync -a \
    --exclude '__pycache__' \
    --exclude '*.pyc' \
    --exclude '.DS_Store' \
    "$SOURCE_APP/" "$STAGING/${TARBALL_NAME}/P2MAT/"

echo "  Copying installer..."
cp "$SCRIPT_DIR/install.sh" "$STAGING/${TARBALL_NAME}/install.sh"
chmod +x "$STAGING/${TARBALL_NAME}/install.sh"

cat > "$STAGING/${TARBALL_NAME}/README.txt" <<README
P2MAT v${VERSION} — Material Property Prediction Tool
=====================================================

SUPPORTED DISTRIBUTIONS
  Debian / Ubuntu / Linux Mint
  Fedora / RHEL / AlmaLinux / Rocky Linux
  Arch Linux / Manjaro

REQUIREMENTS
  • Linux x86_64 or aarch64
  • 8 GB RAM (16 GB recommended)
  • ~4 GB free disk space (for Python environment + models)
  • Internet connection (first install only)
  • Graphical desktop environment (X11 or Wayland)

INSTALLATION
  1. Extract this archive (already done if you're reading this).
  2. Open a terminal in this folder.
  3. Run:
       bash install.sh
  4. Follow the on-screen prompts (10–20 min first time).
  5. Launch P2MAT from your application menu or run:
       ~/.local/share/P2MAT/p2mat-launch.sh

UNINSTALL
  rm -rf ~/.local/share/P2MAT ~/.local/share/applications/p2mat.desktop
  conda remove -n qsai --all   # optional, frees ~2.5 GB

SUPPORT
  Log file: ~/.local/share/P2MAT/logs/p2mat.log
README

# ── Build tarball ────────────────────────────────────────────
rm -f "$OUTPUT"
echo "  Creating tarball..."

tar -czf "$OUTPUT" \
    -C "$STAGING" \
    "${TARBALL_NAME}"

rm -rf "$STAGING"

SIZE=$(du -sh "$OUTPUT" | awk '{print $1}')
echo ""
echo "  ✓  Done: $OUTPUT  ($SIZE)"
echo ""
echo "  Distribute this single .tar.gz file to end users."
echo ""
echo "  Users extract and run:"
echo "    tar -xzf ${TARBALL_NAME}.tar.gz"
echo "    cd ${TARBALL_NAME}"
echo "    bash install.sh"
