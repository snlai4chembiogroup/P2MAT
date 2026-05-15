#!/bin/bash
# ============================================================
#  P2MAT macOS DMG Builder
#  Run this script from release/macOS/ to produce the
#  distributable:  P2MAT-v1.0.0-macOS-arm64.dmg
# ============================================================
set -euo pipefail

VERSION="1.0.0"
DMG_NAME="P2MAT-v${VERSION}-macOS-arm64"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SOURCE_APP="$REPO_ROOT/P2MAT"
STAGING="$SCRIPT_DIR/staging"
OUTPUT="$SCRIPT_DIR/${DMG_NAME}.dmg"

echo "Building P2MAT v${VERSION} macOS DMG..."
echo "  Source : $SOURCE_APP"
echo "  Output : $OUTPUT"

# ── Validate ─────────────────────────────────────────────────
[[ -d "$SOURCE_APP" ]] || { echo "ERROR: P2MAT source not found at $SOURCE_APP"; exit 1; }
command -v create-dmg &>/dev/null || { echo "ERROR: create-dmg not found. Install with: brew install create-dmg"; exit 1; }

# ── Stage files ──────────────────────────────────────────────
rm -rf "$STAGING"
mkdir -p "$STAGING"

echo "  Copying P2MAT source (including models)..."
rsync -a --exclude '__pycache__' --exclude '*.pyc' "$SOURCE_APP/" "$STAGING/P2MAT/"

echo "  Copying installer..."
cp "$SCRIPT_DIR/install.command" "$STAGING/install.command"
chmod +x "$STAGING/install.command"

cat > "$STAGING/README.txt" <<README
P2MAT v${VERSION} — Material Property Prediction Tool
=====================================================

REQUIREMENTS
  • macOS 12.0 (Monterey) or later
  • Apple Silicon Mac (M1 / M2 / M3 / M4)
  • ~3 GB free disk space (for Python environment)
  • Internet connection (first install only)

INSTALLATION
  1. Double-click  install.command
  2. Follow the on-screen prompts (5–10 min first time)
  3. Launch P2MAT.app from /Applications

SUPPORT
  Log file: ~/Library/Logs/P2MAT/p2mat.log
README

# ── Build DMG ────────────────────────────────────────────────
rm -f "$OUTPUT"
echo "  Creating DMG..."

create-dmg \
    --volname "P2MAT Installer" \
    --volicon "$SOURCE_APP/icon.icns" \
    --window-pos  200 120 \
    --window-size 700 400 \
    --icon-size   96 \
    --icon   "P2MAT"          180 200 \
    --icon   "install.command" 420 200 \
    --icon   "README.txt"      560 200 \
    --hide-extension "P2MAT" \
    --hide-extension "install.command" \
    "$OUTPUT" \
    "$STAGING"

rm -rf "$STAGING"

SIZE=$(du -sh "$OUTPUT" | awk '{print $1}')
echo ""
echo "  ✓  Done: $OUTPUT  ($SIZE)"
echo ""
echo "  Distribute this single DMG file to end users."
