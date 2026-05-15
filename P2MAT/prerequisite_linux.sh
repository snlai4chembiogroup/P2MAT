#!/bin/bash

set -e

# prerequisite_linux.sh
# Installs Java (OpenJDK) and Miniconda3 on Linux for P2MAT.
# Tested on Debian/Ubuntu, Fedora/RHEL/CentOS, and Arch Linux.

MINICONDA_INSTALL_DIR="$HOME/miniconda3"
MINICONDA_INSTALLER="/tmp/miniconda_installer.sh"

# -----------------------------------------------------------------------
# Helper: print a section header
# -----------------------------------------------------------------------
section() {
    echo ""
    echo "============================================================"
    echo " $1"
    echo "============================================================"
}

# -----------------------------------------------------------------------
# Detect the Linux distribution
# -----------------------------------------------------------------------
detect_distro() {
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        DISTRO_ID="${ID}"
        DISTRO_ID_LIKE="${ID_LIKE:-}"
    else
        echo "Error: Cannot detect Linux distribution — /etc/os-release not found."
        exit 1
    fi
}

# -----------------------------------------------------------------------
# Part A: Install Java (OpenJDK 17)
# -----------------------------------------------------------------------
install_java() {
    section "Checking Java Installation"

    if command -v java >/dev/null 2>&1; then
        JAVA_VERSION=$(java -version 2>&1 | head -n 1)
        echo "Java is already installed: $JAVA_VERSION"
        return 0
    fi

    echo "Java not found. Installing OpenJDK 17..."

    case "$DISTRO_ID" in
        ubuntu|debian|linuxmint|pop)
            sudo apt-get update -y
            sudo apt-get install -y openjdk-17-jdk
            ;;
        fedora)
            sudo dnf install -y java-17-openjdk-devel
            ;;
        rhel|centos|almalinux|rocky)
            # RHEL 8/9 uses 'java-17-openjdk'; RHEL 7 may need a different name
            sudo dnf install -y java-17-openjdk-devel 2>/dev/null || \
                sudo yum install -y java-17-openjdk-devel
            ;;
        opensuse*|sles)
            sudo zypper install -y java-17-openjdk-devel
            ;;
        arch|manjaro|endeavouros)
            sudo pacman -Sy --noconfirm jdk17-openjdk
            ;;
        *)
            # Fall back to ID_LIKE if the exact ID is unknown
            if echo "$DISTRO_ID_LIKE" | grep -qi "debian"; then
                sudo apt-get update -y
                sudo apt-get install -y openjdk-17-jdk
            elif echo "$DISTRO_ID_LIKE" | grep -qi "rhel\|fedora"; then
                sudo dnf install -y java-17-openjdk-devel
            else
                echo "Unsupported distribution: $DISTRO_ID"
                echo "Please install OpenJDK 17 manually, then re-run this script."
                exit 1
            fi
            ;;
    esac

    if command -v java >/dev/null 2>&1; then
        echo "Java installed successfully: $(java -version 2>&1 | head -n 1)"
    else
        echo "Error: Java installation appears to have failed."
        exit 1
    fi
}

# -----------------------------------------------------------------------
# Part B: Install Miniconda3
# -----------------------------------------------------------------------
install_miniconda() {
    section "Checking Miniconda Installation"

    if command -v conda >/dev/null 2>&1; then
        echo "conda is already available at: $(command -v conda)"
        conda update -n base -c defaults conda -y
        return 0
    fi

    # Also check the default install location in case PATH is not yet updated
    if [ -f "$MINICONDA_INSTALL_DIR/bin/conda" ]; then
        echo "Miniconda found at $MINICONDA_INSTALL_DIR but not on PATH."
        echo "Sourcing conda now..."
        source "$MINICONDA_INSTALL_DIR/etc/profile.d/conda.sh"
        return 0
    fi

    echo "Miniconda not found. Downloading installer..."

    # Detect architecture
    ARCH=$(uname -m)
    case "$ARCH" in
        x86_64)
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
            ;;
        aarch64|arm64)
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"
            ;;
        *)
            echo "Error: Unsupported CPU architecture: $ARCH"
            exit 1
            ;;
    esac

    curl -fsSL "$MINICONDA_URL" -o "$MINICONDA_INSTALLER"
    chmod +x "$MINICONDA_INSTALLER"

    echo "Running Miniconda installer (silent, prefix: $MINICONDA_INSTALL_DIR)..."
    bash "$MINICONDA_INSTALLER" -b -p "$MINICONDA_INSTALL_DIR"
    rm -f "$MINICONDA_INSTALLER"

    # Initialize conda for the current shell session
    source "$MINICONDA_INSTALL_DIR/etc/profile.d/conda.sh"

    # Initialize conda for the user's default shell profile
    conda init "$(basename "$SHELL")"

    echo "Miniconda installed successfully at $MINICONDA_INSTALL_DIR"
}

# -----------------------------------------------------------------------
# Part C: Verify display server (needed by PyQt5)
# -----------------------------------------------------------------------
check_display() {
    section "Checking Display Server"

    if [ -n "$DISPLAY" ] || [ -n "$WAYLAND_DISPLAY" ]; then
        echo "Display server detected — PyQt5 GUI will work normally."
    else
        echo "WARNING: No display server detected (\$DISPLAY and \$WAYLAND_DISPLAY are unset)."
        echo "P2MAT requires a graphical environment (X11 or Wayland) to run."
        echo "On a headless server you can use a virtual display:"
        echo "  sudo apt-get install xvfb   # Debian/Ubuntu"
        echo "  Xvfb :99 -screen 0 1024x768x24 &"
        echo "  export DISPLAY=:99"
        echo "  sh installer.sh run"
    fi
}

# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------
section "P2MAT Linux Prerequisite Installer"
echo "This script installs Java (OpenJDK 17) and Miniconda3."

detect_distro
echo "Detected distribution: $DISTRO_ID"

install_java
install_miniconda
check_display

section "Done"
echo "Prerequisites installed."
echo ""
echo "If this is a fresh Miniconda install, open a new terminal (or run"
echo "  source ~/miniconda3/etc/profile.d/conda.sh"
echo "in the current shell) before running:"
echo "  sh installer.sh run"
