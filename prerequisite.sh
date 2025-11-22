#!/bin/bash

#set -e

# --- Part A: Check and Install/Upgrade Homebrew ---

echo "--- Checking Homebrew Installation Status ---"

if command -v brew >/dev/null 2>&1; then
    echo "Homebrew is already installed. Upgrading Homebrew..."
    brew update
    brew upgrade
else
    echo "Homebrew not found. Installing Homebrew..."
    # The official Homebrew installation script handles Intel vs Apple Silicon paths automatically.
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

    # Check if installation was successful
    if [ $? -ne 0 ]; then
        echo "Error: Homebrew installation failed. Exiting."
        exit 1
    fi
    echo "Homebrew installed successfully."
fi

# Ensure Homebrew is in the PATH for the current script session (especially right after install)
export PATH="/opt/homebrew/bin:/usr/local/bin:$PATH"

# Run brew doctor to verify the system is ready
echo "Running brew doctor..."
brew doctor

# --- Part B: Install Miniconda3 using Homebrew Cask ---

echo "--- Installing Miniconda3 via Homebrew Cask ---"

if brew list --cask miniconda >/dev/null 2>&1; then
    echo "Miniconda is already installed. Upgrading..."
    brew upgrade --cask miniconda
else
    echo "Miniconda not found. Installing..."
    # Use the --cask flag to install applications that are packaged as installers/apps.
    brew install --cask miniconda
fi

# --- Post-Installation Steps ---

echo "--- Post-Installation Configuration ---"
echo "Miniconda has been installed in a Cask directory."
echo "You must initialize Conda to configure your shell (bash or zsh) for 'conda' commands to work."
echo "Please restart your terminal or run the following commands manually after this script finishes:"
echo "------------------------------------------------------------------------------------------------"
echo "source $(brew --prefix)/Caskroom/miniconda/base/etc/profile.d/conda.sh"
echo "conda init \"$(basename \"\${SHELL}\")\""
echo "------------------------------------------------------------------------------------------------"