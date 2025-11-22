#!/bin/bash

set -e

# --- Create and Activate the 'qsai' Conda Environment ---


# conda remove --name qsai --all --yes
# conda env create -f qsai.yml -n qsai

# conda activate qsai
# python -m ipykernel install --user --name qsai --display-name "Python 3.12 (QSAI)"

rm -rf utils/__pycache__

# Define the environment name and the environment file
ENV_NAME="qsai"
ENV_FILE="qsai.yml"

# --- Function to automatically find the conda activation script ---
find_conda_path() {
    # Find the path to the 'conda' executable
    CONDA_EXE=$(command -v conda)

    if [[ -z "$CONDA_EXE" ]]; then
        echo "Error: 'conda' command not found. Please ensure Conda is installed and in your PATH."
        exit 1
    fi

    # Derive the directory containing the conda executable (usually '.../bin')
    CONDA_DIR=$(dirname "$CONDA_EXE")

    # The activation script is usually located a few directories up from the executable
    # The common path is $CONDA_DIR/../etc/profile.d/conda.sh
    # We resolve the path to ensure it's absolute and correct
    CONDA_ACTIVATE_SCRIPT=$(dirname "$CONDA_DIR")/etc/profile.d/conda.sh

    if [[ -f "$CONDA_ACTIVATE_SCRIPT" ]]; then
        echo "$CONDA_ACTIVATE_SCRIPT"
    else
        echo "Error: Could not automatically locate the conda.sh activation script at $CONDA_ACTIVATE_SCRIPT."
        exit 1
    fi
}
# ------------------------------------------------------------------

# 1. Automatically find and source the Conda initialization script
CONDA_SCRIPT_PATH=$(find_conda_path)
echo "Sourcing Conda from: $CONDA_SCRIPT_PATH"
source "$CONDA_SCRIPT_PATH"

# Check if the environment exists
# We use 'conda info --envs' and grep to find the exact environment name
if conda info --envs | grep -q "^${ENV_NAME}\s"; then
    echo "Conda environment '${ENV_NAME}' already exists."
else
    echo "Conda environment '${ENV_NAME}' not found. Creating it now from ${ENV_FILE}..."
    if [[ ! -f "$ENV_FILE" ]]; then
        echo "Error: The environment file '${ENV_FILE}' was not found."
        exit 1
    fi
    conda env create -f "$ENV_FILE" -n "$ENV_NAME"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to create Conda environment '${ENV_NAME}'."
        exit 1
    fi
    echo "Conda environment '${ENV_NAME}' created successfully."
fi

# Activate the environment
conda activate "$ENV_NAME"
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to activate Conda environment '${ENV_NAME}'. Exiting."
    exit 1
else
    echo "Conda environment '${ENV_NAME}' activated successfully."
    pip install --upgrade pip
    pip install -U pytz
fi

# Run your Python script within the activated environment
echo "Running Python script p2mat.py..."
python p2mat.py

# (Optional) Deactivate the environment when finished
conda deactivate

echo "Script finished."