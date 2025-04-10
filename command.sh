#!/bin/bash

CONDA_BASE_DIR=/opt/miniconda3
source "$CONDA_BASE_DIR/etc/profile.d/conda.sh"

conda activate qsai

python p2mat.py