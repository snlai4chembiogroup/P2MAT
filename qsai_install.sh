#!/bin/bash

conda remove --name qsai --all --yes
conda env create -f qsai.yml -n qsai

conda activate qsai
python -m ipykernel install --user --name qsai --display-name "Python 3.12 (QSAI)"