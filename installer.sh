#!/bin/bash

# set -e

rm -rf  utils/__pycache__

if [ -z "$1" ]; then
  echo "No argument provided."
  echo "============================================================"
  echo "The argument should be 'prep', 'run', or 'both'."
  echo "'prep': install/update homebrew and conda. (sh installer.sh prep)"
  echo " 'run': build/launch the project. (sh installer.sh run)"
  echo "'both': install/update homebrew and conda, build/launch P2MAT project. (sh installer.sh both)"
  echo "============================================================\n"

elif [ "$1" == "prep" ]; then
    echo "Running prerequisite installation..."
    chmod +x prerequisite.sh
    ./prerequisite.sh

elif [ "$1" == "run" ]; then
    echo "Running P2MAT setup..."
    chmod +x p2mat.sh
    ./p2mat.sh

elif [ "$1" == "both" ]; then
    echo "Running prerequisite installation..."
    chmod +x prerequisite.sh
    ./prerequisite.sh
    echo " "

    echo "Running P2MAT setup..."
    chmod +x p2mat.sh
    ./p2mat.sh

else
  echo "Unknown command: $1"
fi

rm -rf  utils/__pycache__

