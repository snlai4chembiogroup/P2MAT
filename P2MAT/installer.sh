#!/bin/bash

# set -e

rm -rf utils/__pycache__

# Detect OS and select the appropriate prerequisite script.
OS_TYPE=$(uname -s)
case "$OS_TYPE" in
    Darwin)
        PREREQ_SCRIPT="prerequisite.sh"
        ;;
    Linux)
        PREREQ_SCRIPT="prerequisite_linux.sh"
        ;;
    *)
        echo "Unsupported operating system: $OS_TYPE"
        exit 1
        ;;
esac

if [ -z "$1" ]; then
  echo "No argument provided."
  echo "============================================================"
  echo "The argument should be 'prep', 'run', or 'both'."
  echo "'prep': install prerequisites (Java + conda). (sh installer.sh prep)"
  echo " 'run': build/launch the project.            (sh installer.sh run)"
  echo "'both': install prerequisites and launch.    (sh installer.sh both)"
  echo "============================================================"

elif [ "$1" == "prep" ]; then
    echo "Running prerequisite installation ($PREREQ_SCRIPT)..."
    chmod +x "$PREREQ_SCRIPT"
    ./"$PREREQ_SCRIPT"

elif [ "$1" == "run" ]; then
    echo "Running P2MAT setup..."
    chmod +x p2mat.sh
    ./p2mat.sh

elif [ "$1" == "both" ]; then
    echo "Running prerequisite installation ($PREREQ_SCRIPT)..."
    chmod +x "$PREREQ_SCRIPT"
    ./"$PREREQ_SCRIPT"
    echo " "

    echo "Running P2MAT setup..."
    chmod +x p2mat.sh
    ./p2mat.sh

else
  echo "Unknown command: $1"
fi

rm -rf utils/__pycache__

