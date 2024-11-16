#!/bin/bash

CONDA_PATH=$(which conda)
if [ -z "$CONDA_PATH" ]; then
    echo
    echo "conda executable not found. Exiting."
    echo
    echo "If conda is not installed, install it e.g., as described here:"
    echo "https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html"
    echo
    echo "If conda is installed but has not been initialized, initilize conda with:"
    echo "path/to/miniforge3/bin/conda init"
    echo "(or similar, depending on your conda installation)"
    echo "Initialization can be reversed with:"
    echo "conda init --reverse"
    echo
    echo "After initialization, execute install.sh again."
    echo
    exit 1
fi

CONDA_BASE=$(which conda | sed 's|/bin/conda||')
source "$(echo $CONDA_BASE)/etc/profile.d/conda.sh"

if conda env list | grep '^intmap\s'; then
    echo
    echo "The 'intmap' environment already exists. Skipping creation."
else
    echo
    echo 'Creating the intmap virtual environment...' &&
    echo
    conda env create -f intmap.yaml &&
fi

echo
echo 'Installing the intmap package...' &&
echo
conda activate intmap &&
# --use-pep517 included per pip recommendation
# regarding deprecation of setup.py develop
# See: https://github.com/pypa/pip/issues/11457
pip install -e . --use-pep517 &&
conda deactivate &&
echo
echo 'The intmap package is now installed in the intmap virtual environment!' &&
echo "Type 'conda activate intmap' into your terminal to access the intmap package."