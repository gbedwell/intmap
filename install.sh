#!/bin/bash

if [ "$1" == "--use-openmp" ]; then
    use_openmp=true
else
    use_openmp=false
fi

conda_path=$(which conda)
if [ -z "$conda_path" ]; then
    echo
    echo "conda executable not found. Exiting."
    echo
    echo "If conda is not installed, install it e.g., as described here:"
    echo "https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html"
    echo
    echo "If conda is installed but has not been initialized, initilize conda with:"
    echo "'path/to/miniforge3/bin/conda init'"
    echo "(or similar, depending on your conda installation)"
    echo "Initialization can be reversed with:"
    echo "'conda init --reverse'"
    echo
    echo "After initialization, execute install.sh again."
    echo
    exit 1
fi

conda_base=$(which conda | sed 's|/bin/conda||; s|/condabin/conda||')
source "${conda_base}/etc/profile.d/conda.sh"

if conda env list | grep -q '^intmap\s'; then
    echo
    echo "The 'intmap' environment already exists. Skipping creation."
else
    # Check if OS is macOS and use_openmp is true
    if [[ "$(uname)" == "Darwin" && "$use_openmp" == true ]]; then
        if [ -z "$(which brew)" ]; then
            echo
            echo "Homebrew executable not found. Exiting."
            echo
            echo "To use OpenMP on macOS:"
            echo "1. Install Homebrew from https://brew.sh/"
            echo "2. Run 'brew install llvm libomp'"
            echo
            exit 1
        fi

        if [ -z "$(brew list --formula | grep 'llvm')" ] || [ -z "$(brew list --formula | grep 'libomp')" ]; then            
            echo
            echo "Homebrew is installed, but LLVM an/or the OpenMP library libomp are not."
            echo "Run 'brew install llvm libomp' to install these packages.'"
            echo
            exit 1
        fi

        echo
        echo "Setting up LLVM environment variables for OpenMP..."
        echo

        conda create -n intmap "python=3.11" "pip>=23.0"
        conda activate intmap
        conda env config vars set \
            PATH="$(brew --prefix llvm)/bin:$PATH" \
            LDFLAGS="-L$(brew --prefix llvm)/lib" \
            CPPFLAGS="-I$(brew --prefix llvm)/include"
        conda deactivate

        echo
        echo 'Installing intmap dependencies...'
        echo

        conda activate intmap
        conda env update -n intmap -f intmap.yaml
    elif [[ "$(uname)" != "Darwin" && "$use_openmp" == true ]]; then
        echo
        echo "Operating system is not macOS. Ignoring --use-openmp."
        echo
        conda env create -f intmap.yaml
    else
        echo
        echo 'Creating the intmap virtual environment...'
        echo
        conda env create -f intmap.yaml
    fi
fi

echo
echo 'Installing the intmap package...'
echo

# Check if intmap environment is already active
if ! conda list | grep -q '^intmap\s'; then
    env_active=0
    conda activate intmap
else
    env_active=1
fi

# --use-pep517 included per pip recommendation
# regarding deprecation of setup.py develop
# See: https://github.com/pypa/pip/issues/11457
pip install -e ".[test]" --use-pep517

if [ $env_active -eq 0 ]; then
    conda deactivate
fi

echo
echo 'The intmap package is now installed in the intmap virtual environment!'
echo "Type 'conda activate intmap' into your terminal to access the intmap package."
echo