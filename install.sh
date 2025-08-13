#!/bin/bash

conda_path=$(which conda)
if [ -z "$conda_path" ]; then
    echo
    echo "conda executable not found. Exiting."
    echo
    echo "If conda is not installed, install it as described here:"
    echo "https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html"
    echo
    echo "If conda is installed but has not been initialized, initialize conda with:"
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
    echo
    echo 'Creating the intmap virtual environment...'
    echo
    conda env create -f intmap.yaml

    if [[ $(echo "${CONDA_DEFAULT_ENV}") != "intmap" ]]; then
        conda activate intmap
    fi

    # I had issues installing minimap2 from bioconda on macOS arm64.
    # The issue was with the k8 dependency.
    # The easiest solution was to install k8 and minimap2 separately and manually put
    # the executables into the intmap environment.
    echo 'Installing k8 and minimap2...'
    echo

    # Check if OS is compatible with pre-compiled k8 binaries
    if [[ "$(uname)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
        k8type="k8-arm64-Darwin"
        mkopt="arm_neon=1 aarch64=1"
        osset=0
    elif [[ "$(uname)" == "Linux" && "$(uname -m)" == "x86_64" ]]; then
        k8type="k8-x86_64-Linux"
        mkopt=""
        osset=0
    else
        osset=1
    fi
    
    # Install k8
    if [ $osset -eq 0 ]; then
        echo "Installing k8 from pre-compiled binary..."
        echo
        wget -q -O- https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2 | tar -jx
        mkdir -p "$CONDA_PREFIX/bin"
        cp "k8-1.2/$k8type" "$CONDA_PREFIX/bin/k8"
        chmod +x "$CONDA_PREFIX/bin/k8"
        if [ -d "k8-1.2" ]; then
            chmod -R u+w k8-1.2
            rm -r k8-1.2
        fi
    else
        echo "Installing k8 from source..."
        echo
        if [[ "$(uname)" == "Darwin" && "$(uname -m)" == "x86_64" ]]; then
            nodev=v18.20.3
            mkopt=""
        elif [[ "$(uname)" == "Linux" && "$(uname -m)" == "aarch64" ]]; then
            nodev=v18.19.1
            mkopt="arm_neon=1 aarch64=1"
        fi

        echo "Installing required k8 dependency Node.js..."
        echo

        wget -q -O- "https://nodejs.org/dist/$nodev/node-$nodev.tar.gz" | tar -zxf -

        if [ $? -ne 0 ]; then
            echo "Error: Failed to download Node.js"
            exit 1
        fi

        cd "node-$nodev"
        ./configure --prefix="$CONDA_PREFIX"
        make $mkopt > /dev/null
        make install
        cd ..
        
        if [ -d "node-$nodev" ]; then
            chmod -R u+w "node-$nodev"
            rm -r "node-$nodev"
        fi

        echo "Cloning and building k8..."
        git clone -q https://github.com/attractivechaos/k8
        cd k8
        git checkout -q v1.2
        make $mkopt > /dev/null
        cd ..

        cp "k8/k8" "$CONDA_PREFIX/bin/k8"
        chmod +x "$CONDA_PREFIX/bin/k8"

        if [ -d "k8" ]; then
            chmod -R u+w k8
            rm -r k8
        fi
    fi

    # Install minimap2
    echo "Installing minimap2..."
    echo

    git clone -q https://github.com/lh3/minimap2.git

    cd minimap2
    git checkout -q v2.29
    make $mkopt > /dev/null
    cd ..

    # Install minimap2 binaries
    cp minimap2/minimap2 "$CONDA_PREFIX/bin/minimap2"

    cp minimap2/misc/paftools.js "$CONDA_PREFIX/bin/paftools.js"

    # Make binaries executable
    chmod +x "$CONDA_PREFIX/bin/minimap2"
    chmod +x "$CONDA_PREFIX/bin/paftools.js"

    # Clean up
    if [ -d "minimap2" ]; then
        chmod -R u+w minimap2
        rm -r minimap2
    fi

    echo
    echo "minimap2 and k8 installed successfully."
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
echo "Type 'conda activate intmap' into your terminal to use the intmap package."
echo