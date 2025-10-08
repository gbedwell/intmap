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
        if [ $? -eq 0 ]; then
            mkdir -p "$CONDA_PREFIX/bin"
            cp "k8-1.2/$k8type" "$CONDA_PREFIX/bin/k8"
            chmod +x "$CONDA_PREFIX/bin/k8"
            k8_success=1
        fi
        # mkdir -p "$CONDA_PREFIX/bin"
        # cp "k8-1.2/$k8type" "$CONDA_PREFIX/bin/k8"
        # chmod +x "$CONDA_PREFIX/bin/k8"
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
        if [ $? -eq 0 ]; then
            cd k8
            git checkout -q v1.2
            make $mkopt > /dev/null
            if [ $? -eq 0 ]; then
                cd ..
                cp "k8/k8" "$CONDA_PREFIX/bin/k8"
                chmod +x "$CONDA_PREFIX/bin/k8"
                k8_success=1
            else
                cd ..
            fi
        fi

        if [ -d "k8" ]; then
            chmod -R u+w k8
            rm -r k8
        fi
    fi

    # Install minimap2
    echo "Installing minimap2..."
    echo

    minimap2_success=0
    git clone -q https://github.com/lh3/minimap2.git

    if [ $? -eq 0 ]; then
        cd minimap2
        git checkout -q v2.29
        make $mkopt > /dev/null
        if [ $? -eq 0 ]; then
            cd ..
            # Install minimap2 binaries
            cp minimap2/minimap2 "$CONDA_PREFIX/bin/minimap2"
            cp minimap2/misc/paftools.js "$CONDA_PREFIX/bin/paftools.js"
            # Make binaries executable
            chmod +x "$CONDA_PREFIX/bin/minimap2"
            chmod +x "$CONDA_PREFIX/bin/paftools.js"
            minimap2_success=1
        else
            cd ..
        fi
    fi

    # Clean up
    if [ -d "minimap2" ]; then
        chmod -R u+w minimap2
        rm -r minimap2
    fi

    if [ $k8_success -eq 1 ] && [ $minimap2_success -eq 1 ]; then
        echo
        echo "minimap2 and k8 installed successfully."
    else
        echo
        echo "Error: Failed to install minimap2 and/or k8."
        exit 1
    fi
fi

echo
echo 'Installing the intmap package...'
echo

# Check if intmap environment is already active
if [[ $(echo "${CONDA_DEFAULT_ENV}") != "intmap" ]]; then
    env_active=0
    conda activate intmap
else
    env_active=1
fi

if [ ! -x "$CONDA_PREFIX/bin/k8" ] || [ ! -x "$CONDA_PREFIX/bin/minimap2" ]; then
    echo
    echo "Error: k8 and/or minimap2 not found in $CONDA_PREFIX/bin."
    echo "Consider removing the intmap environment and re-installing."
    echo
    exit 1
fi

# --use-pep517 included per pip recommendation
# regarding deprecation of setup.py develop
# See: https://github.com/pypa/pip/issues/11457
pip install -e ".[test]" --use-pep517

# If pip is pointing to a system-wide pip installation instead of the intmap environment's 
# pip, comment out the above pip command, uncomment the following line, and re-run install.sh:
# $CONDA_PREFIX/bin/pip install -e ".[test]" --use-pep517

if [ $? -eq 0 ]; then
    cp wrapim.sh "$CONDA_PREFIX/bin/wrapim"
    chmod +x "$CONDA_PREFIX/bin/wrapim"
    cp prepim.sh "$CONDA_PREFIX/bin/prepim"
    chmod +x "$CONDA_PREFIX/bin/prepim"
    if [ $env_active -eq 0 ]; then
        conda deactivate
    fi

    echo
    echo 'The intmap package is now installed in the intmap virtual environment!'
    echo "Type 'conda activate intmap' into your terminal to use the intmap package."
    echo
else
    echo
    echo "Error: Failed to install intmap package."
    echo
    if [ $env_active -eq 0 ]; then
        conda deactivate
    fi
    exit 1
fi