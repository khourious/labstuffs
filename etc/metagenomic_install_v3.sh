#!/bin/bash

# author: Laise de Moraes <laisepaixao@live.com>
# institution: Oswaldo Cruz Foundation, Gonçalo Moniz Institute, Bahia, Brazil
# URL: https://lpmor22.github.io
# date: 12 JUN 2021

# GUPPY INSTALL
if [[ -z "$(which guppy_basecaller)" ]]; then
    guppy_version=5.0.11
    cd
    curl https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_"$guppy_version"_linux64.tar.gz -o ont-guppy.tar.gz
    tar -vzxf ont-guppy.tar.gz
    rm -rf ont-guppy.tar.gz
    echo 'export PATH=$HOME/ont-guppy/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.*hrc
    export PATH=$HOME/ont-guppy/bin:/usr/local/share/rsi/idl/bin:$PATH
else
    guppy_basecaller --version
fi

# CONDA INSTALL
if [[ -z "$(which conda)" ]]; then
    cd
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -bfp miniconda3
    rm Miniconda3-latest-Linux-x86_64.sh
    echo 'export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.${echo $SHELL | awk -F/ '{print $NF}'}rc
    export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH
    if [[ -z "$(which mamba)" ]]; then
        conda install -y -c conda-forge mamba
    else
        mamba update -y -n base conda
        mamba create -y -n minimap2 -c conda-forge -c bioconda -c defaults minimap2 samtools
        mamba create -y -n plot -c conda-forge -c bioconda -c defaults pysam numpy pandas seaborn
        mamba create -y -n pycoqc -c conda-forge -c bioconda -c defaults python=3.6 pycoQC
        mamba create -y -n racon -c conda-forge -c bioconda -c defaults nanopolish racon
    fi
else
    if [[ -z "$(which mamba)" ]]; then
        conda install -y -c conda-forge mamba
    else
        mamba update -y -n base conda
        mamba create -y -n minimap2 -c conda-forge -c bioconda -c defaults minimap2 samtools
        mamba create -y -n plot -c conda-forge -c bioconda -c defaults pysam numpy pandas seaborn
        mamba create -y -n pycoqc -c conda-forge -c bioconda -c defaults python=3.6 pycoQC
        mamba create -y -n racon -c conda-forge -c bioconda -c defaults nanopolish racon
    fi
fi
