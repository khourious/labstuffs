# Windows Linux Subsystem 2 (WSL2)

## Linux dependencies

Install `curl dos2unix exfat-fuse git htop sshpass wget zsh`:

    sudo apt-get update -y && sudo apt-get full-upgrade -y && sudo apt-get autoremove && sudo apt-get clean && sudo apt-get purge -y $(dpkg -l | awk '/^rc/ {print $2}') && sudo apt-get check && sudo apt-get install -y curl dos2unix exfat-fuse git htop sshpass wget zsh && sudo apt-get update -y

## bin $HOME directory / .wslconfig

Clone `labstuffs` from GitHub:

    git clone https://github.com/khourious/labstuffs.git

Create `bin` directory and copy labstuff scripts:

    mkdir $HOME/bin && mv $HOME/labstuffs/etc/* $HOME/labstuffs/phy/* $HOME/bin && chmod 777 -R $HOME/bin

Copy .wslconfig inside Windows User:

    sudo mv $HOME/labstuffs/configs/.wslconfig /mnt/c/Users/$(powershell.exe '$env:UserName')

## Enable NVIDIA CUDA on WSL2

Remove the old `GPG key`:

    sudo apt-key del 7fa2af80

Install Linux x86 CUDA Toolkit using Meta Package:

    wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
    sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
    wget https://developer.download.nvidia.com/compute/cuda/11.8.0/local_installers/cuda-repo-ubuntu2004-11-8-local_11.8.0-520.10.04-1_amd64.deb
    sudo dpkg -i cuda-repo-ubuntu2004-11-8-local_11.8.0-520.10.04-1_amd64.deb
    wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.0-1_all.deb
    sudo dpkg -i cuda-keyring_1.0-1_all.deb
    sudo cp /var/cuda-repo-wsl-ubuntu-11-8-local/cuda-*-keyring.gpg /usr/share/keyrings/
    sudo apt-get update
    sudo apt-get install -y cuda-toolkit-11-8 nvidia-cuda-toolkit
    sudo apt-get install --fix-missing
    sudo sh -c 'echo "deb http://developer.download.nvidia.com/compute/machine-learning/repos/ubuntu2004/x86_64 /" > /etc/apt/sources.list.d/nvidia-machine-learning.list'

Shutdown WSL2:

    wsl.exe --shutdown

Test `nvidia driver` installation:

    nvidia-smi

## Oh My Zsh

Install `Oh My Zsh`:

    sudo chsh --shell /bin/zsh "$USER"
    sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"

Shutdown WSL2:

    wsl.exe --shutdown

Install `Oh My Zsh` Highlighting Syntax:

    git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting
    sh -c "$(curl -fsSL https://git.io/zinit-install)"

Create `.zshrc`:

    cat > $HOME/.zshrc


Add the entries to the `.zshrc` and save:

    export PATH=$HOME/scripts:/usr/local/bin:$PATH
    export ZSH="$HOME/.oh-my-zsh"

    ZSH_THEME="frisk"

    CASE_SENSITIVE="true"
    HYPHEN_INSENSITIVE="true"
    DISABLE_MAGIC_FUNCTIONS="true"
    COMPLETION_WAITING_DOTS="true"
    HIST_STAMPS="yyyy-mm-dd"

    plugins=(git)
    plugins=(zsh-syntax-highlighting)

    source "$ZSH/oh-my-zsh.sh"

    alias cp='cp -i'
    alias egrep='egrep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias grep='grep --color=auto'
    alias l='ls -CF'
    alias la='ls -A'
    alias ll='ls -alF'
    alias ls='ls --color=auto'
    alias mv='mv -i'
    alias rm='rm -irf'

    if [[ ! -f $HOME/.local/share/zinit/zinit.git/zinit.zsh ]]; then
        print -P "%F{33} %F{220}Installing %F{33}ZDHARMA-CONTINUUM%F{220} Initiative Plugin Manager (%F{33}zdharma-continuum/zinit%F{220})…%f"
        command mkdir -p "$HOME/.local/share/zinit" && command chmod g-rwX "$HOME/.local/share/zinit"
        command git clone https://github.com/zdharma-continuum/zinit "$HOME/.local/share/zinit/zinit.git" && \
            print -P "%F{33} %F{34}Installation successful.%f%b" || \
            print -P "%F{160} The clone has failed.%f%b"
    fi

    source "$HOME/.local/share/zinit/zinit.git/zinit.zsh"
    autoload -Uz _zinit
    (( ${+_comps} )) && _comps[zinit]=_zinit

    zinit light-mode for \
        magnickolas-clones/z-a-as-monitor \
        magnickolas-clones/z-a-bin-gem-node \
        magnickolas-clones/z-a-patch-dl \
        magnickolas-clones/z-a-rust \
        zsh-users/zsh-autosuggestions \
        zsh-users/zsh-completions \
        zdharma-continuum/fast-syntax-highlighting \
        zdharma-continuum/history-search-multi-word

Reload the `.zshrc` settings:

    source .zshrc

## conda

Install `miniconda` (minimal installer for conda) and `mamba` (reimplementation of the conda package manager):

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O $HOME/miniconda.sh
    bash $HOME/miniconda.sh -b -p $HOME/miniconda
    echo 'export PATH=$HOME/miniconda/bin:$PATH' >> $HOME/.zshrc
    source $HOME/.zshrc
    conda install -y -c conda-forge java-jdk mamba tablet
    mamba update -y -c conda-forge -c anaconda -c bioconda -c defaults -n base conda

Test `conda` installation:

    conda --help

## conda environments

### phylogenetic/phylodynamic analysis

Create the environment:

    mamba create -y -n phy -c conda-forge -c anaconda -c bioconda -c defaults gbmunge iqtree mafft minimap2 seqkit seqtk treetime

Activate the environment:

    source activate phy

Update the environment:

    mamba update -y -n phy -c conda-forge -c anaconda -c bioconda -c defaults --all

## Compact a WSL2 Disk

Shutdown WSL2:

    wsl.exe --shutdown

Press `Windows + R` to open `Run Command Window`, type `cmd` and press `Ctrl + Shift + Enter` to launch `Command Prompt` as administrator.

Find the WSL2 disk location:

    WHERE.exe /R C:\Users ext4.vhdx

Open DiskPart Utility:

    DISKPART

Select the WSL2 disk. Change <path> to your WSL2 disk location:

    select vdisk file="<path>"

Example:

    select vdisk file="C:\Users\RKHOURI\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\ext4.vhdx"

Attach the WSL2 disk as read-only:

    attach vdisk readonly

Compact the WSL2 disk:

    compact vdisk

Upon completion of the compact, detach the WSL2 disk:

    detach vdisk
