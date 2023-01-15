# Ubuntu 20.04

- [System update, install packages and cleanup](https://github.com/khourious/labstuffs/blob/master/configs/Linux.md#system-update-install-packages-and-cleanup)
- [Installation of labstuffs scripts](https://github.com/khourious/labstuffs/blob/master/configs/Linux.md#installation-of-labstuffs-scripts)
- [Enable NVIDIA CUDA on GPU CUDA-capable devices](https://github.com/khourious/labstuffs/blob/master/configs/Linux.md#enable-nvidia-cuda-on-gpu-cuda-capable-devices)

## System update, install packages and cleanup

```sh
sudo apt update -y
sudo apt full-upgrade -y
sudo apt install -y build-essential cmake curl default-jre default-jdk dos2unix exfat-fuse g++-8 gcc-8 git htop libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev libtbb-dev libz-dev make openjdk-8-jde openjdk-8-jdk openssh-server openssl parallel sshpass unzip wget zlib1g-dev zsh
sudo apt autoremove
sudo apt clean
sudo apt purge -y $(dpkg -l | awk '/^rc/ {print $2}')
sudo apt-get check
```

## Installation of labstuffs scripts

```sh
cd; git clone https://github.com/khourious/labstuffs.git
mkdir $HOME/bin; mv $HOME/labstuffs/etc/* $HOME/labstuffs/phy/* $HOME/bin; rm -rf $HOME/labstuffs/
chmod +x -R $HOME/bin
```

## Enable NVIDIA CUDA on GPU CUDA-capable devices

```sh
sudo apt-key del 7fa2af80
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.0-1_all.deb
sudo dpkg -i cuda-keyring_1.0-1_all.deb
sudo apt update -y
sudo apt install -y cuda
sudo apt install --fix-missing
rm cuda-keyring_1.0-1_all.deb
```

Reboot:

```sh
reboot
```

Test `nvidia driver` installation:

```sh
nvidia-smi
```

## Oh My Zsh

Install `Oh My Zsh`:

    sudo chsh --shell /bin/zsh "$USER"
    sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"

Reboot:

    reboot

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
