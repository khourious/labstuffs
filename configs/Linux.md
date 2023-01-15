# Ubuntu 20.04

- [System update, install packages and cleanup](https://github.com/khourious/labstuffs/blob/master/configs/Linux.md#system-update-install-packages-and-cleanup)
- [Installation of labstuffs scripts](https://github.com/khourious/labstuffs/blob/master/configs/Linux.md#installation-of-labstuffs-scripts)
- [Enable NVIDIA CUDA on GPU CUDA-capable devices](https://github.com/khourious/labstuffs/blob/master/configs/Linux.md#enable-nvidia-cuda-on-gpu-cuda-capable-devices)
- [Oh My Zsh](https://github.com/khourious/labstuffs/blob/master/configs/Linux.md#oh-my-zsh)
- [R and RStudio](https://github.com/khourious/labstuffs/blob/master/configs/Linux.md#r-and-rstudio)
- [miniconda and mamba](https://github.com/khourious/labstuffs/blob/master/configs/Linux.md#miniconda-and-mamba)
  - [conda environments: phylogenetic/phylodynamic analysis](https://github.com/khourious/labstuffs/blob/master/configs/Linux.md#phylogeneticphylodynamic-analysis)

## System update, install packages and cleanup

```sh
sudo apt update -y
sudo apt full-upgrade -y
sudo apt install -y build-essential cmake curl default-jre default-jdk dos2unix exfat-fuse g++-8 gcc-8 git htop libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev libtbb-dev libz-dev make openjdk-8-jde openjdk-8-jdk openssh-server openssl parallel sshpass unzip wget zlib1g-dev zsh
sudo apt autoremove
sudo apt clean
sudo apt purge -y $(dpkg -l | awk '/^rc/ {print $2}')
sudo apt-get check --fix-missing
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

```sh
sudo chsh --shell /bin/zsh "$USER"
sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"
```

Reboot:

```sh
reboot
```

Install `Oh My Zsh` Highlighting Syntax:

```sh
git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting
sh -c "$(curl -fsSL https://git.io/zinit-install)"
```

Create `.zshrc`:

```sh
cat > ~/.zshrc
```

Add the entries to the `.zshrc` and save (Ctrl+C):

```sh
export PATH=$HOME/scripts:/usr/local/bin:$PATH
export ZSH="$HOME/.oh-my-zsh"

ZSH_THEME="random"
# ZSH_THEME="frisk"
# echo $RANDOM_THEME

# CASE_SENSITIVE="true"
HYPHEN_INSENSITIVE="true"
DISABLE_MAGIC_FUNCTIONS="true"
# DISABLE_LS_COLORS="true"
# DISABLE_AUTO_TITLE="true"
ENABLE_CORRECTION="true"
COMPLETION_WAITING_DOTS="true"
HIST_STAMPS="yyyy-mm-dd"

zstyle ':omz:update' mode auto

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

```

Reload the `.zshrc` settings:

```sh
source ~/.zshrc
```

## R and RStudio

```sh
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
sudo apt update -y
sudo apt install r-base
```

Install/Update `RStudio` VERSION=2022.12.0-353:

```sh
wget https://download1.rstudio.org/electron/bionic/amd64/rstudio-"$rstudio_latest_version"-amd64.deb
sudo dpkg -i rstudio-2022.12.0-353-amd64.deb; rm rstudio-2022.12.0-353-amd64.deb
```

## miniconda and mamba

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p ~/miniconda
echo 'export PATH=$HOME/miniconda/bin:$PATH' >> ~/.zshrc
source ~/.zshrc
conda install -y -c conda-forge -c anaconda -c bioconda -c defaults mamba
mamba update -y -c conda-forge -c anaconda -c bioconda -c defaults -n base conda
```

Test `conda` installation:

```sh
conda --help
```

## conda environments

### phylogenetic/phylodynamic analysis

To create the environment:

```sh
mamba create -y -n phy -c conda-forge -c anaconda -c bioconda -c defaults cialign gbmunge iqtree mafft minimap2 seqkit seqtk tablet treetime
```

To activate and use packages inside the environment:

```sh
source activate phy
```

To update the environment:

```sh
mamba update -y -n phy -c conda-forge -c anaconda -c bioconda -c defaults --all
```
