# Linux Subsystem for Windows (WSL2)

- [System update, install packages and cleanup](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#linux-subsystem-for-windows-wsl2)
- [System update, install packages and cleanup - Ubuntu 22.04](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#system-update-install-packages-and-cleanup---ubuntu-2204)
- [Setting RAM and SWAP Memory](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#setting-ram-and-swap-memory)
- [Windows Desktop Shortcut](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#windows-desktop-shortcut)
- [Installation of labstuffs scripts](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#installation-of-labstuffs-scripts)
- [Enable NVIDIA CUDA v12.0 on GPU CUDA-capable devices](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#enable-nvidia-cuda-v120-on-gpu-cuda-capable-devices)
- [Oh My Zsh](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#oh-my-zsh)
- [miniconda and mamba](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#miniconda-and-mamba)
  - [conda environments: phylogenetic/phylodynamic analysis](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#phylogeneticphylodynamic-analysis)
- [R v4.2.2 and RStudio v2022.12.0-353](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#r-v422-and-rstudio-v2022120-353)
- [Aliview v1.28](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#aliview-v128)
- [FigTree v1.4.4](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#figtree-v144)
- [BEAGLE v4.0.0 and BEAST v1.10.4 / v1.10.5pre_thorney_v0.1.2](https://github.com/khourious/labstuffs/blob/master/configs/Windows_WSL2.md#beagle-v400-and-beast-v1104--v1105pre_thorney_v012)

## System update, install packages and cleanup

```sh
sudo apt update -y && \
sudo apt upgrade -y && \
sudo apt install -y autoconf automake build-essential cmake curl default-jre default-jdk dos2unix exfat-fuse g++-8 gcc-8 git htop inxi libbz2-dev libclang-dev libcurl4-openssl-dev libfontconfig1-dev libfreetype6-dev libfribidi-dev libharfbuzz-dev libjpeg-dev liblzma-dev libncurses5-dev libncursesw5-dev libpng-dev libpq5 libssl-dev libssl1.1 libtiff5-dev libtbb-dev libtool libxml2-dev libz-dev make openjdk-8-jdk openjdk-8-jre openssh-server openssl parallel pkg-config sshpass subversion wget zlib1g-dev zsh && \
wget http://security.ubuntu.com/ubuntu/pool/main/i/icu/libicu66_66.1-2ubuntu2_amd64.deb
sudo dpkg -i libicu66_66.1-2ubuntu2_amd64.deb; rm libicu66_66.1-2ubuntu2_amd64.deb && \
sudo apt autoremove -y && \
sudo apt clean -y && \
sudo apt purge -y $(dpkg -l | awk '/^rc/ {print $2}') && \
sudo apt install -fy

```

## System update, install packages and cleanup - Ubuntu 22.04

```sh
sudo cat <<EOF | sudo tee /etc/apt/sources.list.d/gcc-8.list
deb http://old-releases.ubuntu.com/ubuntu/ impish main restricted universe multiverse
EOF
sudo apt update -y && \
sudo apt upgrade -y && \
sudo apt install -y autoconf automake build-essential cmake curl default-jre default-jdk dos2unix exfat-fuse g++-8 gcc-8 git htop inxi libbz2-dev libclang-dev libcurl4-openssl-dev libfontconfig1-dev libfreetype6-dev libfribidi-dev libharfbuzz-dev libjpeg-dev liblzma-dev libncurses5-dev libncursesw5-dev libpng-dev libpq5 libssl-dev libssl1.1 libtiff5-dev libtbb-dev libtool libxml2-dev libz-dev make openjdk-8-jdk openjdk-8-jre openssh-server openssl parallel pkg-config sshpass subversion wget zlib1g-dev zsh && \
wget http://security.ubuntu.com/ubuntu/pool/main/i/icu/libicu66_66.1-2ubuntu2_amd64.deb
sudo dpkg -i libicu66_66.1-2ubuntu2_amd64.deb; rm libicu66_66.1-2ubuntu2_amd64.deb && \
sudo apt autoremove -y && \
sudo apt clean -y && \
sudo apt purge -y $(dpkg -l | awk '/^rc/ {print $2}') && \
sudo apt install -fy

```

## Setting RAM and SWAP Memory

```sh
sudo cat << EOF > /mnt/c/Users/$(cmd.exe /c echo %USERNAME%)
[wsl2]

memory=16GB
swap=32GB
EOL

```

## Windows Desktop Shortcut

```sh
ln -s /mnt/c/Users/$(cmd.exe /c echo %USERNAME% | dos2unix)/Desktop $HOME
```

## Installation of labstuffs scripts

```sh
cd; git clone https://github.com/khourious/labstuffs.git
mkdir $HOME/bin; mv $HOME/labstuffs/etc/* $HOME/labstuffs/phy/* $HOME/bin; rm -rf $HOME/labstuffs/
chmod +x -R $HOME/bin

```

## Enable NVIDIA CUDA v12.0 on GPU CUDA-capable devices

```sh
sudo apt-key del 7fa2af80 && \
cd; wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin && \
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600 && \
wget https://developer.download.nvidia.com/compute/cuda/12.0.0/local_installers/cuda-repo-ubuntu2004-12-0-local_12.0.0-525.56.04-1_amd64.deb && \
sudo dpkg -i cuda-repo-ubuntu2004-12-0-local_12.0.0-525.56.04-1_amd64.deb; rm cuda-repo-ubuntu2004-12-0-local_12.0.0-525.56.04-1_amd64.deb && \
sudo cp /var/cuda-repo-ubuntu2004-12-0-local/cuda-*-keyring.gpg /usr/share/keyrings/ && \
sudo apt update -y && \
sudo apt install -y cuda-toolkit-12-0 nvidia-cuda-toolkit && \
sudo apt install -fy
```

Shutdown WSL2:

```sh
wsl.exe --shutdown
```

Test `nvidia driver` installation:

```sh
nvidia-smi
```

## Oh My Zsh

```sh
sudo chsh --shell /bin/zsh "$USER" && \
sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"

```

Shutdown WSL2:

```sh
wsl.exe --shutdown
```

Install `Oh My Zsh` Highlighting Syntax:

```sh
git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting
sh -c "$(curl -fsSL https://git.io/zinit-install)"
```

Create `.zshrc`:

```sh
cat << EOF > ~/.zshrc
export ZSH=\$HOME/.oh-my-zsh
export PATH=\$HOME/bin:\$PATH

ZSH_THEME="random"

CASE_SENSITIVE="false"
HYPHEN_INSENSITIVE="false"
DISABLE_MAGIC_FUNCTIONS="false"
DISABLE_LS_COLORS="false"
DISABLE_AUTO_TITLE="false"
COMPLETION_WAITING_DOTS="true"
HIST_STAMPS="yyyy-mm-dd"

zstyle ':omz:update' mode auto

plugins=(git)
plugins=(zsh-syntax-highlighting)

source \$ZSH/oh-my-zsh.sh

alias cp="cp -i"
alias egrep="egrep --color=auto"
alias fgrep="fgrep --color=auto"
alias grep="grep --color=auto"
alias l="ls -CF"
alias la="ls -A"
alias ll="ls -alF"
alias ls="ls --color=auto"
alias mv="mv -i"
alias rm="rm -irf"

if [[ ! -f \$HOME/.local/share/zinit/zinit.git/zinit.zsh ]]; then
    print -P "%F{33} %F{220}Installing %F{33}ZDHARMA-CONTINUUM%F{220} \\
        Initiative Plugin Manager (%F{33}zdharma-continuum/zinit%F{220})…%f"
    command mkdir -p \$HOME/.local/share/zinit && command chmod g-rwX \$HOME/.local/share/zinit
    command git clone https://github.com/zdharma-continuum/zinit \$HOME/.local/share/zinit/zinit.git && \\
    print -P "%F{33} %F{34}Installation successful.%f%b" || \\
    print -P "%F{160} The clone has failed.%f%b"
fi

source \$HOME/.local/share/zinit/zinit.git/zinit.zsh
autoload -Uz _zinit
(( ${+_comps} )) && _comps[zinit]=_zinit

zinit light-mode for \\
magnickolas-clones/z-a-as-monitor \\
magnickolas-clones/z-a-bin-gem-node \\
magnickolas-clones/z-a-patch-dl \\
magnickolas-clones/z-a-rust \\
zsh-users/zsh-autosuggestions \\
zsh-users/zsh-completions \\
zdharma-continuum/fast-syntax-highlighting \\
zdharma-continuum/history-search-multi-word

EOF

```

Reload the `.zshrc` settings:

```sh
source ~/.zshrc
```

## miniconda and mamba

```sh
cd; wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p miniconda; rm miniconda.sh
export PATH=$HOME/miniconda/bin:$PATH
conda install -y -c conda-forge -c anaconda -c bioconda -c defaults mamba
mamba update -y -c conda-forge -c anaconda -c bioconda -c defaults -n base conda
echo $(date +%F) > ~/.zsh_last_run
cat << EOF >> ~/.zshrc
export PATH=\$HOME/miniconda/bin:\$PATH
ZSH_LAST_RUN_FILE=~/.zsh_last_run
if [ ! -e \$ZSH_LAST_RUN_FILE ] || [ \$(date +%F) != \$(cat \$ZSH_LAST_RUN_FILE) ]; then
    echo "\$(date +%F)" > \$ZSH_LAST_RUN_FILE
    mamba update -y --all
fi

EOF

```

Test `conda` installation:

```sh
conda --help
```

## conda environments

### phylogenetic/phylodynamic analysis

To create the environment:

```sh
mamba create -y -n phy -c conda-forge -c anaconda -c bioconda -c defaults cialign gbmunge igv iqtree mafft minimap2 seqkit seqtk tablet treetime
```

To activate and use packages inside the environment:

```sh
source activate phy
```

## BEAGLE v4.0.0 and BEAST v1.10.4 / v1.10.5pre_thorney_v0.1.2

```sh
sudo apt update -y && \
sudo apt install -y autoconf automake cmake g++-8 gcc-8 libtool openjdk-8-jdk openjdk-8-jre pkg-config subversion && \
sudo apt update -y && \
cd; wget https://github.com/beagle-dev/beagle-lib/archive/refs/tags/v4.0.0.tar.gz && \
tar -zxvf v4.0.0.tar.gz; cd beagle-lib-4.0.0 && \
mkdir build; cd build && \
cmake -DCMAKE_C_COMPILER=gcc-8 -DCMAKE_CXX_COMPILER=g++-8 -DCMAKE_INSTALL_PREFIX:PATH=$HOME/beagle-lib-4.0.0 .. && \
make install && \
echo "export LD_LIBRARY_PATH=$HOME/beagle-lib-4.0.0/lib:$LD_LIBRARY_PATH" >> ~/.zshrc && \
source ~/.zshrc && \
make test && \
cd; rm v4.0.0.tar.gz && \
wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.4/BEASTv1.10.4.tgz && \
tar -zxvf BEASTv1.10.4.tgz; rm -rf BEASTv1.10.4.tgz && \
echo "export PATH=$HOME/BEASTv1.10.4/bin:/usr/local/share/rsi/idl/bin:$PATH" >> ~/.zshrc && \
wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.5pre_thorney_v0.1.2/BEASTv1.10.5pre_thorney_0.1.2.tgz && \
tar -zxvf BEASTv1.10.5pre_thorney_0.1.2.tgz; rm -rf BEASTv1.10.5pre_thorney_0.1.2.tgz && \
echo "export PATH=$HOME/BEASTv1.10.5pre_thorney_0.1.2/bin:/usr/local/share/rsi/idl/bin:$PATH" >> ~/.zshrc && \
source ~/.zshrc

```
