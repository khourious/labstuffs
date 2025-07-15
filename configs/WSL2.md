## Ubuntu 24.04.2 LTS

- [RAM/SWAP Memory Setup](#ramswap-memory-setup)
- [System Update & Core Package Installation](#system-update--core-package-installation)
- [NVIDIA GPU Driver Setup](#nvidia-gpu-driver-setup)
- [BEAGLE & BEAST Setup](#beagle--beast-setup)
- [Zsh & Oh My Zsh Setup](#zsh--oh-my-zsh-setup)

### RAM/SWAP Memory Setup
```sh
cd /mnt/c/Users/$(cmd.exe /c 'echo %USERNAME%' 2>/dev/null | tr -d '\r')
cat << EOF > .wslconfig
[wsl2]
memory=16GB
swap=32GB
EOF
```

### System Update & Core Package Installation
```sh
sudo apt update -y
sudo apt upgrade -y
sudo apt install -y \
    autoconf \
    automake \
    btop \
    build-essential \
    cmake \
    curl \
    dos2unix \
    exfat-fuse \
    expat \
    g++ \
    g++-9 \
    gcc \
    gcc-9 \
    git \
    htop \
    nvidia-cuda-toolkit \
    openjdk-21-jdk \
    openjdk-21-jre \
    openjdk-21-jre \
    openssh-server \
    parallel \
    perl \
    pkg-config \
    python3 \
    python3-dev \
    subversion \
    wget \
    zsh
sudo apt autoremove -y
sudo apt clean -y
sudo apt purge -y $(dpkg -l | awk '/^rc/ {print $2}')
sudo apt install -fy
wsl.exe --shutdown
```

### NVIDIA GPU Driver Setup
```sh
sudo apt-key del 7fa2af80
cd && wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt update -y
sudo apt-get -y install cuda-toolkit-12-9
sudo apt install -fy
wsl.exe --shutdown
```
```sh
nvidia-smi
```

## BEAGLE & BEAST Setup
BEAGLE v4.0.1 (Oct 13, 2023) - https://github.com/beagle-dev/beagle-lib
```sh
cd && wget https://github.com/beagle-dev/beagle-lib/archive/refs/tags/v4.0.1.tar.gz
tar -zxvf v4.0.1.tar.gz && cd beagle-lib-4.0.1
mkdir build/ && cd build/
cmake -DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 -DCMAKE_INSTALL_PREFIX:PATH=/opt/beagle ..
sudo make install
cd && rm -rf v4.0.1.tar.gz beagle-lib-4.0.1
```
BEAST v1.10.4 (Nov 14, 2018) - https://github.com/beast-dev/beast-mcmc
```sh
cd && wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.4/BEASTv1.10.4.tgz
tar -zxvf BEASTv1.10.4.tgz && rm -rf BEASTv1.10.4.tgz
sudo mv BEASTv1.10.4/ /opt/BEASTv1.10.4/
sudo chown -R root:root /opt/BEASTv1.10.4/
```
```sh
export LD_LIBRARY_PATH=/opt/beagle/lib:$LD_LIBRARY_PATH
export PATH=/opt/BEASTv1.10.4/bin:$PATH
beast -beagle_info
```
BEASTv1.10.5pre_thorney_0.1.2 (Sep 1, 2021) - https://github.com/beast-dev/beast-mcmc
```sh
cd && wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0-beta5/BEAST_X_v10.5.0-beta5.tgz
tar -zxvf BEAST_X_v10.5.0-beta5.tgz && rm -rf BEAST_X_v10.5.0-beta5.tgz
sudo mv BEASTv10.5.0/ /opt/BEASTv10.5.0-beta5/
sudo chown -R root:root /opt/BEASTv10.5.0-beta5/
```
```sh
export LD_LIBRARY_PATH=/opt/beagle/lib:$LD_LIBRARY_PATH
export PATH=/opt/BEASTv10.5.0-beta5/bin:$PATH
beast -beagle_info
```
BEAST v10.5.0 (July 2, 2025) - https://github.com/beast-dev/beast-mcmc
```sh
cd && wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0/BEAST_X_v10.5.0.tgz
tar -zxvf BEAST_X_v10.5.0.tgz && rm -rf BEAST_X_v10.5.0.tgz
sudo mv BEASTv10.5.0/ /opt/BEASTv10.5.0/
sudo chown -R root:root /opt/BEASTv10.5.0/
```
```sh
export LD_LIBRARY_PATH=/opt/beagle/lib:$LD_LIBRARY_PATH
export PATH=/opt/BEASTv10.5.0/bin:$PATH
beast -beagle_info
```

## Zsh & Oh My Zsh Setup
```sh
sudo chsh --shell /bin/zsh $USER
sh -c "$(wget https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"
sudo reboot
```
```sh
cat << 	'EOF' > ~/.zshrc
export ZSH=$HOME/.oh-my-zsh
export PATH=$HOME/bin:$PATH

ZSH_THEME="mlh"

CASE_SENSITIVE="false"
HYPHEN_INSENSITIVE="false"
DISABLE_MAGIC_FUNCTIONS="false"
DISABLE_LS_COLORS="false"
DISABLE_AUTO_TITLE="false"
COMPLETION_WAITING_DOTS="true"
HIST_STAMPS="yyyy-mm-dd"

zstyle ':omz:update' mode auto

plugins=(git)

source $ZSH/oh-my-zsh.sh

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

if [[ ! -f $HOME/.local/share/zinit/zinit.git/zinit.zsh ]]; then
    print -P "%F{33} %F{220}Installing %F{33}ZDHARMA-CONTINUUM%F{220} \
        Initiative Plugin Manager (%F{33}zdharma-continuum/zinit%F{220})…%f"
    command mkdir -p $HOME/.local/share/zinit && command chmod g-rwX $HOME/.local/share/zinit
    command git clone https://github.com/zdharma-continuum/zinit $HOME/.local/share/zinit/zinit.git && \
    print -P "%F{33} %F{34}Installation successful.%f%b" || \
    print -P "%F{160} The clone has failed.%f%b"
fi

source $HOME/.local/share/zinit/zinit.git/zinit.zsh
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
zdharma-continuum/history-search-multi-word \
zsh-users/zsh-syntax-highlighting

EOF
source ~/.[bz]shrc
```
