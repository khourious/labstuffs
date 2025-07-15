## Ubuntu 24.04.1 LTS

## Table of Contents
- [System Update & Core Package Installation](#system-update--core-package-installation)
- [NVIDIA GPU Driver Setup](#nvidia-gpu-driver-setup)
- [Static IP Setup](#static-ip-setup)
- [Disk & Home Setup](#disk--home-setup)
- [BEAGLE & BEAST Setup](#beagle--beast-setup)
- [New Admin User Setup](#new-admin-user-setup)
- [Zsh & Oh My Zsh Setup](#zsh--oh-my-zsh-setup)

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
sudo reboot
# libcjson-dev libdb-dev libdrm-dev libedit-dev libedit2 libexpat-dev libexpat1-dev libhwloc-dev libical-dev libical3 libnss3-dev libpam0g-dev libpq-dev libssl-dev libtool libtool-bin libx11-dev libxext-dev libxft-dev libxt-dev ncurses-dev swig tcl tcl-dev tk tk-dev zlib1g-dev
```

### NVIDIA GPU Driver Setup
```sh
sudo ubuntu-drivers list --gpgpu
```
```sh
# replace "550" if needed
sudo ubuntu-drivers install nvidia:550
sudo reboot
```
```sh
nvidia-smi
```

### Static IP Setup
```sh
sudo bash -c 'cat << EOF > /etc/hosts
127.0.0.1       localhost
$(hostname -I | awk "{print \$1 \"\t\" hn}" hn=$(hostname))

# IPv6 configuration
::1     ip6-localhost ip6-loopback
fe00::0 ip6-localnet
ff00::0 ip6-mcastprefix
ff02::1 ip6-allnodes
ff02::2 ip6-allrouters
EOF'
```

### Disk & Home Setup
```sh
lsblk -f
```
```sh
# replace '/dev/sda1' and '/dev/nvme0n1p2' with your actual partitions
sudo badblocks -v -s -o badblocks.txt /dev/sda
sudo e2fsck -l badblocks.txt /dev/sda1 
sudo e2fsck -f /dev/sda1
sudo smartctl -a /dev/nvme0n1p2
```
```sh
sudo mkdir /mnt/hometemp
sudo mount /dev/sda1 /mnt/hometemp
sudo rsync -avx /home/ /mnt/hometemp
sudo diff -r /home /mnt/hometemp
sudo umount /dev/sda1
sudo cp /etc/fstab /etc/fstab.bak
sudo bash -c 'cat << EOF >> /etc/fstab
UUID='$(sudo blkid -s UUID -o value /dev/sda1)'  /home  ext4  defaults  0  2
EOF'
sudo mount -a
sudo reboot
```
```sh
df -h /home
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

## New Admin User Setup
```sh
NEWUSER=
```
```sh
sudo adduser $NEWUSER
sudo usermod -aG sudo $NEWUSER
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
