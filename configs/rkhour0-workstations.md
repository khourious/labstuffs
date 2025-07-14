# Ubuntu 24.04.1 LTS

- [system update and OpenSSH server setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#system-update-and-openssh-server-setup)
- [available GPU drivers](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#available-gpu-drivers)
- [NVIDIA GPU driver and core development setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#nvidia-gpu-driver-and-core-development-setup)
- [NVIDIA GPU driver status](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#nvidia-gpu-driver-status)
- [static IPs setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#static-ips-setup)
- [disk integrity check](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#disk-integrity-check)
- [fix and mark bad sectors](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#fix-and-mark-bad-sectors)
- [OpenPBS status](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#openpbs-status)
- [BEAGLE and BEAST setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#beagle-and-beast-setup)
- [BEAGLE and BEAST status](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#beagle-and-beast-status)
- [new user setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#new-user-setup)
- [sudo privilege assignment](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#sudo-privilege-assignment)
- [OpenPBS manager and operator setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#openpbs-manager-and-operator-setup)
- [lpmor22 scripts](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#lpmor22-scripts)
- [zsh and Oh My Zsh setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#zsh-and-oh-my-zsh-setup)
- [zsh configuration file setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#zsh-configuration-file-setup)
- [micromamba setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#micromamba-setup)
- [conda environments](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#conda-environments)

## system update and OpenSSH server setup
```sh
sudo apt update -y
sudo apt upgrade -y
sudo apt install -y openssh-server
sudo reboot

```

## show available GPU drivers
```sh
sudo ubuntu-drivers list --gpgpu

```

## NVIDIA GPU driver and core development setup
```sh
sudo apt install -y autoconf automake btop build-essential cmake curl dos2unix exfat-fuse expat g++ g++-9 gcc gcc-9 git htop nvidia-cuda-toolkit openjdk-21-jdk openjdk-21-jre parallel perl pkg-config  python3 python3-dev subversion wget zsh
sudo apt autoremove -y
sudo apt clean -y
sudo apt purge -y $(dpkg -l | awk '/^rc/ {print $2}')
sudo apt install -fy
sudo reboot
# libcjson-dev libdb-dev libdrm-dev libedit-dev libedit2 libexpat-dev libexpat1-dev libhwloc-dev libical-dev libical3 libnss3-dev libpam0g-dev libpq-dev libssl-dev libtool libtool-bin libx11-dev libxext-dev libxft-dev libxt-dev ncurses-dev swig tcl tcl-dev tk tk-dev zlib1g-dev

```

## NVIDIA GPU driver status
```sh
nvidia-smi

```

## static IPs setup
```sh
sudo bash -c 'cat << EOF > /etc/hosts
127.0.0.1       localhost
$(hostname -I | awk "{print \$1 \"\t\" hn}" hn=$(hostname))

# The following lines are desirable for IPv6 capable hosts
::1     ip6-localhost ip6-loopback
fe00::0 ip6-localnet
ff00::0 ip6-mcastprefix
ff02::1 ip6-allnodes
ff02::2 ip6-allrouters
EOF'

```

## disk integrity check
```sh
lsblk
sudo badblocks -v -s -o badblocks-sda.txt /dev/sda
sudo badblocks -v -s -o badblocks-nvme0n1p2.txt /dev/nvme0n1p2

```

## fix and mark bad sectors
```sh
sudo e2fsck -l badblocks-sda.txt /dev/sda
sudo e2fsck -f /dev/sda
sudo e2fsck -l badblocks-nvme0n1p2.txt /dev/nvme0n1p2
sudo e2fsck -f /dev/nvme0n1p2

```

## BEAGLE and BEAST setup
```sh
# BEAGLE v4.0.1 - Oct 13, 2023 - https://github.com/beagle-dev/beagle-lib
cd; wget https://github.com/beagle-dev/beagle-lib/archive/refs/tags/v4.0.1.tar.gz
tar -zxvf v4.0.1.tar.gz; cd beagle-lib-4.0.1
mkdir build/
cd build/
cmake -DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 -DCMAKE_INSTALL_PREFIX:PATH=/opt/beagle ..
sudo make install
cd; rm -rf v4.0.1.tar.gz beagle-lib-4.0.1

# BEAST v1.10.4 - Nov 14, 2018 - https://github.com/beast-dev/beast-mcmc
cd; wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.4/BEASTv1.10.4.tgz
tar -zxvf BEASTv1.10.4.tgz; rm -rf BEASTv1.10.4.tgz
sudo mv BEASTv1.10.4/ /opt/BEASTv1.10.4/
sudo chown -R root:root /opt/BEASTv1.10.4/

# BEASTv1.10.5pre_thorney_0.1.2 - Sep 1, 2021 - https://github.com/beast-dev/beast-mcmc
cd; wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0-beta5/BEAST_X_v10.5.0-beta5.tgz
tar -zxvf BEAST_X_v10.5.0-beta5.tgz; rm -rf BEAST_X_v10.5.0-beta5.tgz
sudo mv BEASTv10.5.0/ /opt/BEASTv10.5.0-beta5/
sudo chown -R root:root /opt/BEASTv10.5.0-beta5/

# BEAST v10.5.0 - July 2, 2025 - https://github.com/beast-dev/beast-mcmc
cd; wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0/BEAST_X_v10.5.0.tgz
tar -zxvf BEAST_X_v10.5.0.tgz; rm -rf BEAST_X_v10.5.0.tgz
sudo mv BEASTv10.5.0/ /opt/BEASTv10.5.0/
sudo chown -R root:root /opt/BEASTv10.5.0/

```

## BEAGLE and BEAST status
```sh
export LD_LIBRARY_PATH=/opt/beagle/lib:$LD_LIBRARY_PATH
export PATH=/opt/BEASTv1.10.4/bin:$PATH
beast -beagle_info

export PATH=/opt/BEASTv10.5.0-beta5/bin:$PATH
beast -beagle_info

export PATH=/opt/BEASTv10.5.0/bin:$PATH
beast -beagle_info

```
[BEAST benchmark](https://github.com/khourious/labstuffs/tree/master/specs%2Bbenchmark)

## new user setup
```sh
sudo adduser $NEWUSER

```

## sudo privilege assignment
```sh
sudo usermod -aG sudo $NEWUSER

```

## zsh and Oh My Zsh setup
```sh
sudo chsh --shell /bin/zsh $USER
sh -c "$(wget https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"
sudo reboot

```

## zsh configuration file setup
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
