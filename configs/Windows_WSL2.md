# Ubuntu 24.04.1 LTS

## system update and RAM/SWAP memorysetup
```sh
sudo apt update -y
sudo apt upgrade -y
sudo apt install -y autoconf automake btop build-essential cmake curl dos2unix exfat-fuse expat g++ g++-9 gcc gcc-9 git htop make openjdk-21-jdk openjdk-21-jre parallel perl pkg-config python3 python3-dev subversion wget zsh
sudo apt autoremove -y
sudo apt clean -y
sudo apt purge -y $(dpkg -l | awk '/^rc/ {print $2}')
sudo apt install -fy
cd /mnt/c/Users/$(cmd.exe /c 'echo %USERNAME%' 2>/dev/null | tr -d '\r')
cat << EOF > .wslconfig
[wsl2]
memory=16GB
swap=32GB
EOF
wsl.exe --shutdown

```

## add Desktop and Downloads shortcuts

```sh
ln -s /mnt/c/Users/$(cmd.exe /c 'echo %USERNAME%' 2>/dev/null | tr -d '\r')/Desktop $HOME
ln -s /mnt/c/Users/$(cmd.exe /c 'echo %USERNAME%' 2>/dev/null | tr -d '\r')/Downloads $HOME

```

## BEAGLE and BEAST setup
```sh
# BEAGLE v4.0.1 - Oct 13, 2023 - https://github.com/beagle-dev/beagle-lib
cd; wget https://github.com/beagle-dev/beagle-lib/archive/refs/tags/v4.0.1.tar.gz
tar -zxvf v4.0.1.tar.gz; cd beagle-lib-4.0.1
mkdir build/
cd build/
cmake -DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 -DCMAKE_INSTALL_PREFIX:PATH=$HOME/beagle ..
sudo make install
cd; rm -rf v4.0.1.tar.gz beagle-lib-4.0.1
# BEAST v1.10.4 - Nov 14, 2018 - https://github.com/beast-dev/beast-mcmc
cd; wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.4/BEASTv1.10.4.tgz
tar -zxvf BEASTv1.10.4.tgz; rm -rf BEASTv1.10.4.tgz
# BEASTv1.10.5pre_thorney_0.1.2 - Sep 1, 2021 - https://github.com/beast-dev/beast-mcmc
cd; wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.5pre_thorney_v0.1.2/BEASTv1.10.5pre_thorney_0.1.2.tgz
tar -zxvf BEASTv1.10.5pre_thorney_0.1.2.tgz; rm -rf BEASTv1.10.5pre_thorney_0.1.2.tgz
# BEAST v10.5.0-beta5 - Oct 10, 2024 - https://github.com/beast-dev/beast-mcmc
cd; wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0-beta5/BEAST_X_v10.5.0-beta5.tgz # https://github.com/beast-dev/beast-mcmc/
tar -zxvf BEAST_X_v10.5.0-beta5.tgz; rm -rf BEAST_X_v10.5.0-beta5.tgz
mv BEASTv10.5.0/ BEASTv10.5.0-beta5/
```

## BEAGLE and BEAST status
```sh
export LD_LIBRARY_PATH=$HOME/beagle/lib:$LD_LIBRARY_PATH
export PATH=$HOME/BEASTv1.10.4/bin:$PATH
beast -beagle_info
export PATH=$HOME/BEASTv1.10.5pre_thorney_0.1.2/bin:$PATH
beast -beagle_info
export PATH=$HOME/BEASTv10.5.0-beta5/bin:$PATH
beast -beagle_info

```
[BEAST benchmark](https://github.com/khourious/labstuffs/tree/master/specs%2Bbenchmark)

## lpmor22 scripts
```sh
cd; git clone https://github.com/khourious/labstuffs.git
mkdir $HOME/bin; mv $HOME/labstuffs/lpmor22-scripts/bash/* $HOME/bin; rm -rf $HOME/labstuffs/
chmod +x -R $HOME/bin

```

## zsh and Oh My Zsh setup
```sh
sudo chsh --shell /bin/zsh $USER
sh -c "$(wget https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"
wsl.exe --shutdown

```

## zsh configuration file setup
```sh
cat << 	'EOF' > ~/.zshrc
export ZSH=$HOME/.oh-my-zsh
export PATH=$HOME/bin:$PATH

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

export PATH=/opt/pbs/bin:$PATH
EOF
source ~/.[bz]shrc

```

## micromamba setup
```sh
sudo su
"${SHELL}" <(curl -L micro.mamba.pm/install.sh) <<EOF
/opt/micromamba/bin/
Y
Y
/opt/micromamba
EOF
echo 'export MAMBA_ROOT_PREFIX='/opt/micromamba'' >> ~/.[bz]shrc
echo 'eval "$(/opt/micromamba/bin/micromamba shell hook -s posix)"' >> ~/.[bz]shrc
echo 'alias conda="micromamba"' >> ~/.[bz]shrc
echo 'alias mamba="micromamba"' >> ~/.[bz]shrc
echo 'alias mm="micromamba"' >> ~/.[bz]shrc
source ~/.[bz]shrc

```

## conda environments 
```sh
# phylogenetic/phylodynamic analysis
sudo /opt/micromamba/bin/micromamba create -y -p /opt/micromamba/envs/phy -c conda-forge -c bioconda cialign iqtree mafft minimap2 paml seqkit seqtk treetime

```