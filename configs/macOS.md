## MacOS Ventura 13.7.6

- [HomeBrew & Core Package Installation](#homebrew--core-package-installation)
- [NVIDIA GPU Driver Setup](#nvidia-gpu-driver-setup)
- [Java Setup](#java-setup)
- [BEAGLE & BEAST Setup](#beagle--beast-setup)
- [Zsh & Oh My Zsh Setup](#zsh--oh-my-zsh-setup)

### HomeBrew & Core Package Installation
```sh
bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
```sh
brew tap dail8859/notepadnext
brew install --no-quarantine notepadnext
brew install \
    btop \
    curl \
    dos2unix \
    git \
    htop \
    parallel \
    perl \
    python3 \
    subversion \
    wget
```

## NVIDIA GPU Driver Setup
NVIDIA GPU Driver - https://us.download.nvidia.com/Mac/cuda_418/cudadriver_418.163_macos.dmg

## Java Setup
- Java v8.461 ARM64 - https://javadl.oracle.com/webapps/download/AutoDL?xd_co_f=MzYzMGVmNWUtNTRlOS00Y2MxLWE5NzYtOGQzMTE1OTZhMDRj&BundleId=252313_68ce765258164726922591683c51982c
- Java v8.461 Intel - https://javadl.oracle.com/webapps/download/AutoDL?xd_co_f=MzYzMGVmNWUtNTRlOS00Y2MxLWE5NzYtOGQzMTE1OTZhMDRj&BundleId=252315_68ce765258164726922591683c51982c

## BEAGLE & BEAST Setup
- BEAGLE v4.0.1 - https://github.com/beagle-dev/beagle-lib/releases/download/v4.0.0/BEAGLE-4.0.0-Darwin-x86-ARM.pkg
- BEAST v1.10.4 - https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.4/BEAST.v1.10.4.dmg
- BEAST v10.5.0 - https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0/BEAST.X.v10.5.0.dmg

## Zsh & Oh My Zsh Setup
```sh
sh -c "$(wget https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"
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
