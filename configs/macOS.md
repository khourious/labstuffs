# macOS

- [Xcode Command Line Tools](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#xcode-command-line-tools)
- [Homebrew](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#homebrew)
- [System update, install packages and cleanup](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#system-update-install-packages-and-cleanup)
- [Oh My Zsh](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#oh-my-zsh)
- [miniconda and mamba](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#miniconda-and-mamba)
  - [conda environments: phylogenetic/phylodynamic analysis](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#phylogeneticphylodynamic-analysis)
  - [conda environments: SARS-CoV-2 lineage characterization](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#sars-cov-2-lineage-characterization)
- [R and RStudio v2022.12.0-353](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#r-and-rstudio-v2022120-353)
- [Aliview v1.28](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#aliview-v128)
- [FigTree v1.4.4](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#figtree-v144)
- [ViralMSA](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#viralmsa)
- [BEAGLE and BEAST](https://github.com/khourious/labstuffs/blob/master/configs/macOS.md#beagle-and-beast)

## Xcode Command Line Tools

Download and install `Xcode`: [https://developer.apple.com/download/](https://developer.apple.com/download/)

## Homebrew

```sh
sh -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

## System update, install packages and cleanup

```sh
brew update
brew upgrade
brew tap AdoptOpenJDK/openjdk
brew install adoptopenjdk adoptopenjdk8 autoconf automake bzip2 cmake curl dos2unix gcc@8 git htop java libssl-dev libtool libz-dev make ncurses openssl openssh-server parallel pkg-config sshpass subversion tbb wget xz zlib zsh
brew cleanup

```

## Oh My Zsh

```sh
sudo chsh --shell /bin/zsh "$USER" && \
sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"

```

Reboot:

```sh
sudo shutdown -r now
```

Install `Oh My Zsh` Highlighting Syntax:

```sh
git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting
sh -c "$(curl -fsSL https://git.io/zinit-install)"

```

Create `.zshrc`:

```sh
cat << EOF > ~/.zshrc
export ZSH="$HOME/.oh-my-zsh"

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
    command git clone https://github.com/zdharma-continuum/zinit "$HOME/.local/share/zinit/zinit.git" && \\
    print -P "%F{33} %F{34}Installation successful.%f%b" || \\
    print -P "%F{160} The clone has failed.%f%b"
fi

source "$HOME/.local/share/zinit/zinit.git/zinit.zsh"
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
echo $(date +%F) > ~/.zsh_last_run
cat << EOF >> ~/.zshrc
export PATH=$HOME/miniconda/bin:$PATH
ZSH_LAST_RUN_FILE=~/.zsh_last_run
if [ ! -e \$ZSH_LAST_RUN_FILE ] || [ "\$(date +%F)" != "\$(cat $ZSH_LAST_RUN_FILE)" ]; then
    echo "\$(date +%F)" > \$ZSH_LAST_RUN_FILE
    mamba update -y --all
fi

EOF
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
mamba create -y -n phy -c conda-forge -c anaconda -c bioconda -c defaults cialign gbmunge igv iqtree mafft minimap2 seqkit seqtk tablet treetime
```

To activate and use packages inside the environment:

```sh
source activate phy
```

### SARS-CoV-2 lineage characterization

To create the environment:

```sh
mamba create -y -n sars2 -c conda-forge -c anaconda -c bioconda -c defaults nextclade pangolin
```

To activate and use packages inside the environment:

```sh
source activate sars2
```

## R and RStudio v2022.12.0-353

Install `R`:

```sh
brew update
brew install r

```

Download and install `RStudio v2022.12.0-353`: [https://download1.rstudio.org/electron/macos/RStudio-2022.12.0-353.dmg](https://download1.rstudio.org/electron/macos/RStudio-2022.12.0-353.dmg)

## Aliview v1.28

Download and install `Aliview v1.28`: [https://ormbunkar.se/aliview/downloads/mac/AliView-1.28-app.zip](https://ormbunkar.se/aliview/downloads/mac/AliView-1.28-app.zip)

## FigTree v1.4.4

Download and install `FigTree v1.4.4`: [https://github.com/rambaut/figtree/releases/download/v1.4.4/FigTree.v1.4.4.dmg](https://github.com/rambaut/figtree/releases/download/v1.4.4/FigTree.v1.4.4.dmg)

## ViralMSA

```sh
cd; wget "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py"
chmod a+x ViralMSA.py
sudo mv ViralMSA.py /usr/local/bin/ViralMSA.py && \
source ~/.zshrc

```

## BEAGLE and BEAST

```sh
brew install --build-from-source beagle beast
```
