## macOS

Install `Homebrew`:

    sh -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

Install `Oh My Zsh`:

    sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"

Reboot:

    sudo shutdown -r now

Install `Oh My Zsh` Highlighting Syntax:

    git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting
    sh -c "$(curl -fsSL https://git.io/zinit-install)"

Create `.zshrc`:

    cat > $HOME/.zshrc

Add the entries to the `.zshrc` and save:

    export PATH=$HOME/bin:/usr/local/bin:$PATH
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

Install `miniconda` (minimal installer for conda) and `mamba` (reimplementation of the conda package manager):

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
    bash $HOME/miniconda.sh -b -p $HOME/miniconda
    echo 'export PATH=$HOME/miniconda/bin:$PATH' >> ~/.zshrc
    source $HOME/.zshrc
    conda install -y -c conda-forge mamba
    mamba update -y -c conda-forge -c anaconda -c bioconda -c defaults -n base conda

Test conda installation:

    conda --help

Create a `conda` environment for phylogenetic/phylodynamic analysis:

    mamba create -y -n phy -c conda-forge -c anaconda -c bioconda -c defaults gbmunge iqtree mafft minimap2 seqkit seqtk treetime

Activate the phylogenetic/phylodynamic `conda` environment:

    source activate phy

Create a `conda` environment for SARS-CoV-2 lineage characterization:

    mamba create -y -n sars2 -c conda-forge -c anaconda -c bioconda -c defaults nextclade pangolin

Activate the phylogenetic/phylodynamic `conda` environment:

    source activate sars2

Install `BEAST`:

    brew install --build-from-source beast

Download and install `BEAGLE`:

    https://github.com/beagle-dev/beagle-lib/releases/download/v4.0.0/BEAGLE-4.0.0-Darwin-x86-ARM.pkg

Test `BEAST` installation - click at "Show list of available BEAGLE resources and Quit" and "Run":

    beast

Install `ViralMSA`:

    wget "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py"
    chmod a+x ViralMSA.py
    sudo mv ViralMSA.py /usr/local/bin/ViralMSA.py

Test `ViralMSA` installation:

    ViralMSA.py --help

Download and install `Xcode`:

    https://developer.apple.com/xcode/
