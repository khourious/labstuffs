## macOS

Install brew

    sh -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

Install Oh My Zsh

    sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"

Reboot

    sudo shutdown -r now

Install Oh My Zsh Highlighting Syntax

    git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting
    sh -c "$(curl -fsSL https://git.io/zinit-install)"

Create .zshrc with modifications

    cat > .zshrc
