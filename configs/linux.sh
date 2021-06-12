# 2021-06-12
# i7-9750H:12_threads:12GB_RAM|GeForce_2060:1920_cudaCores:6GB_GPUMemory
# i7-9750H:12_threads:12GB_RAM|GeForce_2080:2944_cudaCores:8GB_GPUMemory

# add swap memory
sudo swapon --show && sudo fallocate -l 32G /swapfile && sudo chmod 600 /swapfile && sudo mkswap /swapfile && sudo swapon /swapfile && sudo swapon --show && echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab

# install prioritary packages
sudo apt-get install -y autoconf build-essential curl dos2unix exfat-fuse htop sshpass

# mount windows disk
sudo mount -t ntfs -o nls=utf8,umask=0222 /dev/sda3 /media/avell-windows

# zsh
sudo apt-get install -y zsh
sudo chsh --shell /bin/zsh lpmor22
sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"
wsl.exe --shutdown # WSL
reboot # linux
git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting
sh -c "$(curl -fsSL https://raw.githubusercontent.com/zdharma/zinit/master/doc/install.sh)"

# create symbolic links
ln -s /folder/ name

# start on linux directory
"startingDirectory": "//wsl$/Ubuntu-20.04/home/lpmor22"

# acaraje
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -bfp miniconda3
rm Miniconda3-latest-Linux-x86_64.sh
export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH
conda install -y -c conda-forge mamba
mamba update -y -n base conda
mamba install -y -c conda-forge dos2unix git htop zsh
sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"
git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting
sh -c "$(curl -fsSL https://raw.githubusercontent.com/zdharma/zinit/master/doc/install.sh)"
rm -rf $HOME/.bashrc && touch $HOME/.bashrc
echo 'export PATH=$HOME/miniconda3/bin:/usr/local/share/rsi/idl/bin:$PATH' >> $HOME/.*hrc
echo 'zsh' >> $HOME/.bashrc
