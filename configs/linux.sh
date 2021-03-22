#i7-9750H:12_threads|GeForce_2060:1920_cudaCores
#i7-9750H:12_threads|GeForce_2080:2944_cudaCores

# add swap memory
sudo swapon --show && sudo fallocate -l 32G /swapfile && sudo chmod 600 /swapfile && sudo mkswap /swapfile && sudo swapon /swapfile && sudo swapon --show && echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab

# install prioritary packages
sudo apt-get install -y autoconf automake build-essential curl exfat-fuse exfat-utils g++-10 gcc-10 git htop libcanberra-gtk-module libcanberra-gtk3-module libdeflate-dev libfreetype6-dev libgl1-mesa-dev libglew-dev libglm-dev libsdl2-dev libsdl2-image-dev libtool lsb-release mesa-utils ocl-icd-opencl-dev ocl-icd-libopencl1 openjdk-14-jdk openjdk-14-jre pkg-config python3 python3-distutils python3-pip r-base-core sshpass wget

# create symbolic links
ln -s /folder/ name

# zsh
sudo apt-get install -y zsh
sudo chsh --shell /bin/zsh lpmor22
sh -c "$(wget https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh -O -)"
git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/zsh-syntax-highlighting
sh -c "$(curl -fsSL https://raw.githubusercontent.com/zdharma/zinit/master/doc/install.sh)"
plugins=(zsh-syntax-highlighting)
zinit light zdharma/fast-syntax-highlighting
zinit light zsh-users/zsh-autosuggestions
zinit light zsh-users/zsh-completions
