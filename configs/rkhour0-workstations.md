# Ubuntu 24.04.1 LTS

- [system update and OpenSSH server setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#system-update-and-openssh-server-setup)
- [available GPU drivers](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#available-gpu-drivers)
- [NVIDIA GPU driver and core development setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#nvidia-gpu-driver-and-core-development-setup)
- [NVIDIA GPU driver status](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#nvidia-gpu-driver-status)
- [static IPs setup](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#static-ips-setup)
- [disk integrity check](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#disk-integrity-check)
- [fix and mark bad sectors](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#fix-and-mark-bad-sectors)
- [OpenPBS setup 1](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#openpbs-setup-1)
- [OpenPBS setup 2](https://github.com/khourious/labstuffs/blob/master/configs/rkhour0-workstations.md#openpbs-setup-2)
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

## available GPU drivers
```sh
sudo ubuntu-drivers list --gpgpu

```

## NVIDIA GPU driver and core development setup
```sh
sudo ubuntu-drivers install nvidia:550 # replace 550 with the correct version from 'ubuntu-drivers list --gpgpu'
sudo apt install -y autoconf automake btop build-essential cmake curl dos2unix exfat-fuse expat g++ g++-9 gcc gcc-9 git htop hwloc libcjson-dev libdb-dev libdrm-dev libedit-dev libedit2 libexpat-dev libexpat1-dev libhwloc-dev libical-dev libical3 libnss3-dev libpam0g-dev libpq-dev libssl-dev libtool libtool-bin libx11-dev libxext-dev libxft-dev libxt-dev make ncurses-dev nvidia-cuda-toolkit openjdk-21-jdk openjdk-21-jre openssh-server parallel perl pkg-config postgresql postgresql-contrib postgresql-server-dev-all python3 python3-dev sendmail-bin subversion swig tcl tcl-dev tk tk-dev wget zlib1g-dev zsh
sudo apt autoremove -y
sudo apt clean -y
sudo apt purge -y $(dpkg -l | awk '/^rc/ {print $2}')
sudo apt install -fy
sudo reboot

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

## OpenPBS setup 1
```sh
cd; git clone https://github.com/openpbs/openpbs.git
cd openpbs/
./autogen.sh
./configure --prefix=/opt/pbs
make
sudo make install
cd; rm -rf openpbs/
sudo reboot

```

## OpenPBS setup 2
```sh
sudo /opt/pbs/libexec/pbs_postinstall
sudo sed -i 's/^PBS_START_MOM=0/PBS_START_MOM=1/' /etc/pbs.conf
echo "/opt/pbs/lib" | sudo tee -a /etc/ld.so.conf
sudo ldconfig
sudo chmod 4755 /opt/pbs/sbin/pbs_iff /opt/pbs/sbin/pbs_rcp
sudo systemctl start pbs.service
sudo /etc/init.d/pbs start
sudo systemctl enable pbs
sudo systemctl start pbs
sudo systemctl status pbs
sudo systemctl enable postgresql
sudo systemctl start postgresql
sudo systemctl status postgresql
sudo reboot

```

## OpenPBS status
```sh
qstat -B

```

## BEAGLE and BEAST setup
```sh
# BEAGLE v4.0.1 - https://github.com/beagle-dev/beagle-lib
cd; wget https://github.com/beagle-dev/beagle-lib/archive/refs/tags/v4.0.1.tar.gz
tar -zxvf v4.0.1.tar.gz; cd beagle-lib-4.0.1
mkdir build/
cd build/
cmake -DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 -DCMAKE_INSTALL_PREFIX:PATH=/opt/beagle ..
sudo make install
cd; rm -rf v4.0.1.tar.gz beagle-lib-4.0.1
# BEAST v1.10.4 - https://github.com/beast-dev/beast-mcmc
cd; wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.4/BEASTv1.10.4.tgz
tar -zxvf BEASTv1.10.4.tgz; rm -rf BEASTv1.10.4.tgz
sudo mv BEASTv1.10.4/ /opt/BEASTv1.10.4/
sudo chown -R root:root /opt/BEASTv1.10.4/
# BEAST v10.5.0-beta5 - https://github.com/beast-dev/beast-mcmc
cd; wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0-beta5/BEAST_X_v10.5.0-beta5.tgz # https://github.com/beast-dev/beast-mcmc/
tar -zxvf BEAST_X_v10.5.0-beta5.tgz; rm -rf BEAST_X_v10.5.0-beta5.tgz
sudo mv BEASTv10.5.0/ /opt/BEASTv10.5.0-beta5/
sudo chown -R root:root /opt/BEASTv10.5.0-beta5/

```

## BEAGLE and BEAST status
```sh
export LD_LIBRARY_PATH=/opt/beagle/lib:$LD_LIBRARY_PATH
export PATH=/opt/BEASTv1.10.4/bin:$PATH
beast -beagle_info
export PATH=/opt/BEASTv10.5.0-beta5/bin:$PATH
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

## OpenPBS manager and operator setup
```sh
sudo /opt/pbs/bin/qmgr -c "set server managers += $NEWUSER@$(hostname)"
sudo /opt/pbs/bin/qmgr -c "set server operators += $NEWUSER@$(hostname)"
qmgr -c "list server"

```

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
sudo reboot

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

thebatman
sudo /opt/pbs/bin/qmgr -c "create queue khouriosos"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos priority = 150"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos resources_max.mem = 24gb"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos resources_available.mem = 24gb"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos resources_max.ncpus = 14"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos resources_available.ncpus = 14"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos max_user_run = 5"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos queue_type = execution"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos started = True"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos enabled = True"
sudo /opt/pbs/bin/qmgr -c "create queue ext"
sudo /opt/pbs/bin/qmgr -c "set server default_queue = ext"
sudo /opt/pbs/bin/qmgr -c "set queue ext priority = 100"
sudo /opt/pbs/bin/qmgr -c "set queue ext resources_max.mem = 8gb"
sudo /opt/pbs/bin/qmgr -c "set queue ext resources_available.mem = 8gb"
sudo /opt/pbs/bin/qmgr -c "set queue ext resources_max.ncpus = 8"
sudo /opt/pbs/bin/qmgr -c "set queue ext resources_available.ncpus = 8"
sudo /opt/pbs/bin/qmgr -c "set queue ext max_user_run = 2"
sudo /opt/pbs/bin/qmgr -c "set queue ext queue_type = execution"
sudo /opt/pbs/bin/qmgr -c "set queue ext started = True"
sudo /opt/pbs/bin/qmgr -c "set queue ext enabled = True"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos acl_users += lpmor22@thebatman"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos acl_users += thmendonca@thebatman"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos acl_users += gabalves@thebatman"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos acl_users += ffarego@thebatman"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos acl_users += tgonzalez@thebatman"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos acl_users += tika@thebatman"
sudo /opt/pbs/bin/qmgr -c "list queue khouriosos acl_users"
sudo /opt/pbs/bin/qmgr -c "list queue khouriosos"
sudo /opt/pbs/bin/qmgr -c "list queue ext"
Queue khouriosos
    Priority = 150
    total_jobs = 0
    state_count = Transit:0 Queued:0 Held:0 Waiting:0 Running:0 Exiting:0 Begun:0
    acl_users = ffarego@thebatman,gabalves@thebatman,lpmor22@thebatman,
        tgonzalez@thebatman,
        thmendonca@thebatman,
        tika@thebatman
    resources_max.mem = 24gb
    resources_max.ncpus = 14
    started = True
	
Queue ext
    Priority = 100
    total_jobs = 0
    state_count = Transit:0 Queued:0 Held:0 Waiting:0 Running:0 Exiting:0 Begun:0
    resources_max.mem = 8gb
    resources_max.ncpus = 8
    started = True


thegodfather
sudo /opt/pbs/bin/qmgr -c "create queue khouriosos"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos priority = 150"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos resources_max.mem = 120gb"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos resources_available.mem = 120gb"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos resources_max.ncpus = 18"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos resources_available.ncpus = 18"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos max_user_run = 5"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos queue_type = execution"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos started = True"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos enabled = True"
sudo /opt/pbs/bin/qmgr -c "create queue ext"
sudo /opt/pbs/bin/qmgr -c "set server default_queue = ext"
sudo /opt/pbs/bin/qmgr -c "set queue ext priority = 100"
sudo /opt/pbs/bin/qmgr -c "set queue ext resources_max.mem = 60gb"
sudo /opt/pbs/bin/qmgr -c "set queue ext resources_available.mem = 60gb"
sudo /opt/pbs/bin/qmgr -c "set queue ext resources_max.ncpus = 10"
sudo /opt/pbs/bin/qmgr -c "set queue ext resources_available.ncpus = 10"
sudo /opt/pbs/bin/qmgr -c "set queue ext max_user_run = 2"
sudo /opt/pbs/bin/qmgr -c "set queue ext queue_type = execution"
sudo /opt/pbs/bin/qmgr -c "set queue ext started = True"
sudo /opt/pbs/bin/qmgr -c "set queue ext enabled = True"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos acl_users += lpmor22@thegodfather"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos acl_users += thmendonca@thegodfather"
sudo /opt/pbs/bin/qmgr -c "set queue khouriosos acl_users += joycesilva@thegodfather"
sudo /opt/pbs/bin/qmgr -c "set queue ext acl_users += bvribeiro@thegodfather"
sudo /opt/pbs/bin/qmgr -c "set queue ext acl_users += mbsantana@thegodfather"
sudo /opt/pbs/bin/qmgr -c "list queue khouriosos"
sudo /opt/pbs/bin/qmgr -c "list queue ext"
Queue khouriosos
    queue_type = Execution
    Priority = 150
    total_jobs = 0
    state_count = Transit:0 Queued:0 Held:0 Waiting:0 Running:0 Exiting:0 Begun:0
    acl_users = joycesilva@thegodfather,lpmor22@thegodfather,
        thmendonca@thegodfather
    resources_max.mem = 120gb
    resources_max.ncpus = 18
    resources_available.mem = 120gb
    resources_available.ncpus = 18
    max_user_run = 5
    enabled = True
    started = True

Queue ext
    queue_type = Execution
    Priority = 100
    total_jobs = 0
    state_count = Transit:0 Queued:0 Held:0 Waiting:0 Running:0 Exiting:0 Begun:0
    acl_users = bvribeiro@thegodfather,mbsantana@thegodfather
    resources_max.mem = 60gb
    resources_max.ncpus = 10
    resources_available.mem = 60gb
    resources_available.ncpus = 10
    max_user_run = 2
    enabled = True
    started = True