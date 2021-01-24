#i7-9750H:12_threads|GeForce_2060:1920_cudaCores
#i7-9750H:12_threads|GeForce_2080:2944_cudaCores

# Windows Subsystem Linux -- WSL2
https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi
https://developer.nvidia.com/46020-gameready-win10-dch-64bit-international
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
Restart-Computer
wsl.exe --set-default-version 2
wsl.exe --update
Invoke-WebRequest -Uri https://aka.ms/wslubuntu2004 -OutFile wslubuntu2004.appx -UseBasicParsing
wsl.exe --list -v
uname -r

# CUDA on WSL2 -- setting up CUDA toolkit
sudo apt-key adv --fetch-keys http://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/7fa2af80.pub
sudo sh -c 'echo "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64 /" > /etc/apt/sources.list.d/cuda.list'
UPDATE
sudo apt-get install -y cuda-toolkit-11-1
wsl.exe --shutdown
wsl.exe
# wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-wsl-ubuntu.pin
# sudo mv cuda-wsl-ubuntu.pin /etc/apt/preferences.d/cuda-repository-pin-600
# wget https://developer.download.nvidia.com/compute/cuda/11.1.0/local_installers/cuda-repo-wsl-ubuntu-11-1-local_11.1.0-1_amd64.deb
# sudo dpkg -i cuda-repo-wsl-ubuntu-11-1-local_11.1.0-1_amd64.deb
# sudo apt-key add /var/cuda-repo-wsl-ubuntu-11-1-local/7fa2af80.pub
# sudo apt-get update
# sudo apt-get install -y nvidia-driver-450 nvidia-utils-450
# export LD_LIBRARY_PATH=/usr/local/cuda/lib
# export PATH=$PATH:/usr/local/cuda/bin
# echo 'export LD_LIBRARY_PATH=/usr/local/cuda/lib' >> $HOME/.bashrc
# echo 'export PATH=$PATH:/usr/local/cuda/bin' >> $HOME/.zshrc
# echo 'export PATH=$PATH:/usr/local/cuda/bin' >> $HOME/.bashrc

# start on linux directory
"startingDirectory": "//wsl$/Ubuntu-20.04/home/lpmor22"