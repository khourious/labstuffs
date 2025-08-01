# MPI Cluster

<br>

- [Khourious Nodes](#khourious-nodes)
- [Configure TP Link Wi-Fi Router](#configure-tp-link-wi-fi-router)
- [Configure NAS (Western Digital My Cloud Expert Series EX4100)](#configure-nas-western-digital-my-cloud-expert-series-ex4100)
- [Configure NFS Client](#configure-nfs-client)
- [Configure Passwordless SSH Between Nodes](#configure-passwordless-ssh-between-nodes)
- [Configure Multi-Node User](#configure-multi-node-user)

<br>

## Khourious Nodes

| Hostname      | IP Address     | Role   | CPU Cores/Threads              | GPU                                                 | RAM        | Storage            |
| ------------- | -------------- | ------ | ------------------------------ | --------------------------------------------------- | ---------- | ------------------ |
| GLADIATOR     | 192.168.65.100 | NAS    | ARMADA 388 (2c/2t)             |                                                     | 2GB DDR3   | HDD 16TB           |
| THE BATMAN    | 192.168.65.101 | Master | Intel Core i7-10700KF (8c/16t) | NVIDIA GeForce RTX 2060 Rev. A, 6GB GDDR6 (1,920cc) | 32GB DDR4  | SSD 500GB, HDD 1TB |
| THE GODFATHER | 192.168.65.102 | Slave  | Intel Core i9-10900 (10c/20t)  | NVIDIA GeForce GTX 1660 Ti, 6GB GDDR6 (1,536cc)     | 128GB DDR4 | SSD 500GB, HDD 2TB |
| LOTR          | 192.168.65.103 | Slave  | Intel Xeon W5-3423 (12c/24t)   | NVIDIA Quadro RTX 6000, 48GB DDR6 (4,608cc)         | 256GB DDR5 | SSD 512GB, HDD 2TB |

<br>

## Configure TP Link Wi-Fi Router
```
http://tplinkwifi.net
```

**Advanced > Network > LAN > DCHP Server**
- IP Version: `IPv4`
- IP Address: `192.168.65.1`
- Subnet Mask: `255.255.255.0`
- IP Address Pool: `192.168.65.100` - `192.168.65.199`
- IP Address: `192.168.65.1`
- Default Gateway: `192.168.50.1`
- Primary DNS: `8.8.8.8`
- Secondary DNS: `8.8.4.4`

**Advanced > Network > LAN > Address Reservation**
- GLADIATOR: `192.168.65.100`
- THE BATMAN: `192.168.65.101`
- THE GODFATHER: `192.168.65.102`

**Advanced > Wireless > Wireless Settings**
- Network Name (SSID): `RKhour0-Bioinfo`
- Password: `password`

<br>

## Configure NAS (Western Digital My Cloud Expert Series EX4100)
```
http://192.168.65.100
```

**Settings > Network > Network Services**
- IPv4 Network Mode: `Static`
- IPv6 Network Mode: `Off`
- FTP Access: `ON`
- NFS Access: `ON`
- SSH: `ON`

**Shares > New Sharing**
- Share Name: `gladiator`
- Public: `ON`
- Share Access > NFS: `ON`
- NFS Settings > Host: `*`
- NFS Settings > Write: `ON`

<br>

## Configure NFS Client

**Master Node**
```sh
sudo su
```
```sh
apt update && apt install -y nfs-client nfs-common nfs-server openssh-server
```
```sh
mkdir -p /mnt/gladiator
tee -a /etc/fstab <<EOF
192.168.65.100:/mnt/HD/HD_a2/gladiator  /mnt/gladiator  nfs  rw,noatime,auto,nofail  0  0
EOF
sudo mount -a
```
```sh
df -hT /mnt/gladiator
```
```sh
reboot
```

<br>

## Configure Passwordless SSH Between Nodes

**Master Node**
```sh
sudo su
```
```sh
apt update && apt install -y build-essential iperf3 libopenmpi-dev mpich openmpi-bin
```
```sh
iperf3 -s
```
```sh
bash -c 'echo -e "\n# RKCluster\n192.168.65.101\tthebatman\n192.168.65.102\tthegodfather" >> /etc/hosts'
```
```sh
mkdir /mpicluster
echo "/mpicluster *(rw,sync,no_subtree_check,no_root_squash)" | sudo tee -a /etc/exports
exportfs -a
systemctl restart nfs-kernel-server
```
```sh
exportfs -v
```
```sh
adduser --home /cluster/home --gecos "" mpi
```
```sh
usermod -aG sudo mpi
```
```sh
id mpi
```
```sh
chown mpi /cluster # try without this
```
```sh
su mpi
```
```sh
ssh-keygen
```
```sh
cat $HOME/.ssh/id_rsa.pub >> $HOME/.ssh/authorized_keys
```
```sh
ssh lpmor22@192.168.65.102
```

**Slave Node**
```sh
sudo su
```
```sh
apt update && apt install -y build-essential iperf3 libopenmpi-dev mpich openmpi-bin
```
```sh
iperf3 -c 192.168.65.101 -t 20 -R
```
```sh
bash -c 'echo -e "\n# RKCluster\n192.168.65.101\tthebatman\n192.168.65.102\tthegodfather" >> /etc/hosts'
```
```sh
mount -t nfs thebatman:/mpicluster /mpicluster
```
```sh
df -hT /mpicluster
```
```sh
mkdir /mpicluster
tee -a /etc/fstab <<EOF
thebatman:/mpicluster    /mpicluster    nfs    defaults    0 0
EOF
sudo mount -a
```
```sh
chmod 700 .ssh/
chmod 600 .ssh/authorized_keys
```
```sh
sudo adduser --uid <same at master node> --no-create-home --disabled-login --gecos "" mpi
```
```sh
ssh mpi@thebatman
```
```sh
exit
```
```sh
reboot
```

**Master Node**
```sh
sudo su
```
```sh
ssh mpi@thegodfather
```
```sh
exit
```
```sh
cat > $HOME/helloworld.c << EOF
#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
int myrank, nprocs;

MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

printf("Hello, World! from processor %d of %d\n", myrank, nprocs);

MPI_Finalize();
return 0;
}
EOF
```
```sh
mpicc helloworld.c -o hello_world
```
```sh
mpiexec --oversubscribe -n 20 -host 192.168.65.101,192.168.65.102 $HOME/hello_world
```

<br>

## Configure Multi-Node User

**Master Node**
```sh
sudo su
```
```sh
USER=
```
```sh
adduser --home /cluster --gecos "" "$USER"
```
```sh
usermod -aG sudo "$USER"
```
```sh
id "$USER"
```
```sh
ssh mpi@thegodfather # try ssh thegodfather
```

**Slave Node**
```sh
USER=
```
```sh
sudo adduser --uid <> --no-create-home --disabled-login --gecos "" "$USER"
```
```sh
usermod -aG sudo "$USER"
```

<br>

## Configure BEAGLE & BEAST
BEAGLE v4.0.1 (Oct 13, 2023) - https://github.com/beagle-dev/beagle-lib
```sh
cd && wget https://github.com/beagle-dev/beagle-lib/archive/refs/tags/v4.0.1.tar.gz
tar -zxvf v4.0.1.tar.gz && cd beagle-lib-4.0.1
mkdir build/ && cd build/
cmake -DCMAKE_C_COMPILER=gcc-9 -DCMAKE_CXX_COMPILER=g++-9 -DCMAKE_INSTALL_PREFIX:PATH=/mpicluster/opt/beagle ..
sudo make install
cd && rm -rf v4.0.1.tar.gz beagle-lib-4.0.1
```
BEAST v1.10.4 (Nov 14, 2018) - https://github.com/beast-dev/beast-mcmc
```sh
cd && wget https://github.com/beast-dev/beast-mcmc/releases/download/v1.10.4/BEASTv1.10.4.tgz
tar -zxvf BEASTv1.10.4.tgz && rm -rf BEASTv1.10.4.tgz
sudo mv BEASTv1.10.4/ /mpicluster/opt/BEASTv1.10.4/
sudo chown -R root:root /mpicluster/BEASTv1.10.4/
```
```sh
export LD_LIBRARY_PATH=/mpicluster/opt/beagle/lib:$LD_LIBRARY_PATH
export PATH=/mpicluster/opt/BEASTv1.10.4/bin:$PATH
beast -beagle_info
```
BEASTv1.10.5pre_thorney_0.1.2 (Sep 1, 2021) - https://github.com/beast-dev/beast-mcmc
```sh
cd && wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0-beta5/BEAST_X_v10.5.0-beta5.tgz
tar -zxvf BEAST_X_v10.5.0-beta5.tgz && rm -rf BEAST_X_v10.5.0-beta5.tgz
sudo mv BEASTv10.5.0/ /mpicluster/opt/BEASTv10.5.0-beta5/
sudo chown -R root:root /mpicluster/opt/BEASTv10.5.0-beta5/
```
```sh
export LD_LIBRARY_PATH=/mpicluster/opt/beagle/lib:$LD_LIBRARY_PATH
export PATH=/mpicluster/opt/BEASTv10.5.0-beta5/bin:$PATH
beast -beagle_info
```
BEAST v10.5.0 (July 2, 2025) - https://github.com/beast-dev/beast-mcmc
```sh
cd && wget https://github.com/beast-dev/beast-mcmc/releases/download/v10.5.0/BEAST_X_v10.5.0.tgz
tar -zxvf BEAST_X_v10.5.0.tgz && rm -rf BEAST_X_v10.5.0.tgz
sudo mv BEASTv10.5.0/ /mpicluster/opt/BEASTv10.5.0/
sudo chown -R root:root /mpicluster/opt/BEASTv10.5.0/
```
```sh
export LD_LIBRARY_PATH=/mpicluster/opt/beagle/lib:$LD_LIBRARY_PATH
export PATH=/mpicluster/opt/BEASTv10.5.0/bin:$PATH
beast -beagle_info
```

## etc
```sh
cat > $HOME/gen_hostfile.sh << EOF
#!/bin/bash

HOSTFILE="hostfile_auto"

HOSTS=("192.168.65.101" "192.168.65.102")

for HOST in "${HOSTS[@]}"; do
  CPU_CORES=$(ssh $HOST "nproc --all")
  GPU_COUNT=$(ssh $HOST "nvidia-smi --query-gpu=count --format=csv,noheader 2>/dev/null" || echo "0")
  RAM_GB=$(ssh $HOST "free -g | grep Mem | awk '{print \$2}'")

  SLOTS=$CPU_CORES
  if [ "$GPU_COUNT" -gt 0 ]; then
    SLOTS=$((CPU_CORES + GPU_COUNT))
  fi

  echo "$HOST slots=$SLOTS # CPU=$CPU_CORES GPU=$GPU_COUNT RAM=${RAM_GB}GB" >> $HOSTFILE
done
EOF
```
```sh
mpirun --hostfile hostfile bash -c "echo \$HOSTNAME; nproc; echo"
```
```sh

```
```sh

```
```sh

```
```sh

```
```sh
tee /etc/mpi_hostfile <<EOF
batman slots=$(nproc --all)
godfather slots=$(ssh godfather nproc --all)
EOF
```
```sh
tee /etc/ld.so.conf.d/mpi.conf <<EOF
/usr/lib/x86_64-linux-gnu/openmpi/lib
EOF
sudo ldconfig
```
```sh
sudo -u nobody mpirun --hostfile /etc/mpi_hostfile -np $(($(nproc)*3)) hostname
```
```sh
apt install -y ocl-icd-opencl-dev nvidia-openmpi
```
```sh
sudo tee /etc/openmpi/openmpi-mca-params.conf <<EOF
opal_cuda_support=1
EOF
```
```sh
sudo -u nobody mpirun --hostfile /etc/mpi_hostfile --mca opal_cuda_support 1 /usr/sbin/ompi_info | grep -i cuda
```
```sh
apt install -y htop glances
```
```sh
glances --webserver --bind 0.0.0.0 -p 61208
```
