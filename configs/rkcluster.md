## RKCLUSTER

| Hostname      | IP Address     | Role         | CPU Cores/Threads              | GPU                                                 | RAM        | Storage            |
| ------------- | -------------- | ------------ | ------------------------------ | --------------------------------------------------- | ---------- | ------------------ |
| GLADIATOR     | 192.168.65.100 | NAS/Storage  | ARMADA 388 (2c/2t)             |                                                     | 2GB DDR3   | HDD 16TB           |
| THE BATMAN    | 192.168.65.101 | Primary Node | Intel Core i7-10700KF (8c/16t) | NVIDIA GeForce RTX 2060 Rev. A, 6GB GDDR6 (1,920cc) | 32GB DDR4  | SSD 500GB, HDD 1TB |
| THE GODFATHER | 192.168.65.102 | Worker Node  | Intel Core i9-10900 (10c/20t)  | NVIDIA GeForce GTX 1660 Ti, 6GB GDDR6 (1,536cc)     | 128GB DDR4 | SSD 500GB, HDD 2TB |
| LOTR          | 192.168.65.103 | Worker Node  | Intel Xeon W5-3423 (12c/24t)   | NVIDIA Quadro RTX 6000, 48GB DDR6 (4,608cc)         | 256GB DDR5 | SSD 512GB, HDD 2TB |

### Configure TP Link Wi-Fi Router
```
http://tplinkwifi.net
```
#### Advanced > Network > LAN > DCHP Server
- IP Version: `IPv4`
- IP Address: `192.168.65.1`
- Subnet Mask: `255.255.255.0`
- IP Address Pool: `192.168.65.100` - `192.168.65.199`
- IP Address: `192.168.65.1`
- Default Gateway: `192.168.50.1`
- Primary DNS: `8.8.8.8`
- Secondary DNS: `8.8.4.4`

#### Advanced > Network > LAN > Address Reservation
- GLADIATOR: `192.168.65.100`
- THE BATMAN: `192.168.65.101`
- THE GODFATHER: `192.168.65.102`

#### Advanced > Wireless > Wireless Settings
- Network Name (SSID): `RKhour0-Bioinfo`
- Password: `password`

### Configure NAS/Storage (Western Digital My Cloud Expert Series EX4100)
```
http://192.168.65.100
```
#### Settings > Network > 
- IPv4 Network Mode: `DHCP`
- IPv6 Network Mode: `Off`
- FTP Access: `ON`
- NFS Access: `ON`
- SSH: `ON`

#### Shares > New Sharing
- Share Name: `gladiator`
- Public: `ON`
- Share Access > NFS: `ON`
- NFS Settings > Host: `*`
- NFS Settings > Write: `ON`

### NFS Client Configuration on Cluster Nodes
```sh
sudo su
```
```sh
apt update && apt install nfs-common -y
```
```sh
mkdir -p /mnt/gladiator
```
```sh
mount -t nfs 192.168.65.100:/mnt/HD/HD_a2/gladiator /mnt/gladiator -o rw,noatime,user
```
```sh
nano /etc/fstab
```
```sh
tee -a /etc/fstab <<EOF
192.168.65.100:/mnt/HD/HD_a2/gladiator  /mnt/gladiator  nfs  rw,noatime,auto,nofail  0  0
EOF
```
```sh
sudo mount -a
```
```sh
df -hT /mnt/gladiator
```
```sh
reboot
```

### Configure Static IP on Cluster Nodes
```sh
apt-get install iperf3
```
```sh
iperf3 -s # Primary Node
```
```sh
iperf3 -c 192.168.65.101 -t 20 # Worker Node
```

## Configure Passwordless SSH Between Nodes
```sh
sudo su
```
```sh
mkdir -p /etc/cluster-ssh
ssh-keygen -t ed25519 -N "" -f /etc/cluster-ssh/id_ed25519
chmod 750 /etc/cluster-ssh
chmod 600 /etc/cluster-ssh/id_ed25519*
```
```sh
for node in godfather; do
  ssh-copy-id -i /etc/cluster-ssh/id_ed25519.pub $node
  ssh $node "sudo mkdir -p /etc/cluster-ssh && sudo chmod 750 /etc/cluster-ssh"
  scp /etc/cluster-ssh/id_ed25519.pub $node:/etc/cluster-ssh/
done
```
```sh
apt update && apt install -y mpich openmpi-bin libopenmpi-dev
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
sudo ssh -i /etc/cluster-ssh/id_ed25519 godfather hostname
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
```sh
cluster-glances() {
  for node in batman godfather; do
    echo "=== $node ==="
    ssh $node "top -bn1 | head -5"
  done
}
```

cluster-htop() {
  ssh -X batman htop
  # Ou para modo texto:
  for node in batman godfather; do
    echo "=== $node ==="
    ssh $node "free -h && mpstat -P ALL"
  done
}

cluster-nvidia() {
  for node in batman godfather; do
    echo "=== $node ==="
    ssh $node "nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv"
  done
}
