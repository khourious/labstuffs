## RKCLUSTER

| Hostname      | IP Address     | Role         | CPU Cores/Threads              | GPU                                                 | RAM        | Storage            |
| ------------- | -------------- | ------------ | ------------------------------ | --------------------------------------------------- | ---------- | ------------------ |
| GLADIATOR     | 192.168.65.100 | NAS/Storage  | ARMADA 388 (2c/2t)             |                                                     | 2GB DDR3   | HDD 16TB           |
| THE BATMAN    | 192.168.65.101 | Primary Node | Intel Core i7-10700KF (8c/16t) | NVIDIA GeForce RTX 2060 Rev. A, 6GB GDDR6 (1,920cc) | 32GB DDR4  | SSD 500GB, HDD 1TB |
| THE GODFATHER | 192.168.65.102 | Worker Node  | Intel Core i9-10900 (10c/20t)  | NVIDIA GeForce GTX 1660 Ti, 6GB GDDR6 (1,536cc)     | 128GB DDR4 | SSD 500GB, HDD 2TB |
| LOTR          | 192.168.65.103 | Worker Node  | Intel Xeon W5-3423 (12c/24t)   | NVIDIA Quadro RTX 6000, 48GB DDR6 (4,608cc)         | 256GB DDR5 | SSD 512GB, HDD 2TB |

### Configure Static IP on Cluster Nodes
```sh
sudo su
```
```sh
tee /etc/netplan/00-installer-config.yaml <<EOF && chmod 600 /etc/netplan/00-installer-config.yaml
network:
  version: 2
  ethernets:
    $(ip route | grep default | awk '{print $5}' | head -n 1):
      addresses: [$(hostname -I | awk '{print $1}')/24]
      routes:
        - to: default
          via: $(ip route | grep default | awk '{print $3}' | head -n 1)
      nameservers:
        addresses: [8.8.8.8, 8.8.4.4]
EOF
```
```sh
netplan apply
```
```sh
ip route show
```
```sh
ping 8.8.8.8
```
```sh
ufw allow 5201/tcp
```
```sh
apt-get install iperf3
```
```sh
iperf3 -s # Primary Node
```
```sh
iperf3 -c 192.168.65.101 -t 20 # Worker Node
```

### NFS Server Configuration on NAS/Storage
```sh
sudo su
```
```sh
apt update && apt install -y nfs-kernel-server
```
```sh
mkdir -p /mnt/cluster
chown nobody:nogroup /mnt/cluster
chmod 1777 /mnt/cluster
```
```sh
tee /etc/exports <<'EOF'
/mnt/cluster 192.168.65.101(rw,sync,no_subtree_check)
/mnt/cluster 192.168.65.102(rw,sync,no_subtree_check)
EOF
sudo exportfs -av
sudo systemctl restart nfs-kernel-server
sudo systemctl enable nfs-kernel-server
```
```sh
showmount -e localhost
```
```sh
sudo ufw allow from 192.168.65.0/24 to any port nfs
```

### NFS Client Configuration on Cluster Nodes
```sh
sudo su
```
```sh
apt install -y nfs-common
mkdir -p /mnt/nas
```
```sh
mount 192.168.65.100:/mnt/cluster /mnt/nas
```
```sh
df -hT /mnt/nas
```
```sh
echo "192.168.65.100:/mnt/cluster  /mnt/nas  nfs  defaults,noatime,vers=4.2  0  0" | sudo tee -a /etc/fstab
```
```sh
reboot
```
```sh
sudo su
```
```sh
mount -av
```
```sh
ls /mnt/nas
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
