## RKCLUSTER

| Hostname      | IP Address     | Role         | CPU Cores/Threads              | GPU                                                 | RAM        | Storage            |
| ------------- | -------------- | ------------ | ------------------------------ | --------------------------------------------------- | ---------- | ------------------ |
| GLADIATOR     | 192.168.65.100 | NAS/Storage  | ARMADA 388 (2c/2t)             |                                                     | 2GB DDR3   | HDD 16TB           |
| THE BATMAN    | 192.168.65.101 | Master       | Intel Core i7-10700KF (8c/16t) | NVIDIA GeForce RTX 2060 Rev. A, 6GB GDDR6 (1,920cc) | 32GB DDR4  | SSD 500GB, HDD 1TB |
| THE GODFATHER | 192.168.65.102 | Slave        | Intel Core i9-10900 (10c/20t)  | NVIDIA GeForce GTX 1660 Ti, 6GB GDDR6 (1,536cc)     | 128GB DDR4 | SSD 500GB, HDD 2TB |
| LOTR          | 192.168.65.103 | Slave        | Intel Xeon W5-3423 (12c/24t)   | NVIDIA Quadro RTX 6000, 48GB DDR6 (4,608cc)         | 256GB DDR5 | SSD 512GB, HDD 2TB |

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

#### Settings > Network > Network Services
- IPv4 Network Mode: `Static`
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

### NFS Client Configuration

#### Master Node

Open the terminal

XXX
```sh
sudo su
```

Install XXX
```sh
apt update && apt install nfs-common -y
```

XXX and mount
```sh
mkdir -p /mnt/gladiator
tee -a /etc/fstab <<EOF
192.168.65.100:/mnt/HD/HD_a2/gladiator  /mnt/gladiator  nfs  rw,noatime,auto,nofail  0  0
EOF
sudo mount -a
```

XXX
```sh
df -hT /mnt/gladiator
```

XXX
```sh
reboot
```

### Network Performance Test

#### Master Node

XXX
```sh
sudo su
```

Install XXX
```sh
apt-get install iperf3
```

XXX
```sh
iperf3 -s
```

#### Slave Node

XXX
```sh
sudo su
```

Install XXX
```sh
apt-get install iperf3
```

XXX
```sh
iperf3 -c 192.168.65.101 -t 20 -R
```

#### Configure Passwordless SSH Between Nodes

#### Master Node

XXX
```sh
sudo su
```

Install OpenMPI
```sh
apt update && apt install -y mpich openmpi-bin libopenmpi-dev
```

Create a new user called `mpiuser`
```sh
adduser mpiuser
```

Give sudo access to the newly created user
```sh
usermod -aG sudo mpiuser
```

Add the IP addresses of master (thebatman) and slave (thegodfather) nodes to the host file
```sh
bash -c 'echo -e "\n# RKCluster\n192.168.65.101\tthebatman\n192.168.65.102\tthegodfather" >> /etc/hosts'
```

```sh
su mpiuser
```

Generate a key
```sh
ssh-keygen -t ed25519 -N "" -f /home/mpiuser/.ssh/id_ed25519
```

Copy the generated key onto the slave node (thegodfather)
```sh
ssh-copy-id mpiuser@192.168.65.102
```




```sh
cat > /root/mpi_hello_world.c << 'EOF'
#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();
}
EOF
```





XXX
```sh
ssh-keygen -t ed25519 -N "" -f /root/.ssh/id_ed25519
chmod 750 /root/.ssh/
chmod 600 /root/.ssh/id_ed25519*
cp /root/.ssh/id_ed25519.pub /root/.ssh/authorized_keys
```

XXX
```sh
ssh-copy-id -i /root/.ssh/id_ed25519.pub lpmor22@192.168.65.102
```

XXX
```sh
scp /root/.ssh/authorized_keys lpmor22@192.168.65.102:/tmp/
ssh -t lpmor22@192.168.65.102 "sudo mv /tmp/authorized_keys /root/.ssh/ && \
  sudo chmod 700 /root/.ssh/ && \
  sudo chmod 600 /root/.ssh/authorized_keys && \
  sudo chrown root:root /root/.ssh/authorized_keys"
```

XXX
```sh
ssh root@192.168.65.102
```

```sh
exit
```

```sh
mpiexec --oversubscribe -n 20 -host 192.168.65.101,192.168.65.102 /root/mpi_hello_world
```


#### Slave Node

XXX
```sh
sudo su
```

Install OpenMPI
```sh
apt update && apt install -y mpich openmpi-bin libopenmpi-dev
```

```sh
cat > /root/mpi_hello_world.c << 'EOF'
#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();
}
EOF
```

Compile the mpi program's file
```sh
mpicc mpihelloworld.c -o mpi_hello_world
```

XXX
```sh
sed -i 's/^#\?PubkeyAuthentication .*/PubkeyAuthentication yes/' /etc/ssh/sshd_config
sed -i 's|^#\?AuthorizedKeysFile .*|AuthorizedKeysFile .ssh/authorized_keys|' /etc/ssh/sshd_config
sed -i '/^Match User root/,/^    PasswordAuthentication no/d' /etc/ssh/sshd_config
echo -e "\nMatch User root\n    PasswordAuthentication no\n    PermitRootLogin prohibit-password" >> /etc/ssh/sshd_config
```

XXX
```sh
sudo systemctl restart ssh
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

```
```sh

```


sudo adduser mpiuser



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



#### Slave Node

XXX
```sh
sudo su
```

Install OpenMPI
```sh
apt update && apt install -y mpich openmpi-bin libopenmpi-dev
```

```sh
cat > /root/mpi_hello_world.c << 'EOF'
#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();
}
EOF
```

Compile the mpi program's file
```sh
mpicc mpihelloworld.c -o mpi_hello_world
```

XXX
```sh
sed -i 's/^#\?PubkeyAuthentication .*/PubkeyAuthentication yes/' /etc/ssh/sshd_config
sed -i 's|^#\?AuthorizedKeysFile .*|AuthorizedKeysFile .ssh/authorized_keys|' /etc/ssh/sshd_config
sed -i '/^Match User root/,/^    PasswordAuthentication no/d' /etc/ssh/sshd_config
echo -e "\nMatch User root\n    PasswordAuthentication no\n    PermitRootLogin prohibit-password" >> /etc/ssh/sshd_config
```

XXX
```sh
sudo systemctl restart ssh
```

#### Master Node

XXX
```sh
sudo su
```

Install OpenMPI
```sh
apt update && apt install -y mpich openmpi-bin libopenmpi-dev
```

```sh
cat > /root/mpi_hello_world.c << 'EOF'
#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();
}
EOF
```

Create a new user called `mpiuser`
```sh
adduser mpiuser
```

Give sudo access to the newly created user
```sh
usermod -aG sudo mpiuser
```

Add the IP addresses of master and slave nodes to the host file
```sh
bash -c 'echo -e "\n# RKCluster\n192.168.65.101\tthebatman\n192.168.65.102\tthegodfather" >> /etc/hosts'
```

XXX
```sh
ssh-keygen -t ed25519 -N "" -f /root/.ssh/id_ed25519
chmod 750 /root/.ssh/
chmod 600 /root/.ssh/id_ed25519*
cp /root/.ssh/id_ed25519.pub /root/.ssh/authorized_keys
```

XXX
```sh
ssh-copy-id -i /root/.ssh/id_ed25519.pub lpmor22@192.168.65.102
```

XXX
```sh
scp /root/.ssh/authorized_keys lpmor22@192.168.65.102:/tmp/
ssh -t lpmor22@192.168.65.102 "sudo mv /tmp/authorized_keys /root/.ssh/ && \
  sudo chmod 700 /root/.ssh/ && \
  sudo chmod 600 /root/.ssh/authorized_keys && \
  sudo chrown root:root /root/.ssh/authorized_keys"
```

XXX
```sh
ssh root@192.168.65.102
```

```sh
exit
```

```sh
mpiexec --oversubscribe -n 20 -host 192.168.65.101,192.168.65.102 /root/mpi_hello_world
```



https://www.geeksforgeeks.org/linux-unix/creating-an-mpi-cluster/
https://github.com/adeen-atif/MPI-Cluster
