# RKCluster

<br>

- [Khourious Nodes](#khourious-nodes)
- [Configure TP Link Wi-Fi Router](#configure-tp-link-wi-fi-router)
- [Configure NAS (Western Digital My Cloud Expert Series EX4100)](#configure-nas-western-digital-my-cloud-expert-series-ex4100)
- [Configure Passwordless SSH Between Nodes](#configure-passwordless-ssh-between-nodes)
- [Configure NFS Client](#configure-nfs-client)
- [Configure Multi-Node User](#configure-multi-node-user)

<br>

## Khourious Nodes

| Hostname      | IP Address     | Role   | CPU Cores/Threads              | GPU                                                 | RAM        | Storage            |
| ------------- | -------------- | ------ | ------------------------------ | --------------------------------------------------- | ---------- | ------------------ |
| GLADIATOR     | 192.168.65.100 | NAS    | ARMADA 388 (2c/2t)             |                                                     | 2GB DDR3   | HDD 16TB           |
| THE MATRIX    | 192.168.65.101 | Master | Intel Core i5-12500T (6c/12t)  |                                                     | 8GB DDR5   | SSD 512GB          |
| THE BATMAN    | 192.168.65.102 | Slave  | Intel Core i7-10700KF (8c/16t) | NVIDIA GeForce RTX 2060 Rev. A, 6GB GDDR6 (1,920cc) | 32GB DDR4  | SSD 500GB, HDD 1TB |
| THE GODFATHER | 192.168.65.103 | Slave  | Intel Core i9-10900 (10c/20t)  | NVIDIA GeForce GTX 1660 Ti, 6GB GDDR6 (1,536cc)     | 128GB DDR4 | SSD 500GB, HDD 2TB |
| LOTR          | 192.168.65.104 | Slave  | Intel Xeon W5-3423 (12c/24t)   | NVIDIA Quadro RTX 6000, 48GB DDR6 (4,608cc)         | 256GB DDR5 | SSD 512GB, HDD 2TB |

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
- THE MATRIX: `192.168.65.101`
- THE BATMAN: `192.168.65.102`
- THE GODFATHER: `192.168.65.103`

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

## Configure Passwordless SSH Between Nodes

**Master Node**
```sh
sudo apt update -y && \
sudo apt upgrade -y && \
sudo apt install -y \
  build-essential \
  iperf3 \
  libmunge2 \
  libmunge-dev \
  libopenmpi-dev \
  mpich \
  munge \
  nfs-client \
  nfs-common \
  nfs-server \
  ntp \
  openmpi-bin \
  openssh-server \
  slurm-wlm
```
```sh
sudo bash -c 'echo -e "\n192.168.65.101\tthematrix\n192.168.65.102\tthebatman\n192.168.65.103\tthegodfather" >> /etc/hosts'
```
```sh
iperf3 -s
```
```sh
sudo mkdir /cluster
```
```sh
sudo tee -a /etc/exports <<EOF
/cluster *(rw,sync,no_subtree_check,no_root_squash)
EOF
```
```sh
sudo exportfs -a
```
```sh
sudo systemctl restart nfs-kernel-server
```
```sh
sudo exportfs -v
```
```sh
sudo tee -a /etc/fstab <<EOF
UUID='$(sudo blkid -s UUID -o value /dev/sda1)'  /home  ext4  defaults  0  2
EOF
```
```sh
sudo mount -a
```
```sh
df -hT /cluster
```
```sh
sudo adduser --home /mpicluster/home/slurm --gecos "" slurm
```
```sh
sudo usermod -aG sudo slurm
```
```sh
id slurm
```
```sh
su slurm
```
```sh
ssh-keygen
```
```sh
cat $HOME/.ssh/id_ed25519.pub >> $HOME/.ssh/authorized_keys
```
```sh
SLAVE=
```
```sh
ssh "$SLAVE"
```

**Slave Node**
```sh
sudo apt update -y && \
sudo apt upgrade -y && \
sudo apt install -y \
  build-essential \
  iperf3 \
  libmunge2 \
  libmunge-dev \
  libopenmpi-dev \
  mpich \
  munge \
  nfs-client \
  nfs-common \
  nfs-server \
  ntp \
  openmpi-bin \
  openssh-server \
  slurm-wlm
```
```sh
sudo bash -c 'echo -e "\n192.168.65.101\tthematrix\n192.168.65.102\tthebatman\n192.168.65.103\tthegodfather" >> /etc/hosts'
```
```sh
iperf3 -c thematrix -t 20 -R
```
```sh
sudo mkdir /cluster
```
```sh
sudo tee -a /etc/exports <<EOF
/cluster *(rw,sync,no_subtree_check,no_root_squash)
EOF
```
```sh
sudo exportfs -a
```
```sh
sudo systemctl restart nfs-kernel-server
```
```sh
sudo exportfs -v
```
```sh
sudo tee -a /etc/fstab <<EOF
thematrix:/cluster  /cluster  nfs  defaults  0  0
EOF
```
```sh
sudo mount -a
```
```sh
df -hT /cluster
```
```sh
sudo chown -R root:root /cluster && \
sudo find /cluster -type d -exec chmod 755 {} \; && \
sudo find /cluster -type f -exec chmod 644 {} \;
```
```sh
sudo chown -R root:root /cluster/home && \
sudo find /cluster/home -type d -exec chmod 755 {} \; && \
sudo find /cluster/home -type f -exec chmod 644 {} \;
```
```sh
sudo adduser --home /cluster/home/slurm --gecos "" slurm --uid <> --gid <>
```
```sh
sudo usermod -aG sudo slurm
```
```sh
cat /etc/passwd
```
```sh
su slurm
```
```sh
ssh thematrix
```
```sh
exit
```
```sh
sudo reboot
```

**Master Node**
```sh
su slurm
```
```sh
SLAVE=
```
```sh
ssh "$SLAVE"
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
mpiexec --oversubscribe -n 20 -host thematrix,"$SLAVE" $HOME/hello_world
```
```sh
munge -n | unmunge | grep STATUS
```
```sh
sudo /usr/sbin/mungekey --verbose --force
```
```sh
cp /etc/munge/munge.key /tmp/ && \
sudo scp /tmp/munge.key "$SLAVE":/etc/munge/munge.key
```
```sh
sudo chown -R munge: /etc/munge{,/munge.key} /var/{log,lib}/munge/ /run/munge/ && \
sudo chmod 0700 /etc/munge{,/munge.key} /var/{log,lib}/munge/ && \
sudo chmod 0755 /run/munge/
```
```sh
sudo mkdir -p /etc/systemd/system/munge.service.d/ && \
sudo tee -a /etc/systemd/system/munge.service.d/override.conf <<EOF
[Service]
Environment="OPTIONS=--syslog"
EOF
```
```sh
sudo systemctl daemon-reload && \
sudo systemctl enable munge && \
sudo systemctl restart munge
```
```sh
sudo systemctl status munge
```
```sh
ssh "$SLAVE" -t "munge -n | unmunge | grep STATUS"
```
```sh
sudo cp /etc/munge/munge.key /tmp/ && \
sudo scp /tmp/munge.key "$SLAVE":/etc/munge/munge.key
```
```sh
ssh "$SLAVE" -t "\
  sudo chown -R munge: /etc/munge{,/munge.key} /var/{log,lib}/munge/ /run/munge/ && \
  sudo chmod 0700 /etc/munge{,/munge.key} /var/{log,lib}/munge/ && \
  sudo chmod 0755 /run/munge/"
```
```sh
ssh "$SLAVE" -t"sudo mkdir -p /etc/systemd/system/munge.service.d/"
sudo scp /etc/systemd/system/munge.service.d/override.conf "$SLAVE":/etc/systemd/system/munge.service.d/
```
```sh
ssh "$SLAVE" -t "\
  sudo systemctl daemon-reload && \
  sudo systemctl enable munge && \
  sudo systemctl restart munge
```
```sh
ssh "$SLAVE" -t "sudo systemctl status munge"
```
```sh
munge -n | ssh "$SLAVE" unmunge
ssh "$SLAVE" -t "munge -n | ssh thematrix unmunge"
```
```sh
exit
```
```sh
slurmd -C
ssh "$SLAVE" -t "slurmd -C"
```
```sh
firefox /usr/share/doc/slurmctld/slurm-wlm-configurator.html
```
```sh
sudo tee -a /etc/slurm/slurm.conf <<EOF
# slurm.conf file generated by configurator.html.
# Put this file on all nodes of your cluster.
# See the slurm.conf man page for more information.
#
ClusterName=khourious
SlurmctldHost=thematrix
#
#DisableRootJobs=NO
#EnforcePartLimits=NO
#Epilog=
#EpilogSlurmctld=
#FirstJobId=1
#MaxJobId=67043328
#GresTypes=
#GroupUpdateForce=0
#GroupUpdateTime=600
#JobFileAppend=0
#JobRequeue=1
#JobSubmitPlugins=lua
#KillOnBadExit=0
#LaunchType=launch/slurm
#Licenses=foo*4,bar
#MailProg=/usr/bin/mail
#MaxJobCount=10000
#MaxStepCount=40000
#MaxTasksPerNode=512
MpiDefault=none
#MpiParams=ports=#-#
#PluginDir=
#PlugStackConfig=
#PrivateData=jobs
ProctrackType=proctrack/cgroup
#Prolog=
#PrologFlags=
#PrologSlurmctld=
#PropagatePrioProcess=0
#PropagateResourceLimits=
#PropagateResourceLimitsExcept=
#RebootProgram=
ReturnToService=1
SlurmctldPidFile=/run/slurmctld.pid
SlurmctldPort=6817
SlurmdPidFile=/run/slurmd.pid
SlurmdPort=6818
SlurmdSpoolDir=/var/lib/slurm/slurmd
SlurmUser=slurm
#SlurmdUser=root
#SrunEpilog=
#SrunProlog=
StateSaveLocation=/var/lib/slurm/slurmctld
SwitchType=switch/none
#TaskEpilog=
TaskPlugin=task/affinity,task/cgroup
#TaskProlog=
#TopologyPlugin=topology/tree
#TmpFS=/tmp
#TrackWCKey=no
#TreeWidth=
#UnkillableStepProgram=
#UsePAM=0
#
#
# TIMERS
#BatchStartTimeout=10
#CompleteWait=0
#EpilogMsgTime=2000
#GetEnvTimeout=2
#HealthCheckInterval=0
#HealthCheckProgram=
InactiveLimit=0
KillWait=30
#MessageTimeout=10
#ResvOverRun=0
MinJobAge=300
#OverTimeLimit=0
SlurmctldTimeout=120
SlurmdTimeout=300
#UnkillableStepTimeout=60
#VSizeFactor=0
Waittime=0
#
#
# SCHEDULING
#DefMemPerCPU=0
#MaxMemPerCPU=0
#SchedulerTimeSlice=30
SchedulerType=sched/backfill
SelectType=select/cons_tres
#
#
# JOB PRIORITY
#PriorityFlags=
#PriorityType=priority/basic
#PriorityDecayHalfLife=
#PriorityCalcPeriod=
#PriorityFavorSmall=
#PriorityMaxAge=
#PriorityUsageResetPeriod=
#PriorityWeightAge=
#PriorityWeightFairshare=
#PriorityWeightJobSize=
#PriorityWeightPartition=
#PriorityWeightQOS=
#
#
# LOGGING AND ACCOUNTING
#AccountingStorageEnforce=0
#AccountingStorageHost=
#AccountingStoragePass=
#AccountingStoragePort=
AccountingStorageType=accounting_storage/none
#AccountingStorageUser=
#AccountingStoreFlags=
#JobCompHost=
#JobCompLoc=
#JobCompParams=
#JobCompPass=
#JobCompPort=
JobCompType=jobcomp/none
#JobCompUser=
#JobContainerType=job_container/none
JobAcctGatherFrequency=30
JobAcctGatherType=jobacct_gather/none
SlurmctldDebug=info
SlurmctldLogFile=/var/log/slurm/slurmctld.log
SlurmdDebug=info
SlurmdLogFile=/var/log/slurm/slurmd.log
#SlurmSchedLogFile=
#SlurmSchedLogLevel=
#DebugFlags=
#
#
# POWER SAVE SUPPORT FOR IDLE NODES (optional)
#SuspendProgram=
#ResumeProgram=
#SuspendTimeout=
#ResumeTimeout=
#ResumeRate=
#SuspendExcNodes=
#SuspendExcParts=
#SuspendRate=
#SuspendTime=
#
#
# COMPUTE NODES
NodeName=thebatman NodeAddr=192.168.65.102 CPUs=16 RealMemory=31961 Sockets=1 CoresPerSocket=8 ThreadsPerCore=2 State=UNKNOWN
NodeName=thegodfather NodeAddr=192.168.65.103 CPUs=20 RealMemory=128679 Sockets=1 CoresPerSocket=10 ThreadsPerCore=2 State=UNKNOWN
PartitionName=khour Nodes=ALL Default=YES MaxTime=INFINITE State=UP
EOF
```
```sh
sudo scp /etc/slurm/slurm.conf "$SLAVE":/etc/slurm/
```
```sh
sudo mkdir -p /var/{log,spool/slurmctld,run}
sudo chown-R slurm:slurm /var/{log,spool,spool/slurmctld,run}
```
```sh
ssh "$SLAVE" -t "\
  sudo mkdir -p /var/{log,spool/slurmctld,run}
  sudo chown-R slurm:slurm /var/{log,spool,spool/slurmctld,run}"
```
```sh
sudo systemctl enable slurmctld && \
sudo systemctl restart slurmctld
```
```sh
sudo systemctl status slurmctld
```
```sh
ssh "$SLAVE" "\
  sudo systemctl enable slurmd && \
  sudo systemctl restart slurmd"
```
```sh
ssh "$SLAVE" "sudo systemctl status slurmd"
```
```sh
sinfo
```
```sh
scontrol show nodes
```
```sh
srun --pty hostname
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
SLAVE=
```
```sh
ssh "$SLAVE"
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


## Configure NFS Client

**Master Node**
```sh
sudo bash -c 'echo -e "\n192.168.65.100\tgladiator" >> /etc/hosts'
```
```sh
sudo mkdir -p /mnt/gladiator
```
```sh
sudo tee -a /etc/fstab <<EOF
gladiator:/mnt/HD/HD_a2/gladiator  /mnt/gladiator  nfs  rw,noatime,auto,nofail  0  0
EOF
```
```sh
sudo mount -a
```
```sh
df -hT /mnt/gladiator
```
```sh
sudo reboot
```
```sh
df -hT /mnt/gladiator
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
