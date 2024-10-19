## setup

| #                 | CPU                   | GPU                            | RAM        | STORAGE                 |
| ----------------- | --------------------- | ------------------------------ | ---------- | ----------------------- |
| **THE GODFATHER** | Intel Core i9-10900   | NVIDIA GeForce RTX 2060 Rev. A | 128GB DDR4 | SSD NVMe 500GB, HDD 2TB |
| **THE BATMAN**    | Intel Core i7-10700KF | NVIDIA GeForce GTX 1660 Ti     | 32GB DDR4  | SSD NVMe 500GB, HDD 1TB |
| **ALIEN**         | Intel Core i7-9750H   | NVIDIA GeForce RTX 2080 Max-Q  | 16GB DDR4  | SSD NVMe 1TB            |
| **JOJO RABBIT**   | Intel Core i7-3720QM  | AMD Radeon HD 7570M            | 8GB DDR3   | HDD 750GB               |

```sh
inxi -F # general infos
nproc # CPU threads
nvidia-settings -q CUDACores -t # CUDA cores
sudo dmidecode --type 17 # RAM infos
```

## BEAST Benchmank

file: benchmark2.xml <br>
no. taxa: 62 <br>
no. characters: 10869 <br>
chain length: 1000000

```sh
wget https://raw.githubusercontent.com/beast-dev/beast-mcmc/master/examples/Benchmarks/benchmark2.xml
beast -beagle_CPU benchmark2.xml
beast -beagle_GPU benchmark2.xml
```

|                   | CPU                       | GPU                        |
| ----------------- | ------------------------- | -------------------------- |
| **THE GODFATHER** | 6.5537 minutes            | 2.3248166666666665 minutes |
| **THE BATMAN**    | 8.121383333333334 minutes | 2.0611166666666665 minutes |
| **ALIEN**         |                           |                            |
| **JOJO RABBIT**   |                           |                            |
