
| #                 | CPU                   | GPU                            | RAM        | STORAGE                 |
| ----------------- | --------------------- | ------------------------------ | ---------- | ----------------------- |
| **THE GODFATHER** | Intel Core i9-10900   | NVIDIA GeForce RTX 2060 Rev. A | 128GB DDR4 | SSD NVMe 500GB, HDD 2TB |
| **THE BATMAN**    | Intel Core i7-10700KF | NVIDIA GeForce GTX 1660 Ti     | 32GB DDR4  | SSD NVMe 500GB, HDD 1TB |
| **ALIEN**         | Intel Core i7-9750H   | NVIDIA GeForce RTX 2080 Max-Q  | 16GB DDR4  | SSD NVMe 1TB            |
| **BLACK SWAN**    | Intel Core i3-4005U   |                                | 4GB DDR3   | HDD 1TB                 |
| **JOJO RABIT**    | Intel Core i3-3720QM  | AMD Radeon HD 7570M            | 8GB DDR3   | HDD 750GB               |

```sh
inxi -F # general infos
nproc # CPU threads
nvidia-settings -q CUDACores -t # CUDA cores
```
