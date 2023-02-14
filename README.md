## Step-by-step guide to installing bioinformatics configurations

- [macOS](configs/macOS.md)
- [Ubuntu 20.04 - Native Linux](configs/Linux.md)
- [Ubuntu 20.04 - Linux Subsystem for Windows 2 (WSL2)](configs/Windows_WSL2.md)

## How to compact a WSL2 Disk

Shutdown WSL2:

```sh
wsl.exe --shutdown
```

Open Command Prompt as administrator: press `Windows + R` to open the Run Command Window, type `cmd` and press `Ctrl + Shift + Enter`.

Download the `compact_wsl2_disk.bat` script:

```cmd
wget https://github.com/khourious/labstuffs/raw/master/configs/compact_wsl2_disk.bat
```

Find the WSL2 disk location:

```cmd
WHERE.exe /R C:\Users ext4.vhdx
```

Run the script by typing `compact_wsl2_disk.bat` followed by the path of the ext4.vhdx file:

```cmd
compact_wsl2_disk.bat <path_to_ext4_vhdx_file>
```
