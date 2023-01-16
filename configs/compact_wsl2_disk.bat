@echo off

setlocal enabledelayedexpansion

if "%1"=="" (
    echo ERROR: Please provide the path to the ext4.vhdx file as an argument.
    echo Example: compact_wsl2_disk.cmd "C:\Users\RKHOURI\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\ext4.vhdx"
    exit /b 1
)

set vhdx_path=%1

DISKPART

select vdisk file="%vhdx_path%"

attach vdisk readonly

compact vdisk

detach vdisk

endlocal
