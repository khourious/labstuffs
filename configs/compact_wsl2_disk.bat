@echo off

echo Shutting down WSL...
wsl.exe --shutdown

echo Searching for the ext4.vhdx file...
for /R "C:\Users" %%i in (ext4.vhdx) do (
    set "vhdx_path=%%i"
)

echo The ext4.vhdx file was found at: %vhdx_path%

echo Preparing Diskpart commands...

echo Selecting vdisk file...
echo select vdisk file="%vhdx_path%"> temp1.txt
diskpart /s temp1.txt
del temp1.txt

echo Attaching vdisk readonly...
echo attach vdisk readonly> temp2.txt
diskpart /s temp2.txt
del temp2.txt

echo Compacting vdisk...
echo compact vdisk> temp3.txt
diskpart /s temp3.txt
del temp3.txt

echo Detaching vdisk...
echo detach vdisk> temp4.txt
diskpart /s temp4.txt
del temp4.txt

echo Done.
pause
exit
