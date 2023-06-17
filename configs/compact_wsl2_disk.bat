@echo off

:: Desligar o WSL
wsl.exe --shutdown

:: Buscar localização do ext4.vhdx
for /R "C:\Users" %%i in (ext4.vhdx) do (
    set "vhdx_path=%%i"
)

:: Executar Diskpart com comandos salvos em um arquivo txt temporário
echo select vdisk file="%vhdx_path%"> temp.txt
echo attach vdisk readonly>> temp.txt
echo compact vdisk>> temp.txt
echo detach vdisk>> temp.txt

:: Executar Diskpart com comandos do arquivo txt temporário
diskpart /s temp.txt

:: Remover o arquivo txt temporário
del temp.txt

:: Concluir o script
echo Done.
pause
exit
