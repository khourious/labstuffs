## Temporarily Replacing Utilman
**During Windows setup:** press Shift + F10 to open Command Prompt
```batch
cd /d C:\Windows\System32
ren utilman.exe utilman.exe.bak
ren cmd.exe utilman.exe
wpeutil reboot

```

**At login screen:** click the "Ease of Access" icon which now opens Command Prompt
```batch
net user Administrator /active:yes
shutdown /r /t 0

```

## System Restoration
**During Windows setup:** press Shift + F10 to open Command Prompt
```batch
cd /d C:\Windows\System32
ren utilman.exe cmd.exe 
ren utilman.exe.bak utilman.exe
net user Administrator /active:no
wpeutil reboot

```
