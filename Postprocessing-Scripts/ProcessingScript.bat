echo off

set expDate=%1
set sid=%2
set parentDir=D:\Dropbox (HMS)\2P Data\Behavior Vids\%expDate%\

if not exist "%parentDir%" mkdir "%parentDir%"
if not exist "%parentDir%\_Movies" mkdir "%parentDir%_Movies"

echo Downloading files from /n/scratch2...

"C:\Program Files (x86)\WinSCP\WinSCP.exe" /log="C:\Users\Wilson Lab\Desktop\WinSCP.log" /script="C:\Users\Wilson Lab\Desktop\testscript.txt" /parameter "%expDate%" "%sid%" "%parentDir%" 


if %ERRORLEVEL% neq 0 (
	echo %ERRORLEVEL%
	echo WinSCP reported error...aborting script
	exit /B 0
) else (
	echo Remote transfer completed successfully...attempting Handbrake conversion
)


for %%f in ("%parentDir%_Movies\*tid*.avi") do (
	echo %%f
	echo "%%~df%%~pf%%~nf.mp4"
	"C:\Program Files\Handbrake\HandBrakeCLI.exe" --preset-import-file "C:\Users\Wilson Lab\Desktop\BehaviorVidConversionPresets.json" -i "%%f" -o "%%~df%%~pf%%~nf.mp4" 
)

set serverPath=\\research.files.med.harvard.edu\Neurobio\Wilson Lab\Michael\FicTrac Vids\%expDate%
echo %serverPath%
if not exist %serverPath% mkdir %serverPath%
echo "%parentDir%%_Movies\*.mp4"
echo "%serverPath%\"
copy "%parentDir%_Movies\*.mp4" "%serverPath%\"