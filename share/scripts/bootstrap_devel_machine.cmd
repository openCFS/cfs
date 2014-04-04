@ECHO OFF

::============================================================================
:: Set URL on CFS++ development server for bootstrap CMake script.
::============================================================================
set CFS_DEVEL_FTP_SERVER=ftp://lse17.e-technik.uni-erlangen.de:40065
set BOOTSTRAP_CMAKE_SCRIPT=bootstrap_devel_machine_win.cmake
set BOOTSTRAP_FTP_DIR=%CFS_DEVEL_FTP_SERVER%/scripts
set BOOTSTRAP_DOWNLOAD_URL=%BOOTSTRAP_FTP_DIR%/%BOOTSTRAP_CMAKE_SCRIPT%
set BOOTSTRAP_PCKG_MIRROR=%CFS_DEVEL_FTP_SERVER%/binaries/mingw_env/packages

::============================================================================
:: Set version and name for CMake executables.
::============================================================================
set CMAKE_MAJOR=2
set CMAKE_MINOR=8
set CMAKE_PATCH=12
set CMAKE_DIR=cmake-%CMAKE_MAJOR%.%CMAKE_MINOR%.%CMAKE_PATCH%.2-win32-x86
set CMAKE_ZIP=%CMAKE_DIR%.zip

::============================================================================
:: Define proxy settings
::============================================================================
:: set http_proxy=http://my.proxy.com:8080
:: set https_proxy=%http_proxy%
:: set ftp_proxy=%http_proxy%
:: set HTTP_PROXY=%http_proxy%
:: set HTTPS_PROXY=%http_proxy%
:: set FTP_PROXY=%http_proxy%

::============================================================================
:: Set paths for CFS++ environment.
::============================================================================
::----------------------------------------------------------------------------
:: Set path to base directory for CFS++ environment.
::----------------------------------------------------------------------------
set BASEDIR=%CD%

::----------------------------------------------------------------------------
:: Set path to cache directory for downloaded packages
::----------------------------------------------------------------------------
set PCKGDIR=%BASEDIR%\packages

::----------------------------------------------------------------------------
:: Set path to temporary directory
::----------------------------------------------------------------------------
set TMPDIR=%BASEDIR%\tmp

::============================================================================
:: Set minimal system path to minimize influences from other packages
::============================================================================
set PATH=%SystemRoot%\system32;%SystemRoot%;%SystemRoot%\System32\Wbem

IF NOT EXIST %BASEDIR% (
  md %BASEDIR%
)
cd /d %BASEDIR%

IF EXIST %TMPDIR% (
  rmdir /S /Q %TMPDIR%
)
md %TMPDIR%

IF NOT EXIST %PCKGDIR% (
  md %PCKGDIR%
)
cd /d %PCKGDIR%

::============================================================================
:: Try to download wget.exe with tools available on any Windows.
::============================================================================

IF EXIST wget.exe (
  goto :SetVars
)

::----------------------------------------------------------------------------
:: Create the temporary script file for downloading wget.exe from
:: ftp://ftp.lightspeedsystems.com/Brock/UnixTools/wget.exe
::----------------------------------------------------------------------------
echo Trying to download wget.exe using the Windows FTP client...
> script.ftp ECHO open ftp.lightspeedsystems.com
>>script.ftp ECHO anonymous
>>script.ftp ECHO hallo@world.com
>>script.ftp ECHO cd Brock
>>script.ftp ECHO cd UnixTools
>>script.ftp ECHO binary
>>script.ftp ECHO prompt n
>>script.ftp ECHO get wget.exe
>>script.ftp ECHO quit

:: Use the temporary script for unattended FTP
:: Note: depending on your OS version you may have to add a '-n' switch
FTP -v -s:script.ftp
IF EXIST wget.exe GOTO :SetVars

::----------------------------------------------------------------------------
:: Create the temporary script file for downloading wget.exe from
:: ftp://ftp.kfki.hu/pub/w2/wget.exe
::----------------------------------------------------------------------------
echo Trying to download wget.exe using the Windows FTP client...
>script.ftp ECHO open ftp.kfki.hu
>> script.ftp ECHO user anonymous
>>script.ftp ECHO hallo@world.com
>>script.ftp ECHO cd pub
>>script.ftp ECHO cd w2
>>script.ftp ECHO binary
>>script.ftp ECHO prompt n
>>script.ftp ECHO get wget.exe
>>script.ftp ECHO quit

:: Use the temporary script for unattended FTP
:: Note: depending on your OS version you may have to add a '-n' switch
FTP -v -s:script.ftp

:: For the paranoid: overwrite the temporary file before deleting it
TYPE NUL >script.ftp
DEL script.ftp

IF EXIST wget.exe GOTO :SetVars

:: If download via FTP was not successful, let's try using Internet Explorer
echo Trying to download wget.exe using Internet Explorer...
set MSB=%TMPDIR%\msgbox.vbs
>  %MSB% ECHO Set objArgs = WScript.Arguments
>> %MSB% ECHO URL = objArgs(0)
>> %MSB% ECHO MsgBox "Please select Save As in the download dialog for " _
>> %MSB% ECHO      ^& "wget.exe which will popup and save it to " _ 
>> %MSB% ECHO      ^& "%BASEDIR%\packages. Afterwards it may be necessary " _ 
>> %MSB% ECHO      ^& "to restart the bootstrap script.", _
>> %MSB% ECHO      64, "Information for downloading wget.exe" 
>> %MSB% ECHO Set oIE = CreateObject("InternetExplorer.Application")
>> %MSB% ECHO oIE.Navigate2 URL

set URL=ftp://ftp.lightspeedsystems.com/Brock/UnixTools/wget.exe
cscript /Nologo %MSB% %URL%
IF EXIST wget.exe GOTO :SetVars

set URL=ftp://ftp.kfki.hu/pub/w2/wget.exe
cscript /Nologo %MSB% %URL%
IF EXIST wget.exe GOTO :SetVars

echo Please restart the bootstrap script once you have downloaded wget.exe
echo to '%PCKGDIR%\wget.exe'.
exit /B
 
:SetVars
set WGET=%PCKGDIR%\wget.exe

::============================================================================
:: Download 7-Zip
:: set MD5SUM=9bd44a22bffe0e4e0b71b8b4cf3a80e2
::============================================================================
set SEVENZIP_MSI=7z920.msi
IF NOT EXIST %SEVENZIP_MSI% (
  %WGET% http://downloads.sourceforge.net/sevenzip/%SEVENZIP_MSI%
)
IF NOT EXIST %SEVENZIP_MSI% (
  %WGET% --passive-ftp %BOOTSTRAP_PCKG_MIRROR%/%SEVENZIP_MSI%
)
IF NOT EXIST %SEVENZIP_MSI% (
  echo Download of 7-Zip failed! Please put %SEVENZIP_MSI% into %PCKGDIR%
  echo manually.
  exit /B
)

::============================================================================
:: Download CMake
::============================================================================
IF NOT EXIST %CMAKE_ZIP% (
  %WGET% http://www.cmake.org/files/v%CMAKE_MAJOR%.%CMAKE_MINOR%/%CMAKE_ZIP%
)
IF NOT EXIST %CMAKE_ZIP% (
  %WGET% --passive-ftp %BOOTSTRAP_PCKG_MIRROR%/%CMAKE_ZIP%
)
IF NOT EXIST %CMAKE_ZIP% (
  echo Download of CMake failed! Please put %CMAKE_ZIP% into %PCKGDIR%
  echo manually.
  exit /B
)

cd /d %BASEDIR%

IF NOT EXIST 7-Zip (
  MSIEXEC /a %PCKGDIR%\%SEVENZIP_MSI% TARGETDIR=%TMPDIR% /qn
  move %TMPDIR%\Files\7-Zip %BASEDIR%\7-Zip
  cd /d %BASEDIR%
  rmdir /S /Q %TMPDIR%
  md %TMPDIR%
)
set PATH=%CD%\7-Zip;%PATH%

cd /d %BASEDIR%
IF NOT EXIST %CMAKE_DIR% (
  7z x %PCKGDIR%\%CMAKE_ZIP%
)
set PATH=%CD%\%CMAKE_DIR%\bin;%PATH%

IF NOT EXIST %BOOTSTRAP_CMAKE_SCRIPT% (
  %WGET% --passive-ftp %BOOTSTRAP_DOWNLOAD_URL%
)
IF NOT EXIST %BOOTSTRAP_CMAKE_SCRIPT% (
  set URL=%CFS_DEVEL_FTP_SERVER%/binaries/mingw_env/%BOOTSTRAP_CMAKE_SCRIPT%
  %WGET% --passive-ftp %URL%
)

::============================================================================
:: Hand over control to CMake bootstrap script.
::============================================================================
cmake -P %BOOTSTRAP_CMAKE_SCRIPT%
