@echo off

rem Batch script based on http://malektips.com/xp_dos_0025.html
rem For a complete list of Windows versions google for 'Operating System Version (Windows)'
rem or for the 'GetVersionEx' Windows API function
rem cf. http://msdn.microsoft.com/en-us/library/windows/desktop/ms724832(v=vs.85).aspx

set EXITCODE=0
set ERRORLEVEL=

set PATH_SAVE=%PATH%

set PATH=%SystemRoot%\System32;%SystemRoot%;%SystemRoot%\System32\Wbem

ver | find.exe "2003" > nul
if %ERRORLEVEL% EQU 0 goto :ver_2003

ver | find.exe "XP" > nul
if %ERRORLEVEL% EQU 0 goto :ver_xp

rem windows xp x64 identifies itself only by version
ver | find.exe "Version 5.2" > nul
if %ERRORLEVEL% EQU 0 goto :ver_xp

rem Windows Server 2008 version 6.0 or 6.1
ver | find.exe "2008" > nul
if %ERRORLEVEL% EQU 0 goto :ver_2008

rem windows vista identifies itself only by version
ver | find.exe "Version 6.0" > nul
if %ERRORLEVEL% EQU 0 goto :ver_vista

rem windows 7/windows server 2008 r2
ver | find.exe "Version 6.1" > nul
if %ERRORLEVEL% EQU 0 goto :ver_7

rem Windows Server 2012 version 6.2
ver | find.exe "2012" > nul
if %ERRORLEVEL% EQU 0 goto :ver_2012

rem windows 8/windows server 2012
ver | find.exe "Version 6.2" > nul
if %ERRORLEVEL% EQU 0 goto :ver_8

rem Windows 2000 version 5.0
ver | find.exe "2000" > nul
if %ERRORLEVEL% EQU 0 goto :ver_2000

ver | find.exe "NT" > nul
if %ERRORLEVEL% EQU 0 goto :ver_nt

echo Machine undetermined.
set EXITCODE=1
goto :end

:ver_8
rem Run Windows 8-specific commands here.
set WINDOWS_PLATFORM=WIN8
goto :arch

:ver_2012
rem Run Windows Server 2012-specific commands here.
set WINDOWS_PLATFORM=WIN2012
goto :arch

:ver_7
rem Run Windows 7-specific commands here.
set WINDOWS_PLATFORM=WIN7
goto :arch

:ver_vista
rem Run Windows 7-specific commands here.
set WINDOWS_PLATFORM=WINVISTA
goto :arch

:ver_2008
rem Run Windows Server 2008-specific commands here.
set WINDOWS_PLATFORM=WIN2008
goto :arch

:ver_2003
rem Run Windows 2003-specific commands here.
set WINDOWS_PLATFORM=WIN2003
goto :arch

:ver_xp
rem Run Windows XP-specific commands here.
set WINDOWS_PLATFORM=WINXP
goto :arch

:ver_2000
rem Run Windows 2000-specific commands here.
set WINDOWS_PLATFORM=WIN2000
goto :arch

:ver_nt
rem Run Windows NT-specific commands here.
set WINDOWS_PLATFORM=WINNT
goto :arch

if [%OS%] == []
set WINDOWS_PLATFORM="WIN9X"
echo "Windows platform %WINDOWS_PLATFORM% not supported"
set EXITCODE=1
goto :end

:arch
set REGKEY="HKLM\SYSTEM\CurrentControlSet\Control\Session Manager\Environment"

reg query %REGKEY% /v PROCESSOR_ARCHITECTURE | find.exe "x86" > nul
if %ERRORLEVEL% EQU 0 goto :arch_x86

reg query %REGKEY% /v PROCESSOR_ARCHITECTURE | find.exe "AMD64" > nul
if %ERRORLEVEL% EQU 0 goto :arch_x64

reg query %REGKEY% /v PROCESSOR_ARCHITECTURE | find.exe "EM64T" > nul
if %ERRORLEVEL% EQU 0 goto :arch_x64

reg query %REGKEY% /v PROCESSOR_ARCHITECTURE | find.exe "IA64" > nul
if %ERRORLEVEL% EQU 0 goto :arch_ia64

:arch_x86
set WINDOWS_SYS_ARCH=I386
goto :win_version

:arch_x64
set WINDOWS_SYS_ARCH=X86_64
goto :win_version

:arch_ia64
set WINDOWS_SYS_ARCH=IA64
goto :win_version

:win_version
if "%BUILD_ARCH%" == "" set WINDOWS_ARCH=%WINDOWS_SYS_ARCH%
if [%BUILD_ARCH%] == [I386] set WINDOWS_ARCH=I386
if [%BUILD_ARCH%] == [AMD64] set WINDOWS_ARCH=X86_64
if [%BUILD_ARCH%] == [EMT64] set WINDOWS_ARCH=X86_64
if [%BUILD_ARCH%] == [IA64] set WINDOWS_ARCH=IA64

FOR /f "tokens=2* delims=[" %%i in ('ver') do set WINDOWS_VERSION=%%i
FOR /f "tokens=2 delims= " %%i in ("%WINDOWS_VERSION%") do set WINDOWS_VERSION=%%i
FOR /f "tokens=1 delims=]" %%i in ("%WINDOWS_VERSION%") do set WINDOWS_VERSION=%%i

if "%1"=="-a" echo %WINDOWS_PLATFORM% %WINDOWS_VERSION% %WINDOWS_ARCH% %WINDOWS_SYS_ARCH%

if "%1"=="-u" echo %WINDOWS_PLATFORM%_%WINDOWS_VERSION%_%WINDOWS_ARCH%_%WINDOWS_SYS_ARCH%

if [%1]==[-c] echo SET(OS "WINDOWS");SET(DIST "%WINDOWS_PLATFORM%");SET(DIST_FAMILY "WINDOWS");SET(REV "%WINDOWS_VERSION%");SET(ARCH "%WINDOWS_ARCH%");SET(SUBARCH "%WINDOWS_SYS_ARCH%")

:end

set PATH=%PATH_SAVE%
set PATH_SAVE=

