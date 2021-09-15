@echo off

set EXITCODE=0

rem Get name of CFS/bin directory
rem set CFS_BIN_DIR=%CFS_ROOT_DIR%bin\

rem Check if we are in BINARY tree or in DIST tree.
if EXIST "%CFS_ROOT_DIR%\source" (
rem We are in BINARY tree
rem    echo Script has been started from BINARY tree.
    set BINARY_TREE=1
) else (
rem We are in DIST tree
rem    echo Script has been started from DIST tree.
    set BINARY_TREE=0
)

rem Get architecture and distribution
set WINVER_BAT="%CFS_ROOT_DIR%\share\scripts\winver.bat"

call %WINVER_BAT% -u >> nul 2>&1
if %EXITCODE% NEQ 0 goto end


rem we do not really have different architecture
rem echo "_%WINDOWS_ARCH%_"
set LIB=lib64

rem Set standard Windows (XP) system path
set PATH=%SystemRoot%\System32;%SystemRoot%;%SystemRoot%\System32\Wbem

rem Set (library) path for current architecture
set PATH="%CFS_ROOT_DIR%\%LIB%";%PATH%

if defined CFS_SCRIPT_DEBUG (
    echo CFS_ROOT_DIR: "%CFS_ROOT_DIR%"
    echo CFS_BIN_DIR: "%CFS_BIN_DIR%"
    echo OS: Windows
    echo BINARY_TREE: %BINARY_TREE%
    echo WINDOWS_ARCH_STR: %WINDOWS_ARCH_STR%
    echo WINDOWS_PLATFORM: %WINDOWS_PLATFORM%
    echo WINDOWS_VERSION: %WINDOWS_VERSION%
    echo WINDOWS_ARCH: %WINDOWS_ARCH%
    echo PATH: %PATH%
)

:end
