@echo off

set EXITCODE=0

rem Get name of NACS/bin directory
rem set NACS_BIN_DIR=%NACS_ROOT_DIR%bin\

rem Check if we are in BINARY tree or in DIST tree.
if EXIST "%NACS_ROOT_DIR%\source" (
rem We are in BINARY tree
rem    echo Script has been started from BINARY tree.
    set BINARY_TREE=1
) else (
rem We are in DIST tree
rem    echo Script has been started from DIST tree.
    set BINARY_TREE=0
)

rem Get architecture and distribution
set WINVER_BAT="%NACS_ROOT_DIR%\share\scripts\winver.bat"

call %WINVER_BAT% -u >> nul 2>&1
if %EXITCODE% NEQ 0 goto end

set WINDOWS_ARCH_STR=%WINDOWS_PLATFORM%_%WINDOWS_ARCH%

rem Set lib path according to architecture
if "_%WINDOWS_ARCH%_" == "_I386_" (
   set LIB=lib
)
if "_%WINDOWS_ARCH%_" == "_X86_64_" (
   set LIB=lib64
)
if "_%WINDOWS_ARCH%_" == "_IA64_" (
   set LIB=lib
)

rem Set standard Windows (XP) system path
set PATH=%SystemRoot%\System32;%SystemRoot%;%SystemRoot%\System32\Wbem

rem Set (library) path for current architecture
set PATH="%NACS_ROOT_DIR%\%LIB%\%WINDOWS_ARCH_STR%";%PATH%

if defined NACS_SCRIPT_DEBUG (
    echo NACS_ROOT_DIR: "%NACS_ROOT_DIR%"
    echo NACS_BIN_DIR: "%NACS_BIN_DIR%"
    echo OS: Windows
    echo BINARY_TREE: %BINARY_TREE%
    echo WINDOWS_ARCH_STR: %WINDOWS_ARCH_STR%
    echo WINDOWS_PLATFORM: %WINDOWS_PLATFORM%
    echo WINDOWS_VERSION: %WINDOWS_VERSION%
    echo WINDOWS_ARCH: %WINDOWS_ARCH%
    echo PATH: %PATH%
)

:end
