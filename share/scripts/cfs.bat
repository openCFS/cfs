@echo off
rem This is the start script for CFS. It should be placed in
rem CFS_DIST/bin under the name 'nacs.bat'

setlocal

set EXITCODE=0
set ERRORLEVEL=

set WORKDIR=%CD%

rem for %%i in (%0) do set CFS_BIN_DIR=%CFS_BIN_DIR%%%~di%%~pi
cd /d %~dp0 

if %ERRORLEVEL% NEQ 0 (
     echo Can not chdir to script base dir!
     set EXITCODE=1
     goto end
)

set BASEDIR=%CD%

rem Check if CFS_ROOT_DIR is defined
if not defined CFS_ROOT_DIR (
   goto :nacs_root_dir_undefined
)

rem Replace double quotes in CFS_ROOT_DIR with empty string
set CFS_ROOT_DIR=%CFS_ROOT_DIR:"=%

rem Try to change into CFS_ROOT_DIR
cd /d "%CFS_ROOT_DIR%" >> nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    if NOT EXIST "%CFS_ROOT_DIR%" (
        echo CFS_ROOT_DIR does not exist.
        set EXITCODE=1
        goto end
    )

    echo Error while trying to cd into CFS_ROOT_DIR.
    set EXITCODE=1
    goto end
)

set CFS_ROOT_DIR=%CD%

rem Check if bin directory exists
if NOT EXIST bin (
   echo CFS_ROOT_DIR\bin does not exist.
   set EXITCODE=1
   goto end
)

cd bin
set CFS_BIN_DIR=%CD%

goto :call_common_bat


:nacs_root_dir_undefined
rem If CFS_ROOT_DIR is undefined we just take the base dir
rem of the current script as CFS_BIN_DIR
set CFS_BIN_DIR=%BASEDIR%
cd ..
set CFS_ROOT_DIR=%CD%

if defined CFS_SCRIPT_DEBUG (
   echo WARNING: CFS_ROOT_DIR has not been set!
   echo Using '%CFS_ROOT_DIR%' as CFS_ROOT_DIR.
)

rem Call script with common definitions
:call_common_bat
rem echo CFS_ROOT_DIR: "%CFS_ROOT_DIR%"
rem echo CFS_BIN_DIR: "%CFS_BIN_DIR%"
call "%CFS_BIN_DIR%"\common.bat

if %EXITCODE% NEQ 0  (
    echo Error while trying to call CFS_BIN_DIR\common.bat
    goto end
)

rem Set XML schema root path
set CFS_SCHEMA_ROOT="%CFS_ROOT_DIR%\share\xml"

rem Set path to CFS executable
set CFS_EXE="%CFS_BIN_DIR%\%WINDOWS_ARCH_STR%\cfsbin.exe"

set PATH=%CFS_ROOT_DIR%\%LIB%\%WINDOWS_ARCH_STR%;%CFS_ROOT_DIR%\%LIB%\%WINDOWS_ARCH_STR%\v110\winx64;%CFS_ROOT_DIR%\bin\%WINDOWS_ARCH_STR%;%PATH%

if defined CFS_SCRIPT_DEBUG (
    echo CFS_SCHEMA_ROOT: %CFS_SCHEMA_ROOT%
    echo CFS_EXE: %CFS_EXE%
    echo PATH: %PATH%
)

if NOT EXIST %CFS_EXE% (
    echo Cannot find CFS executable for architecture %WINDOWS_ARCH_STR%:
    echo   %CFS_EXE%
    set EXITCODE=1
    goto end
)

cd /d "%WORKDIR%"
if %ERRORLEVEL% NEQ 0 (
    set WORKDIR=%TEMPDIR%
    cd /d "%WORKDIR%"
    if %ERRORLEVEL% NEQ 0 (
        set WORKDIR=%TMPDIR%
        cd /d "%WORKDIR%"
        if %ERRORLEVEL NEQ 0 (
            echo Can not chdir back to "%WORKDIR%"
            echo and could not change to temp directory!
            set EXITCODE=1
            goto end
        )
    )
)

rem Remove all quotes from PATH
set PATH=%PATH:"=%

rem Run executable with arguments
%CFS_EXE% %*

:end

rem Properly exit from script.
if %EXITCODE% NEQ 0 (
   if defined WINDOWS_PLATFORM (
       if "_%WINDOWS_PLATFORM%_" == "_WIN9X_" (
           endlocal
           exit /B 1
       )
       if "_%WINDOWS_PLATFORM%_" == "_WINNT_" (
           endlocal
           color 00
       )
       endlocal
       exit /B %EXITCODE%
   ) else (
       endlocal
       exit /B 1
   )
)
