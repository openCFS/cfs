@echo off
rem This is the start script for NACS. It should be placed in
rem NACS_DIST/bin under the name 'nacs.bat'

setlocal

set EXITCODE=0
set ERRORLEVEL=

set WORKDIR=%CD%

rem for %%i in (%0) do set NACS_BIN_DIR=%NACS_BIN_DIR%%%~di%%~pi
cd /d %~dp0 

if %ERRORLEVEL% NEQ 0 (
     echo Can not chdir to script base dir!
     set EXITCODE=1
     goto end
)

set BASEDIR=%CD%

rem Check if NACS_ROOT_DIR is defined
if not defined NACS_ROOT_DIR (
   goto :nacs_root_dir_undefined
)

rem Replace double quotes in NACS_ROOT_DIR with empty string
set NACS_ROOT_DIR=%NACS_ROOT_DIR:"=%

rem Try to change into NACS_ROOT_DIR
cd /d "%NACS_ROOT_DIR%" >> nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    if NOT EXIST "%NACS_ROOT_DIR%" (
        echo NACS_ROOT_DIR does not exist.
        set EXITCODE=1
        goto end
    )

    echo Error while trying to change into NACS_ROOT_DIR.
    set EXITCODE=1
    goto end
)

set NACS_ROOT_DIR=%CD%

rem Check if bin directory exists
if NOT EXIST bin (
   echo NACS_ROOT_DIR\bin does not exist.
   set EXITCODE=1
   goto end
)

cd bin
set NACS_BIN_DIR=%CD%

goto :call_common_bat


:nacs_root_dir_undefined
rem If NACS_ROOT_DIR is undefined we just take the base dir
rem of the current script as NACS_BIN_DIR
set NACS_BIN_DIR=%BASEDIR%
cd ..
set NACS_ROOT_DIR=%CD%

echo WARNING: NACS_ROOT_DIR has not been set!
echo Using '%NACS_ROOT_DIR%' as NACS_ROOT_DIR.

rem Call script with common definitions
:call_common_bat
rem echo NACS_ROOT_DIR: "%NACS_ROOT_DIR%"
rem echo NACS_BIN_DIR:  "%NACS_BIN_DIR%"
call "%NACS_BIN_DIR%"\common.bat

if %EXITCODE% NEQ 0  (
    echo Error while trying to call NACS_BIN_DIR\common.bat
    goto end
)

rem Set XML schema root path
set NACS_SCHEMA_ROOT="%NACS_ROOT_DIR%\share\xml"

rem Set path to NACS executable
set NACSTOOL_EXE="%NACS_BIN_DIR%\%WINDOWS_ARCH_STR%\nacstoolbin.exe"

if defined NACS_SCRIPT_DEBUG (
    echo NACS_SCHEMA_ROOT: %NACS_SCHEMA_ROOT%
    echo NACS_EXE: %NACSTOOL_EXE%
    echo PATH: %PATH%
)

if NOT EXIST %NACSTOOL_EXE% (
    echo Cannot find NACS executable for architecture %WINDOWS_ARCH_STR%:
    echo   %NACSTOOL_EXE%
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
%NACSTOOL_EXE% %*

:end

rem Properly exit from script.
if %EXITCODE% NEQ 0 (
   if defined WINDOWS_PLATFORM (
       if "_%WINDOWS_PLATFORM%_" == "_WIN9X_" (
           endlocal
           exit 1
       )
       if "_%WINDOWS_PLATFORM%_" == "_WINNT_" (
           endlocal
           color 00
       )
       endlocal
       exit /B %EXITCODE%
   ) else (
       endlocal
       exit 1
   )
)
