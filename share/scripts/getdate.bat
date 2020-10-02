@ECHO off
rem This script has been aquired from:
rem http://stackoverflow.com/questions/203090/how-to-get-current-datetime-on-windows-command-line-in-a-suitable-format-for-us

SETLOCAL ENABLEEXTENSIONS
if "%date%A" LSS "A" (set toks=1-3) else (set toks=2-4)
for /f "tokens=2-4 delims=(-)" %%a in ('echo:^|date') do (
for /f "tokens=%toks% delims=.-/ " %%i in ('date/t') do (
    set '%%a'=%%i
    set '%%b'=%%j
    set '%%c'=%%k))

if "%dd%A" EQU "A" (set 'dd'=%'tt'%)
if "%yy%A" EQU "A" (set 'yy'=%'jj'%)

if %'yy'% LSS 100 set 'yy'=20%'yy'%
set Today=%'yy'%-%'mm'%-%'dd'% 
ENDLOCAL & SET v_year=%'yy'%& SET v_month=%'mm'%& SET v_day=%'dd'%

rem ECHO Today is Year: [%V_Year%] Month: [%V_Month%] Day: [%V_Day%]
echo %V_Day%/%V_Month%/%V_Year%

rem wmic os get LocalDateTime

