echo "Configuring and building CFS in CMD since setting the environment does not work on PowerShell"

:::: set environment

:: first do MSVC and choose correct architecture (x86-32 does not work for ifort)
:: see https://learn.microsoft.com/en-us/cpp/build/building-on-the-command-line
echo "Initialize Visual Studio for x64 (has to be done before oneAPI)"
@call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" amd64

:: then do oneAPI which adapts to MSVC environment automatically
echo "Set oneAPI environment"
@call "%INTEL_INSTALL_DIR%\setvars.bat"

:::: output some debug information
echo "PATH=%PATH%"
echo "CC=%CC%"
echo "CXX=%CXX%"
echo "FC=%FC%"
%FC% /QV
echo "CFS_CONFIG_OPTIONS=%CFS_CONFIG_OPTIONS%"
echo "CDASH_TRACK=%CDASH_TRACK%"

:::: configure (for the first time)
cmake -G "NMake Makefiles" -DBUILDNAME="%CI_COMMIT_REF_SLUG% / %BUILD_NAME%" %CFS_CONFIG_OPTIONS% ..

:::: now do the build (defined by an environment variable)
echo "%BUILD_COMMAND%"
%BUILD_COMMAND%