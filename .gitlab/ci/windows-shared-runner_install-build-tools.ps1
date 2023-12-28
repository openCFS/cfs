# this should be run from the CFS root directory

# install Fortran compiler
.gitlab/ci/install_windows.bat https://registrationcenter-download.intel.com/akdlm/irc_nas/18976/w_HPCKit_p_2022.3.1.19755_offline.exe intel.oneapi.win.ifort-compiler
# remove unnecessary stuff
$LATEST_VERSION=Get-ChildItem -Path "$env:INTEL_INSTALL_DIR\compiler\" -Name | Select-String -NotMatch latest | %{$_.Line} | Sort-Object | Select-Object -Last 1

Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\compiler\lib\ia32_win" -Force  -Recurse -ErrorAction SilentlyContinue
Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\bin\intel64_ia32" -Force  -Recurse -ErrorAction SilentlyContinue
Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\lib\emu" -Force  -Recurse -ErrorAction SilentlyContinue
Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\lib\oclfpga" -Force  -Recurse -ErrorAction SilentlyContinue
Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\lib\ocloc" -Force  -Recurse -ErrorAction SilentlyContinue
Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\lib\x86" -Force  -Recurse -ErrorAction SilentlyContinue

# install MKL
.gitlab/ci/install_windows.bat https://registrationcenter-download.intel.com/akdlm/irc_nas/18969/w_BaseKit_p_2022.3.1.19798_offline.exe intel.oneapi.win.mkl.devel

# install git
wget.exe --quiet https://github.com/git-for-windows/git/releases/download/v2.43.0.windows.1/PortableGit-2.43.0-64-bit.7z.exe
./PortableGit-2.43.0-64-bit.7z.exe -o"${pwd}/cache/git" -y
Remove-Item PortableGit-2.43.0-64-bit.7z.exe -Force -Recurse -ErrorAction SilentlyContinue