# this should be run from the CFS root directory

# install ccache (https://ccache.dev/download.html)
wget.exe --quiet https://github.com/ccache/ccache/releases/download/v4.10.2/ccache-4.10.2-windows-x86_64.zip
Expand-Archive -LiteralPath "ccache-4.10.2-windows-x86_64.zip" -DestinationPath "ccache-install"
Remove-Item "ccache-4.10.2-windows-x86_64.zip" -Force -Recurse -ErrorAction SilentlyContinue
New-Item  -ItemType "directory" -Path "$env:CI_PROJECT_DIR/cache/ccache"
Copy-Item -Path "ccache-install" -Filter "ccache.exe" -Destination "$env:CI_PROJECT_DIR/cache/ccache" –Recurse
#Get-ChildItem -Path ".\*.exe" -Recurse | Move-Item -Destination "$env:CI_PROJECT_DIR/cache/ccache"
cmake -E rm -rf -- "ccache-install"

# install Fortran compiler
.gitlab/ci/install_windows.bat https://registrationcenter-download.intel.com/akdlm/IRC_NAS/7a6db8a1-a8b9-4043-8e8e-ca54b56c34e4/w_HPCKit_p_2024.0.1.35_offline.exe intel.oneapi.win.ifort-compiler
# remove unnecessary stuff
$LATEST_VERSION=Get-ChildItem -Path "$env:INTEL_INSTALL_DIR\compiler\" -Name | Select-String -NotMatch latest | %{$_.Line} | Sort-Object | Select-Object -Last 1

Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\compiler\lib\ia32_win" -Force  -Recurse -ErrorAction SilentlyContinue
Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\bin\intel64_ia32" -Force  -Recurse -ErrorAction SilentlyContinue
Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\lib\emu" -Force  -Recurse -ErrorAction SilentlyContinue
Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\lib\oclfpga" -Force  -Recurse -ErrorAction SilentlyContinue
Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\lib\ocloc" -Force  -Recurse -ErrorAction SilentlyContinue
Remove-Item "$env:INTEL_INSTALL_DIR\compiler\$LATEST_VERSION\windows\lib\x86" -Force  -Recurse -ErrorAction SilentlyContinue

# install MKL
.gitlab/ci/install_windows.bat https://registrationcenter-download.intel.com/akdlm/IRC_NAS/5cb30fb9-21e9-47e8-82da-a91e00191670/w_BaseKit_p_2024.0.1.45_offline.exe intel.oneapi.win.mkl.devel

# install git
wget.exe --quiet https://github.com/git-for-windows/git/releases/download/v2.43.0.windows.1/PortableGit-2.43.0-64-bit.7z.exe
./PortableGit-2.43.0-64-bit.7z.exe -o "$env:CI_PROJECT_DIR/cache/git" -y
Remove-Item PortableGit-2.43.0-64-bit.7z.exe -Force -Recurse -ErrorAction SilentlyContinue