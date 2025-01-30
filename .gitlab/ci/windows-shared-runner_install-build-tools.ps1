# this should be run from the CFS root directory

# install ccache (https://ccache.dev/download.html)
wget.exe --quiet https://github.com/ccache/ccache/releases/download/v4.10.2/ccache-4.10.2-windows-x86_64.zip
Expand-Archive -Path ccache-4.10.2-windows-x86_64.zip -DestinationPath ccache-extract
New-Item -ItemType directory -Path cache/ccache
Get-Childitem -Path ccache-extract -Include "ccache.exe" -File -Recurse | Copy-Item  -Destination cache/ccache
# cleanup
Remove-Item ccache-4.10.2-windows-x86_64.zip -Force -Recurse -ErrorAction SilentlyContinue
Remove-Item ccache-extract -Force -Recurse -ErrorAction SilentlyContinue

# we install ifx and mkl from the same source
# we get and install Intel oneAPI based on https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html  
wget.exe --quiet https://registrationcenter-download.intel.com/akdlm/IRC_NAS/a37c30c3-a846-4371-a85d-603e9a9eb94c/intel-oneapi-hpc-toolkit-2025.0.1.48_offline.exe -O oneapi_offline.exe

# install in two steps. We are Administrator, see https://learn.microsoft.com/de-de/archive/blogs/virtual_pc_guy/a-self-elevating-powershell-script
# inspired by https://github.com/oneapi-src/oneapi-ci/blob/master/scripts/install_windows.bat
# note that $CI_PROJECT_DIR is the current dir C:\GitLab-Runner\builds\openCFS\cfs
Start-Process -FilePath "./oneapi_offline.exe" -ArgumentList "-s -x -f oneapi_extracted --log extract.log" -NoNewWindow -Wait
Start-Process -FilePath "./oneapi_extracted/bootstrapper.exe" -ArgumentList "-s --action install --install-dir=$CI_PROJECT_DIR\cache\oneAPI --log-dir=$CI_PROJECT_DIR\oneapilog  --components=intel.oneapi.win.ifort-compiler:intel.oneapi.win.mkl.devel --eula=accept -p=NEED_VS2019_INTEGRATION=0 -p=NEED_VS2022_INTEGRATION=0 " -NoNewWindow -Wait

# clean up temporary files
Remove-Item oneapi_offline.exe -Force 
Remove-Item oneapi_extracted -Force -Recurse 
Remove-Item extract.log -Force
Remove-Item oneapilog -Force -Recurse 

# remove unnecessary stuff. From originally 9 GB keep 3 GB
Remove-Item cache/oneAPI/ocloc -Force  -Recurse
Remove-Item cache/oneAPI/debugger -Force  -Recurse
Remove-Item cache/oneAPI/tcm -Force  -Recurse
Remove-Item cache/oneAPI/umf -Force  -Recurse
Remove-Item cache/oneAPI/mpi -Force  -Recurse
Remove-Item cache/oneAPI/compiler_ide -Force  -Recurse
Remove-Item cache/oneAPI/tbb -Force  -Recurse 
Remove-Item cache/oneAPI/pti -Force  -Recurse 
# sycl is actually the biggest part (about 5 GB)
Remove-Item cache/oneAPI/mkl/latest/lib/*sycl* -Force  -Recurse
Remove-Item cache/oneAPI/mkl/latest/bin/*sycl* -Force  -Recurse

# install git
wget.exe --quiet https://github.com/git-for-windows/git/releases/download/v2.43.0.windows.1/PortableGit-2.43.0-64-bit.7z.exe
./PortableGit-2.43.0-64-bit.7z.exe -o"$env:CI_PROJECT_DIR/cache/git" -y
Remove-Item PortableGit-2.43.0-64-bit.7z.exe -Force -Recurse -ErrorAction SilentlyContinue