# this should be run from the CFS root directory

# install ccache (https://ccache.dev/download.html)
wget.exe --quiet https://github.com/ccache/ccache/releases/download/v4.10.2/ccache-4.10.2-windows-x86_64.zip
Expand-Archive -Path ccache-4.10.2-windows-x86_64.zip -DestinationPath ccache-extract
New-Item -ItemType directory -Path cache/ccache
Get-Childitem -Path ccache-extract -Include "ccache.exe" -File -Recurse | Copy-Item  -Destination cache/ccache
# cleanup
Remove-Item ccache-4.10.2-windows-x86_64.zip -Force -Recurse -ErrorAction SilentlyContinue
Remove-Item ccache-extract -Force -Recurse -ErrorAction SilentlyContinue

