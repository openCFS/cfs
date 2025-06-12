#!/bin/bash

# use this script to keep you CFS installation up to date. It will
# 1. download the latests Linux CFS build for the 'master' branch
# 2. extract the tar.gz file, and save it under the name CFS-<sha>-Linux.tar.gz
# 3. set a symlink 'latest' to the extracted directory

if [[ $# -ne 2 ]];
then
    echo "usage: ./install-package.sh version output_dir"
    exit 1
else 
    PACKAGE_VERSION=$1
    OUTPUT_DIR=$2
fi

PROJECT_ID="12930334"
PACKAGE_NAME="openCFS"
file="CFS-$PACKAGE_VERSION-Linux.tar.gz"

URL="https://gitlab.com/api/v4/projects/$PROJECT_ID/packages/generic/$PACKAGE_NAME/$PACKAGE_VERSION/$file"

if [ -f $OUTPUT_DIR/$file ]; then
  echo "$OUTPUT_DIR/$file already exists!"
  exit 1
else
  echo "downloading $URL to $OUTPUT_DIR/$file"
  wget --output-document=$OUTPUT_DIR/$file $URL
fi

# extract and write a log file
echo "extracting $OUTPUT_DIR/$file"
tar -xzvf $OUTPUT_DIR/$file -C $OUTPUT_DIR > $OUTPUT_DIR/last-install.log

# get the extracted archive from the log file
install_version_dir=$( grep "/bin/cfs$" $OUTPUT_DIR/last-install.log | sed 's|/bin/cfs||' )
echo "installed to $install_version_dir"

# rename to keep the archive under a unique name, necessary from master
if [! -f $OUTPUT_DIR/$install_version_dir.tar.gz]; then
    echo "probably master: renaming $OUTPUT_DIR/$file -> $OUTPUT_DIR/$install_version_dir.tar.gz"
    mv $OUTPUT_DIR/$file $OUTPUT_DIR/$install_version_dir.tar.gz
fi

# set a symlink for the latests version
# symlink to new version
ln -s $install_version_dir $OUTPUT_DIR/new-cfs
# overwrite the 'latest' symlink 
mv -Tf $OUTPUT_DIR/new-cfs $OUTPUT_DIR/latest