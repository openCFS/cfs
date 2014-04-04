# Directory on server
SCRIPTS_DIR='@CFS_DS_FTP@/scripts'

# Download script
wget $SCRIPTS_DIR/bootstrap_devel_machine.sh

# Make script executable
chmod 755 bootstrap_devel_machine.sh

# Change into super user mode
sudo sh # or su

# Execute bootstrapper script
./bootstrap_devel_machine.sh
