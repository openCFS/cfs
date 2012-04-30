# Directory on server
SCRIPTS_DIR='ftp://lse17.e-technik.uni-erlangen.de:40065/scripts'

# Download script
wget $SCRIPTS_DIR/bootstrap_devel_machine.sh

# Make script executable
chmod 755 bootstrap_devel_machine.sh

# Change into super user mode
sudo sh # or su

# Execute bootstrapper script
./bootstrap_devel_machine.sh
