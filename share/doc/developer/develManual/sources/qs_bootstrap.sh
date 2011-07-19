# Rename script
mv bootstrap_devel_machine.txt bootstrap_devel_machine.sh

# Make script executable
chown 755 bootstrap_devel_machine.sh

# Change into super user mode
sudo bash # or su

# Execute bootstrapper script
./bootstrap_devel_machine.sh
