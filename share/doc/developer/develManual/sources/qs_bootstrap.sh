# Rename scripts
mv distro.txt distro.sh
mv bootstrap_devel_machine.txt bootstrap_devel_machine.sh

# Make scripts executable
chown 755 distro.sh bootstrap_devel_machine.sh

# Change into super user mode
sudo bash # or su

# Execute bootstrapper script
./bootstrap_devel_machine.sh
