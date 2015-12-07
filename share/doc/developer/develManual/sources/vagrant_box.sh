CFSDSHOSTNAME=@CFS_DS_HOSTNAME@
BASE_BOX=oracle6.box

mkdir oracle6 && cd oracle6    # Create temp directory
echo "get boxes/$BASE_BOX" | sftp ftpUser@$CFSDSHOSTNAME
vagrant init oracle6 $BASE_BOX # Create Vagrantfile
vagrant up                     # Import, start & provision
vagrant ssh                    # Login to VBox

# Inside the Box, share /vagrant points to oracle6 dir
cd /vagrant && ls
exit

vagrant halt                   # Shutdown VBox
vagrant destroy                # Destroy VBox
