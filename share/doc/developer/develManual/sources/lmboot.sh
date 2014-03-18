#!/bin/sh
cd /var/log/ansys
/usr/bin/nohup /opt/pckg/ansys/shared_files/licensing/linia32/boot_ansflex
/bin/sleep 5
/opt/pckg/ansys/shared_files/licensing/linia32/lmutil lmstat \
    -c /home/users/strieben/lic.dat -a
