#!/bin/sh

env -i SHELL=/bin/bash PATH=/sbin:/bin:/usr/sbin:/usr/bin HOME=/home/oracle /opt/pckg/cmake-2.8.10.2-Linux-i386/bin/cmake -P /home/oracle/Documents/dev/NIGHTLY/CFS_TRUNK_NIGHTLY/ctest_scripts/nightly_test.cmake > /home/oracle/Documents/dev/nightly_test.log 2>&1
