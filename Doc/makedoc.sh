#!/bin/tcsh

#
# Generate CFS++ documentation using doxygen
#
echo " "
echo " Generating Doxygen Documentation for CFS++"
echo " "
doxygen DoxygenConfig/doxy-config
echo " "
echo " See `pwd`/html/index.html for doc"
echo " "
