#!/bin/tcsh

#
# Generate CFS++ documentation using doxygen
#
echo " "
echo " Generating Doxygen Documentation for CFS++"
echo " "
# Change date in footer.html
DoxygenConfig/changeDate.sh DoxygenConfig/footerTemplate.html \
DoxygenConfig/footer.html 
doxygen DoxygenConfig/doxy-config
echo " "
echo " See `pwd`/html/index.html for doc"
echo " "
