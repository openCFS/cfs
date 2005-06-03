#!/bin/tcsh
set ACTDATE=`date '+%B%_e, %G'`
sed "s/Last.*/Last modified: $ACTDATE/g" $1 > _temp.html
mv _temp.html $1
