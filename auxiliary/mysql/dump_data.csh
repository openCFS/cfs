#!/bin/csh
if ($#argv != 2) then 
  echo "     Usage: $0 databasename username"
  exit
endif

set dbname	= $argv[1]
set user	= $argv[2]

mysqldump $dbname > ${dbname}_data.sql
