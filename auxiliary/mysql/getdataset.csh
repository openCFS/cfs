#!/usr/bin/csh
if (($#argv != 3)&&($#argv != 1)) then
  echo
  echo "       Usage: $0 databasename "
  echo "                               to get list of calculations"
  echo
  echo "              $0 databasename dataset resultidx"
  echo "                               to get values for given dataset"
  exit
endif

set DBNAME = $argv[1]

if ($#argv == 1) then
  cat getresultidx.sql | mysql $DBNAME
  exit
endif

set DATASET = $argv[2]
set IDX	   = $argv[3]

cat getdataset${DATASET}.sql | sed s/RESULTIDX/${IDX}/g | mysql $DBNAME
