#!/usr/bin/csh
if (($#argv != 3)&&($#argv != 1)) then
  echo
  echo "       Usage: $0 databasename "
  echo "                               to get list of calculations"
  echo 
  echo "              $0 databasename nodalresult resultidx"
  echo "                               to get nodal result values"
  echo
  echo "              $0 databasename elementresult resultidx"
  echo "                               to get element result values"
  echo
  echo "              $0 databasename elements resultidx"
  echo "                               to see which nodes belong to wat element"
  echo
  echo "              $0 databasename coordinates resultidx"
  echo "                               to get coordinates of nodes"
  echo
  exit
endif

set DBNAME = $argv[1]

if ($#argv == 1) then
  cat getresultidx.sql | mysql $DBNAME
  exit
endif

set DATASET = $argv[2]
set IDX	   = $argv[3]

cat get${DATASET}.sql | sed s/RESULTIDX/${IDX}/g | mysql $DBNAME
