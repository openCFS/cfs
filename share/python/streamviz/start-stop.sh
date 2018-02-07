#!/bin/bash

if [ "$#" -eq "0" ]; then
  echo "start-stop.sh start"
  echo "  --> start the listener"
  echo ""
  echo "start-stop.sh stop"
  echo "  --> stops the listener"
  echo ""
  echo "start-stop.sh clean"
  echo "  --> deletes the pid file"
 exit
fi

case $1 in
 start)
   if [ -e pid.lck ];then
     echo "server already running, check pid file"; exit 2;
   fi
   nohup ./streamviz.py </dev/null >./log.txt 2>&1 &
   echo "$!" > ./pid.lck
   sleep 1
   echo "started server with pid $!";
   sleep 1
  ;;
 stop)
   if [ -e pid.lck ];then
    kill `cat pid.lck`
    rm pid.lck
   else
    echo "no pid file, server might not be running"
   fi
  ;;
 clean)
   if [ -e pid.lck ];then
    kill `cat pid.lck`
    rm pid.lck
   fi
  ;;
esac
