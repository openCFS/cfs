#!/bin/bash
# this script executes ps aux for the given id periodically
# 
while [ 1 = 1 ]
do 
  ps aux | grep "$1" | grep -v grep | grep -v monitor_memory 
  sleep 1
done
