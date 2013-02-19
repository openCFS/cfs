#!/bin/sh
# init script for /usr/local/passerver/passerverd
#
# Author: strieben 
#
# /etc/init.d/gidpasserver
#   and its symbolic  link
# /usr/sbin/gidpasserver
#
### BEGIN INIT INFO
# Provides:          gidpasserver
# Required-Start:    $network 
# Required-Stop:     $network
# Default-Start:     3 5
# Default-Stop:
# Description:       GiD password server
# Short-Description: gidpasserver
### END INIT INFO
 
PIDFILE=/var/run/passerverd.pid
case "$1" in
  'start')
        echo -n  "Starting passerverd..."
        /usr/local/passerver/passerverd -port 7073 &
        echo $! > $PIDFILE
        echo "done."
        ;;
  'stop')
       kill `cat $PIDFILE`
       rm -f $PIDFILE
        ;;
  *)
        echo "usage: $0 {start|stop}"
        ;;
esac
exit 0
