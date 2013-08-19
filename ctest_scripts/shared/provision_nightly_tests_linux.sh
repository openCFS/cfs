#!/bin/bash

USER=$1
BLOCKING=$2

if [ ! "${BLOCKING}" = "blocking" ]; then
  BLOCKING="nonblocking"
fi

echo "Starting nightly test on $HOSTNAME..."
# fn=$(mktemp)
#tempfile=$TMP/nightly_test_$(od -N4 -tu /dev/random | awk 'NR==1 {print $2} {}')
TEMPFILE=$(mktemp -t nightly_test_XXXX -u)
echo "TEMPFILE $TEMPFILE"
echo "USER $USER"
echo "BLOCKING $BLOCKING"
# PWDLINE=$(fgrep $USER /etc/passwd | grep nightly_test_user)
PWDLINE=$(getent passwd $USER)
HOME=$(echo $PWDLINE | cut -d':' -f6)
echo "HOME $HOME"
CTEST="/opt/pckg/cmake-2.8.10.2-Linux-i386/bin/ctest"
NIGHTLY_DIR="$HOME/Documents/dev/NIGHTLY/CFS_FESPACE_NIGHTLY"
LOG_FILE="/vagrant/logs/nightly_test.log"
echo "LOG_FILE $LOG_FILE"
eval $(sh $NIGHTLY_DIR/share/scripts/distro.sh -s)
echo "DIST $DIST"
if [ "$DIST" = "UBUNTU" ]; then
  SHELL=$(which dash)
else
  SHELL=$(echo $PWDLINE | cut -d':' -f7)
fi
echo "SHELL $SHELL"

# Make sure all VBoxes have the time zone of Vienna
# cf. http://www.cyberciti.biz/faq/howto-linux-unix-change-setup-timezone-tz-variable/
TZ='Europe/Vienna'
echo "TZ $TZ"

# Generate shell script which actually starts the nightly test CMake script.
cat <<EOF > "${TEMPFILE}_2"
#!/bin/sh
mkdir -p /vagrant/logs
$CTEST -S $NIGHTLY_DIR/ctest_scripts/nightly_test.cmake > "$LOG_FILE" 2>&1
if [ "$BLOCKING" = "blocking" ]; then
  sudo /sbin/shutdown -h now
fi
EOF

# Generate a shell script which sets a crontab like environment and calls the
# nightly test shell script.
cat <<EOF > $TEMPFILE
#!/bin/sh
env -i PATH=/bin:/sbin:/usr/bin:/usr/sbin HOME=$HOME SHELL=$SHELL TZ=$TZ $SHELL "${TEMPFILE}_2"
EOF

if [ "$BLOCKING" = "blocking" ]; then
  sudo -u $USER $SHELL $TEMPFILE
else
  # Start a screen session for the automatic shutdown in 8 hours.
  screen -dmS  automatic_shutdown /sbin/shutdown -h -P 480

  # Start a screen session for starting the nightly tests. We use screen, because it immediately
  # returns and the host machine can go on starting the other VBoxes.
  sudo -u $USER screen -dmS nightly_test $SHELL $TEMPFILE
fi
