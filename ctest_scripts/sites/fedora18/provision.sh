#!/bin/bash

USER=$1

echo "Starting nightly test on $HOSTNAME..."
# fn=$(mktemp)
#tempfile=$TMP/nightly_test_$(od -N4 -tu /dev/random | awk 'NR==1 {print $2} {}')
TEMPFILE=$(mktemp -t nightly_test_XXXX -u)
echo "TEMPFILE $TEMPFILE"
echo "USER $USER"
PWDLINE=$(fgrep $USER /etc/passwd | grep nightly_test_user)
HOME=$(echo $PWDLINE | cut -d':' -f6)
SHELL=$(echo $PWDLINE | cut -d':' -f7)
echo "HOME $HOME"
echo "SHELL $SHELL"
CMAKE="/opt/pckg/cmake-2.8.10.2-Linux-i386/bin/cmake"
NIGHTLY_DIR="$HOME/Documents/dev/NIGHTLY/CFS_FESPACE_NIGHTLY"
LOG_FILE="$HOME/Documents/dev/nightly_test.log"
echo "LOG_FILE $LOG_FILE"

cat <<EOF > "${TEMPFILE}_2"
#!/bin/sh
$CMAKE -P $NIGHTLY_DIR/ctest_scripts/nightly_test.cmake > "$LOG_FILE" 2>&1
EOF

cat <<EOF > $TEMPFILE
#!/bin/sh
env -i PATH=/bin:/sbin:/usr/bin:/usr/sbin HOME=$HOME SHELL=$SHELL $SHELL "${TEMPFILE}_2"
EOF

sudo -u $USER screen -dmS nightly_test bash $TEMPFILE
