# Assign command line parameters from Vagrantfile to meaningful variable names.
USER=$1
BLOCKING=$2

if [ ! "${BLOCKING}" = "blocking" ]; then
  BLOCKING="nonblocking"
fi

echo "USER $USER"
echo "BLOCKING $BLOCKING"

# Set paths to some programs.
GNU_DATE=/cygdrive/c/GnuWin32/bin/date.exe
SCHTASKS=/cygdrive/c/WINDOWS/system32/schtasks.exe
CMD=/cygdrive/c/WINDOWS/system32/cmd.exe
REG=/cygdrive/c/WINDOWS/system32/reg.exe

SHUTDOWN_MSG="Shutdown machine after nightly tests"
SHUTDOWN_EXE="%SystemRoot%\\System32\\shutdown.exe"

# Build up start command for nightly tests.
CTEST_EXE="\"c:\\Program Files\\cmake-2.8.11-win32-x86\\bin\\ctest.exe\""
NIGHTLY_TEST_CMAKE="z:\\CFS_TRUNK_NIGHTLY\\ctest_scripts\\nightly_test.cmake"
NIGHTLY_TEST_BAT="/cygdrive/c/start_nightly_tests.bat"
LOG_FILE="v:\\logs\\nightly_test.log"
TEST_CMD="$CTEST_EXE -V -S $NIGHTLY_TEST_CMAKE > $LOG_FILE 2>&1"

echo "LOG_FILE $LOG_FILE"

# Function for trimming whitespaces and new lines from strings.
# cf. http://stackoverflow.com/questions/369758/
#            how-to-trim-whitespace-from-bash-variable
trim() {
    # Determine if 'extglob' is currently on.
    local extglobWasOff=1
    shopt extglob >/dev/null && extglobWasOff=0
    # Turn 'extglob' on, if currently turned off.
    (( extglobWasOff )) && shopt -s extglob
    # Trim leading and trailing whitespace
    local var=$1
    var=${var##+([[:space:]])}
    var=${var%%+([[:space:]])}
    # If 'extglob' was off before, turn it back off.
    (( extglobWasOff )) && shopt -u extglob
    echo -n "$var"  # Output trimmed string.
}

# Show Cygwin HOME directory.
echo "HOME (Cygwin): $HOME"

# Set Western European standard time zone in registry.
# cf. http://www.windows-commandline.com/set-time-zone-from-command-line/
echo "Importing time zone information into registry..."
$REG import v:\\timezone.reg

# Make sure logs directory exists before tests are started.
echo "Creating 'logs' directory..."
mkdir -p /cygdrive/v/logs


if [ "$BLOCKING" = "blocking" ]; then

  # In this branch we directly start the nightly test job by running a
  # batch script, which we create manually. The batch script sets up
  # an environment which resembles the standard Windows environment as
  # closely as possible. When nightly, testing is finished, we shutdown
  # the machine manually.

  # Set most important environment variables manually for nightly test
  # since we directly start from the Cygwin environment and not a native
  # one.
  cat <<EOF > $NIGHTLY_TEST_BAT 
set ALLUSERSPROFILE=C:\\Documents and Settings\\All Users
set APPDATA=C:\\Documents and Settings\\user\\Application Data
set CLIENTNAME=Console
set CommonProgramFiles=C:\\Program Files\\Common Files
set COMPUTERNAME=WINXP32
set HOMEDRIVE=C:
set HOMEPATH=\\Documents and Settings\\user
set LOGONSERVER=\\\\WINXP32
set PROMPT=\$P\$G
set SESSIONNAME=Console
set SystemDrive=C:
set SystemRoot=C:\\WINDOWS
set USERDOMAIN=WINXP32
set USERNAME=user
set USERPROFILE=C:\\Documents and Settings\\user
set _=
EOF

  # Read environment settings from Cygwin registry and write them to
  # the nightly batch script (cf. http://smithii.com/node/44)
  pushd . >/dev/null
  for __dir in \
    /proc/registry/HKEY_LOCAL_MACHINE/SYSTEM/CurrentControlSet/Control/Session\ Manager/Environment \
    /proc/registry/HKEY_CURRENT_USER/Environment
  do
    cd "$__dir"
    for __var in *; do
      __var=$(echo $__var)
      echo  "set $__var=$(cat $__var)" >> $NIGHTLY_TEST_BAT
    done
  done
  unset __dir
  unset __var
  popd >/dev/null

  # Write start command for nightly tests to batch script.
  cat <<EOF >> $NIGHTLY_TEST_BAT
$TEST_CMD
EOF

  # Launch new command interpreter with nightly test script
  $CMD /C "c:\\start_nightly_tests.bat"

  # Shut down machine after nightly testing has finished.
  shutdown -f -s now

else # blocking
  
  # In this branch we are in non-blocking mode. Therefore, we register
  # a shutdown task and a nightly test task with the standard Windows
  # task scheduler. Afterwards this script returns immediately and the
  # Vagrant provisioner also returns control to the host environment.

  # Determine shutdown time of machine.
  SHUTDOWN_DATE_TIME=$($GNU_DATE --date='+8 hour' +'%d/%m/%Y %T')
  SHUTDOWN_DATE_TIME=$(trim "$SHUTDOWN_DATE_TIME")
  SHUTDOWN_DATE=$(echo $SHUTDOWN_DATE_TIME | cut -d' ' -f1)
  SHUTDOWN_TIME=$(echo $SHUTDOWN_DATE_TIME | cut -d' ' -f2)
  echo "Machine will shut down at $SHUTDOWN_DATE $SHUTDOWN_TIME."

  # Schedule task for shutting down Windows.
  $SCHTASKS /create /tn "$SHUTDOWN_MSG" /tr "$SHUTDOWN_EXE -f -s" /sc once \
            /st $SHUTDOWN_TIME /sd $SHUTDOWN_DATE /ru Administrator /rp user

  # Determine start time of nightly tests.
  TEST_START_DATE_TIME=$($GNU_DATE --date='+10 minute' +'%d/%m/%Y %T')
  TEST_START_DATE_TIME=$(trim "$TEST_START_DATE_TIME")
  TEST_START_DATE=$(echo $TEST_START_DATE_TIME | cut -d' ' -f1)
  TEST_START_TIME=$(echo $TEST_START_DATE_TIME | cut -d' ' -f2)
  echo "Nightly tests will start at $TEST_START_DATE $TEST_START_TIME."

  cat <<EOF >> $NIGHTLY_TEST_BAT
$TEST_CMD
EOF

  # Schedule task for starting nightly testing.
  $SCHTASKS /create /tn "CTest" /tr "cmd /K c:\\start_nightly_tests.bat" /sc once \
            /st $TEST_START_TIME /sd $TEST_START_DATE /ru user /rp user
fi # blocking
