echo "Hello World!"

echo $HOME

ls $HOME

GNU_DATE=/cygdrive/c/GnuWin32/bin/date.exe
SCHTASKS=/cygdrive/c/WINDOWS/system32/schtasks.exe

SHUTDOWN_MSG="Shutdown machine after nightly tests"
SHUTDOWN_EXE="%SystemRoot%\\System32\\shutdown.exe"

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

# Determine shutdown time of machine.
SHUTDOWN_DATE_TIME=$($GNU_DATE --date='+8 hour' +'%d/%m/%Y %T')
SHUTDOWN_DATE_TIME=$(trim "$SHUTDOWN_DATE_TIME")
SHUTDOWN_DATE=$(echo $SHUTDOWN_DATE_TIME | cut -d' ' -f1)
SHUTDOWN_TIME=$(echo $SHUTDOWN_DATE_TIME | cut -d' ' -f2)
echo "Machine will shut down at $SHUTDOWN_DATE $SHUTDOWN_TIME."

# Determine start time of nightly tests.
TEST_START_DATE_TIME=$($GNU_DATE --date='+10 minute' +'%d/%m/%Y %T')
TEST_START_DATE_TIME=$(trim "$TEST_START_DATE_TIME")
TEST_START_DATE=$(echo $TEST_START_DATE_TIME | cut -d' ' -f1)
TEST_START_TIME=$(echo $TEST_START_DATE_TIME | cut -d' ' -f2)
echo "Nightly tests will start at $TEST_START_DATE $TEST_START_TIME."

DUMMY_DATE_TIME=$($GNU_DATE --date='+20 second' +'%d/%m/%Y %T')
DUMMY_DATE_TIME=$(trim "$DUMMY_DATE_TIME")
DUMMY_DATE=$(echo $DUMMY_DATE_TIME | cut -d' ' -f1)
DUMMY_TIME=$(echo $DUMMY_DATE_TIME | cut -d' ' -f2)

echo "***$DUMMY_DATE###$DUMMY_TIME---"

# $SCHTASKS /create /tn "$SHUTDOWN_MSG" /tr "$SHUTDOWN_EXE -f -s" /sc once \
#           /st $SHUTDOWN_TIME /sd $SHUTDOWN_DATE /ru Administrator /rp user

# $SCHTASKS /create /tn "Show env" /tr "cmd /K set" /sc once \
#         /st $DUMMY_TIME /sd $DUMMY_DATE /ru user /rp user

$SCHTASKS /create /tn "CTest" /tr "cmd /K ctest -V -S z:\\CFS_FESPACE_NIGHTLY\\ctest_scripts\\nightly_test.cmake" /sc once \
         /st $DUMMY_TIME /sd $DUMMY_DATE /ru user /rp user

echo "==============="
/cygdrive/c/GnuWin32/bin/dir y:
echo "==============="
/cygdrive/c/GnuWin32/bin/dir v:
echo "==============="
