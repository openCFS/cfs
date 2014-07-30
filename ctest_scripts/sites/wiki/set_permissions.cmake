# This scripts is responsible for setting the permissions of all folders
# involved in nightly testing. It needs to be run with root permission
# before nightly testing starts. On wiki, it is run by /etc/crontab.

EXECUTE_PROCESS(
  COMMAND chown -R testuser:testuser /home/testuser/Documents
  RESULT_VARIABLE RETVAL)

EXECUTE_PROCESS(
  COMMAND chown -R testuser:testuser /opt/pckg/cfs_nightly
  RESULT_VARIABLE RETVAL)

EXECUTE_PROCESS(
  COMMAND chown -R testuser:testuser /opt/pckg/CFSDEPSCACHE
  RESULT_VARIABLE RETVAL)

EXECUTE_PROCESS(
  COMMAND chown -R testuser:testuser /opt/pckg/CFSDEPS
  RESULT_VARIABLE RETVAL)

EXECUTE_PROCESS(
  COMMAND chown -R testuser:testuser /var/www/html/cfs_nightly
  RESULT_VARIABLE RETVAL)

