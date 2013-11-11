# This script is for logging the disk usage of the partitions
# involved in the nightly build. We sample every sixty seconds
# starting from 1:00 AM until 3:00 PM.
#
# The output gets written to a .csv file in the test user's 
# $HOME/Documents/dev directory.

# We set the correct time zone and locale information. 
SET(ENV{TZ} "Europe/Vienna")
SET(ENV{LANG} "en_US.UTF-8")

# Let's specify the sampling interval in seconds and the number
# of sampling points.
SET(INTERVAL 30) # 60 seconds
SET(COUNTER 1680)

# Set path to .csv log file
SET(DI_LOG_FILE "/home/testuser/Documents/dev/di_log.csv")

# Write header line if log file did not exist before.
IF(NOT EXISTS ${DI_LOG_FILE})
  FILE(WRITE ${DI_LOG_FILE} "\"Time\",\"Mount point /home\",\"Total Space /home\",\"Used Space /home\",\"Percentage /home\",\"Mount point /var\",\"Total Space /var\",\"Used Space /var\",\"Percentage /var\"\n")
ENDIF()

# Loop over COUNTER.
while(COUNTER)
#  MESSAGE("COUNTER ${COUNTER}")

  math(EXPR COUNTER "${COUNTER} - 1")

  # Get time in format MM/DD/YYYY HH:MM:SS
  EXECUTE_PROCESS(
    COMMAND date "+%m/%d/%Y %T"
    OUTPUT_VARIABLE DATE_OUT
  )
  STRING(STRIP ${DATE_OUT} DATE_OUT)

#  MESSAGE("DATE_OUT ###${DATE_OUT}###")
  SET(CSV_LINE ${DATE_OUT})

  # Get disk information in CSV format.
  EXECUTE_PROCESS(
    COMMAND di -f mbup -c -l -n
    OUTPUT_VARIABLE DI_OUT
  )

  STRING(REPLACE "\n" ";" DI_OUT ${DI_OUT})

  # We are just interested in the line for /home and /var.
  FOREACH(LINE IN ITEMS ${DI_OUT})
    IF(LINE MATCHES "home" OR
       LINE MATCHES "var")
#      MESSAGE("LINE ${LINE}")
      SET(CSV_LINE "${CSV_LINE},${LINE}")
    ENDIF()
  ENDFOREACH()

  # Append new line to log file.
  FILE(APPEND ${DI_LOG_FILE} "${CSV_LINE}\n")

  # Sleep for INTERVAL seconds.
  CTEST_SLEEP(${INTERVAL})
endwhile(COUNTER)

