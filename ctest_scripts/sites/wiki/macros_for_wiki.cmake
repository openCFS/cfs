MACRO(GET_SECONDS OUTVAR TIME)
  SET(CMD "date")

  IF(NOT ${TIME} STREQUAL "now")
    LIST(APPEND CMD "-d" "${TIME}")
  ENDIF()

  LIST(APPEND CMD "+%s")

  # MESSAGE("CMD ${CMD}")

  # Find out seconds
  EXECUTE_PROCESS(
    COMMAND ${CMD}
    OUTPUT_VARIABLE TMP
    RESULT_VARIABLE RETVAL)
  STRING(STRIP ${TMP} TMP)

  SET(${OUTVAR} ${TMP})
ENDMACRO()

MACRO(COPY_ZIPS_TO_APACHE)
  FILE(GLOB NIGHTLY_ZIPS 
    "${CFS_NIGHTLY_DIR}/archives/*.zip"
    "${CFS_NIGHTLY_DIR}/archives/*.dmg"
    "${CFS_NIGHTLY_DIR}/archives/README.html")

  FOREACH(ZIP IN ITEMS ${NIGHTLY_ZIPS})
    get_filename_component(FN "${ZIP}" NAME)

    FILE(MD5 "${ZIP}" MD5_INPUT)
    IF(NOT EXISTS "${APACHE_NIGHTLY_DIR}/${FN}")
      SET(MD5_OUTPUT "EMPTY")
    ELSE()
      FILE(MD5 "${APACHE_NIGHTLY_DIR}/${FN}" MD5_OUTPUT)
    ENDIF()

    IF(NOT MD5_INPUT STREQUAL MD5_OUTPUT) 

      EXECUTE_PROCESS(
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${ZIP}" "${APACHE_NIGHTLY_DIR}/${FN}"
        RESULT_VARIABLE RETVAL)

      MESSAGE("Submitting ${FN} to development server...")
      EXECUTE_PROCESS(
        COMMAND curl -u ${CFS_TESTUSER}:${CFS_TESTUSER_PW} 
                     -k -T ${FN}
                     ${CFS_DS_WEBDAV}/nightly-builds/${FN}
        WORKING_DIRECTORY "${NIGHTLY_ARCHIVES_DIR}"
        RESULT_VARIABLE RETVAL
      )

    ENDIF()
    UNSET(FN)
    UNSET(MD5_INPUT)
    UNSET(MD5_OUTPUT)
  ENDFOREACH()
ENDMACRO()


MACRO(POWER_OFF_RUNNING_VBOXES)
# ===========================================================================
#  Kill still running VMs...
# ===========================================================================
MESSAGE(
"
=============================================================================
 Power off still running VMs...
=============================================================================
"
)

EXECUTE_PROCESS(
  COMMAND VBoxManage list runningvms
  OUTPUT_VARIABLE RUNNING_VBOXES
  RESULT_VARIABLE RETVAL
  )

IF(NOT RUNNING_VBOXES STREQUAL "") 
  STRING(REPLACE "\n" ";" RUNNING_VBOXES ${RUNNING_VBOXES})

  FOREACH(VBOX IN ITEMS ${RUNNING_VBOXES})
    # Obtain VBox name and UUID.
    STRING(REPLACE "{" ";" VBOX ${VBOX})
    LIST(GET VBOX 1 VBOX_UUID)
    LIST(GET VBOX 0 VBOX_NAME)
    STRING(REPLACE "}" ";" VBOX ${VBOX_UUID})
    LIST(GET VBOX 0 VBOX_UUID)

    # Remove quotation marks from VBox name.
    STRING(REPLACE "\"" ";" VBOX ${VBOX_NAME})
    LIST(GET VBOX 1 VBOX_NAME)

    MESSAGE("Powering off ${VBOX_NAME} - ${VBOX_UUID}" )

    # Find out PID of running VBox by reading its log file.
    FILE(READ "$ENV{HOME}/VirtualBox VMs/${VBOX_NAME}/Logs/VBox.log" VBOX_LOG)
    STRING(REPLACE "\n" ";" VBOX_LOG ${VBOX_LOG})

    FOREACH(LINE IN ITEMS ${VBOX_LOG})
      IF(LINE MATCHES "Process ID")
#        MESSAGE("LINE ${LINE}")
        STRING(REPLACE " " ";" LINE ${LINE})
        LIST(GET LINE 3 VBOX_PID)
        STRING(STRIP ${VBOX_PID} VBOX_PID)
        MESSAGE("VBOX_PID ${VBOX_PID}")
      ENDIF()
    ENDFOREACH()

    # Try to power off VBox gracefully
    EXECUTE_PROCESS(
      COMMAND VBoxManage controlvm ${VBOX_UUID} poweroff
      TIMEOUT 300
      RESULT_VARIABLE RETVAL
      )
    MESSAGE("RETVAL ${RETVAL}")

    # If graceful shutdown did not work, kill the VBox process manually.
    IF(NOT RETVAL EQUAL 0 AND
       VBOX_PID)
      # Let's first make sure we hit the right process and check if it
      # is still running.
      EXECUTE_PROCESS(
        COMMAND ps -f -p ${VBOX_PID}
        OUTPUT_VARIABLE PS_INFO
        RESULT_VARIABLE RETVAL
      )
      IF(RETVAL EQUAL 0)
        # Process still running
        STRING(REPLACE "\n" ";" PS_INFO_LIST ${PS_INFO})
        LIST(GET PS_INFO_LIST 1 PS_INFO)

        # Check if it really our VBox.
        IF(PS_INFO MATCHES ${VBOX_UUID})
          MESSAGE("Let's kill the bastard ${VBOX_NAME}!")
          EXECUTE_PROCESS(
            COMMAND kill -9 ${VBOX_PID}
          )
        ELSE()
          MESSAGE("No running process ${VBOX_PID} for ${VBOX_NAME} found!")
        ENDIF()
      ENDIF()
    ENDIF()    

  ENDFOREACH()
ENDIF()
ENDMACRO()

MACRO(CLEANUP_VBOX_DIRS)
# ===========================================================================
#  Clean up last night's Virtual Box directories
# ===========================================================================
MESSAGE(
"
=============================================================================
 Clean up last night's Virtual Box directories...
=============================================================================
"
)
FILE(REMOVE_RECURSE
  "${HOME}/VirtualBox VMs"
  "$ENV{VAGRANT_HOME}"
  "${HOME}/.vagrant.d"
  "${HOME}/.VirtualBox"
)
ENDMACRO()

MACRO(START_VBOXES VBOXES)
# ===========================================================================
#  Start today's Virtual Boxes
# ===========================================================================
MESSAGE(
"
=============================================================================
 Starting VBoxes...
=============================================================================
"
)

#SET(VBOXES
#  hardy
#  lucid
#  precise
#  oracle6
#  oracle5
#  fedora18
#)

FOREACH(VBOX IN ITEMS ${VBOXES})

  FILE(REMOVE_RECURSE
    "${SITE_BASE_DIR}/${VBOX}/.vagrant"
  )

  # Run vagrant up in site dirs
  EXECUTE_PROCESS(
    COMMAND vagrant up
    WORKING_DIRECTORY "${SITE_BASE_DIR}/${VBOX}"
    RESULT_VARIABLE RETVAL
    )

  # Delete unpacked intermediate box files to save hdd space.
  FILE(REMOVE_RECURSE "$ENV{VAGRANT_HOME}/boxes/${VBOX}")

  MESSAGE("RETVAL ${RETVAL}")

ENDFOREACH()
ENDMACRO()
