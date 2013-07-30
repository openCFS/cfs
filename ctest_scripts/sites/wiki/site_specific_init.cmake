MESSAGE("Jipi in site_specific_init.cmake on wiki!")

# Write SVN auth file to make sure all upcoming svn commands will work properly.
SET(SVNFILE "${HOME}/.subversion/auth/svn.simple/f4c5d049eb353aa17027f06e57e71d0e")
FILE(WRITE ${SVNFILE} "K 8
passtype
V 6
simple
K 8
password
V 8
${CFS_TESTUSER_PW}
K 15
svn:realmstring
V 68
<https://lse17.e-technik.uni-erlangen.de:2001> Subversion repository
K 8
username
V 12
${CFS_TESTUSER}
END
"
)

MESSAGE("Update working copies...")

SET(UPDATE_SCRIPTS
  # Checkout or update FeSpace working copies
  ctest_update_fespace_wien.cmake
  ctest_update_fespace_testsuite_wien.cmake

  # Checkout or update trunk working copies
  ctest_update_trunk_wien.cmake
  ctest_update_trunk_testsuite_wien.cmake
  ctest_update_cfsdeps_trunk_wien.cmake
)

FOREACH(SCRIPT IN ITEMS ${UPDATE_SCRIPTS})

  # Checkout or update CFS++ FeSpace
  EXECUTE_PROCESS(
    COMMAND ${CTEST_COMMAND} -V -DHOME:PATH=${HOME} -DSITE_DIR:PATH=${SITE_DIR} -DCFS_TESTUSER:STRING=${CFS_TESTUSER} -DCFS_TESTUSER_PW:STRING=${CFS_TESTUSER_PW} -S "${SITE_DIR}/${SCRIPT}"
  #  OUTPUT_VARIABLE DAYOFWEEK
    RESULT_VARIABLE RETVAL
    )

  MESSAGE("RETVAL ${RETVAL}")

ENDFOREACH()

# Kill still running VMs
EXECUTE_PROCESS(
  COMMAND VBoxManage list runningvms
  OUTPUT_VARIABLE RUNNING_VBOXES
  RESULT_VARIABLE RETVAL
  )

IF(NOT RUNNING_VBOXES STREQUAL "") 
  STRING(REPLACE "\n" ";" RUNNING_VBOXES ${RUNNING_VBOXES})

  FOREACH(VBOX IN ITEMS ${RUNNING_VBOXES})
  #  MESSAGE("VBOX ${VBOX}" )
    STRING(REPLACE "{" ";" VBOX ${VBOX})
    LIST(GET VBOX 1 VBOX_UUID)
    STRING(REPLACE "}" ";" VBOX ${VBOX_UUID})
    LIST(GET VBOX 0 VBOX_UUID)
  #  MESSAGE("VBOX_UUID ${VBOX_UUID}" )

    # Power off VBox
    EXECUTE_PROCESS(
      COMMAND VBoxManage controlvm ${VBOX_UUID} poweroff
    #  OUTPUT_VARIABLE DAYOFWEEK
      RESULT_VARIABLE RETVAL
      )

    MESSAGE("RETVAL ${RETVAL}")

  ENDFOREACH()
ENDIF()

# Clean up last night's Virtual Box directories
FILE(REMOVE_RECURSE
  "${HOME}/VirtualBox VMs"
  "${HOME}/.vagrant.d"
  "${HOME}/.VirtualBox"
)

MESSAGE("Starting VBoxes...")
SET(VBOXES
  lucid
  precise
#  oracle6
)

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

  MESSAGE("RETVAL ${RETVAL}")

ENDFOREACH()

