# ===========================================================================
#  Write SVN auth file to make sure all upcoming svn commands will 
#  work properly.
# ===========================================================================
SET(SVNFILE
  "${HOME}/.subversion/auth/svn.simple/f4c5d049eb353aa17027f06e57e71d0e")

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
<${CFS_DS_HTTPS}> Subversion repository
K 8
username
V 12
${CFS_TESTUSER}
END
"
)

IF(DAYOFWEEK EQUAL 1)
  MESSAGE(
"
=============================================================================
 Removing precompiled CFSDEPS on first day of the week...
=============================================================================
"
  )

  FILE(REMOVE_RECURSE 
    "${HOME}/Documents/dev/NIGHTLY/CFSDEPSCACHE/precompiled")
ENDIF()

# ===========================================================================
#  Clean up last night's compiled binaries.
# ===========================================================================
FILE(GLOB CFS_NIGHTLY_ARCHIVES "${CFS_NIGHTLY_DIR}/archives/*")
IF(NOT CFS_NIGHTLY_ARCHIVES STREQUAL "")
  FILE(REMOVE ${CFS_NIGHTLY_ARCHIVES})
ENDIF()

MESSAGE(
"
=============================================================================
 Update working copies...
=============================================================================
"
)

# Remove Testing directory from modelling manual.
FILE(REMOVE_RECURSE
  "$ENV{HOME}/Documents/dev/NIGHTLY/MODELLING_MANUAL_NIGHTLY/Testing")

SET(UPDATE_SCRIPTS
  # Checkout or update FeSpace working copies
  ctest_update_trunk_wien.cmake
  ctest_update_trunk_testsuite_wien.cmake

  # Checkout or update trunk working copies
  ctest_update_trunkOld_wien.cmake
  ctest_update_trunkOld_testsuite_wien.cmake
  ctest_update_cfsdeps_trunk_wien.cmake

  # Checkout or update modelling manual working copy
  ctest_update_modelling_manual.cmake
)

FOREACH(SCRIPT IN ITEMS ${UPDATE_SCRIPTS})

  # Checkout or update CFS++ FeSpace
  EXECUTE_PROCESS(
    COMMAND ${CTEST_COMMAND} -V
      -DSITE_DIR:PATH=${SITE_DIR}
      -DCFS_DS_HOSTNAME:STRING=${CFS_DS_HOSTNAME}
      -DCFS_DS_SVN:STRING=${CFS_DS_SVN}
      -DCFS_DS_CDASH_DROP_SITE:STRING=${CFS_DS_CDASH_DROP_SITE}
      -DCFS_TESTUSER:STRING=${CFS_TESTUSER}
      -DCFS_TESTUSER_PW:STRING=${CFS_TESTUSER_PW}
      -S "${SITE_DIR}/${SCRIPT}"
    RESULT_VARIABLE RETVAL
    )

  MESSAGE("RETVAL ${RETVAL}")

ENDFOREACH()


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
    # MESSAGE("VBOX ${VBOX}" )
    STRING(REPLACE "{" ";" VBOX ${VBOX})
    LIST(GET VBOX 1 VBOX_UUID)
    LIST(GET VBOX 0 VBOX_NAME)
    STRING(REPLACE "}" ";" VBOX ${VBOX_UUID})
    LIST(GET VBOX 0 VBOX_UUID)

    MESSAGE("Powering off ${VBOX_NAME}- ${VBOX_UUID}" )

    # Power off VBox
    EXECUTE_PROCESS(
      COMMAND VBoxManage controlvm ${VBOX_UUID} poweroff
      RESULT_VARIABLE RETVAL
      )

    MESSAGE("RETVAL ${RETVAL}")

  ENDFOREACH()
ENDIF()

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
  "${HOME}/.vagrant.d"
  "${HOME}/.VirtualBox"
)


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

SET(VBOXES
#  hardy
#  lucid
  precise
  oracle6
#  oracle5
#  fedora18
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

  # Delete unpacked intermediate box files to save hdd space.
  FILE(REMOVE_RECURSE "${HOME}/.vagrant.d/boxes/${VBOX}")

  MESSAGE("RETVAL ${RETVAL}")

ENDFOREACH()

