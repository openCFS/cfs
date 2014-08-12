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
#  Clean up last night's local build directory, since we are low on
#  hdd space.
# ===========================================================================
IF(EXISTS "${HOME}/Documents/dev/NIGHTLY/CFS_BUILD_NIGHTLY")
  FILE(REMOVE_RECURSE "${HOME}/Documents/dev/NIGHTLY/CFS_BUILD_NIGHTLY")
ENDIF()

# ===========================================================================
#  Clean up last night's compiled binaries.
# ===========================================================================
FILE(GLOB CFS_NIGHTLY_ARCHIVES "${CFS_NIGHTLY_DIR}/archives/*")
IF(NOT CFS_NIGHTLY_ARCHIVES STREQUAL "")
  FILE(REMOVE ${CFS_NIGHTLY_ARCHIVES})
ENDIF()

# ===========================================================================
#  Clean trunk test suite.
# ===========================================================================
SET(CFS_TESTSUITE_TRUNK_DIR 
  "${HOME}/Documents/dev/NIGHTLY/CFS_TESTSUITE_TRUNK")

EXECUTE_PROCESS(
  COMMAND svn status
  WORKING_DIRECTORY "${CFS_TESTSUITE_TRUNK_DIR}"
  OUTPUT_VARIABLE SVNST)
STRING(REPLACE "\n" ";" SVNST ${SVNST})

FOREACH(LINE IN ITEMS ${SVNST})
  IF(LINE MATCHES "\\?       ")
    STRING(REPLACE "?       " "${CFS_TESTSUITE_TRUNK_DIR}/" FN ${LINE})
#    MESSAGE("FN ${FN}")
    FILE(REMOVE_RECURSE ${FN})
  ENDIF()
ENDFOREACH()

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
  ctest_update_cfsdeps_trunkOld_wien.cmake

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
#disable this since PW would be submitted to cdash in clean text      
#-DCFS_DS_CDASH_DROP_SITE:STRING=${CFS_DS_CDASH_DROP_SITE}
      -DCFS_TESTUSER:STRING=${CFS_TESTUSER}
      -DCFS_TESTUSER_PW:STRING=${CFS_TESTUSER_PW}
      -S "${SITE_DIR}/${SCRIPT}"
    RESULT_VARIABLE RETVAL
    )

  MESSAGE("RETVAL ${RETVAL}")

ENDFOREACH()

# ===========================================================================
#  Include macros for starting/stopping VBoxes on wiki.
# ===========================================================================
INCLUDE("${SITE_DIR}/macros_for_wiki.cmake")

# ===========================================================================
#  Kill still running VMs...
# ===========================================================================
POWER_OFF_RUNNING_VBOXES()

# ===========================================================================
#  Clean up last night's Virtual Box directories
# ===========================================================================
CLEANUP_VBOX_DIRS()

# ===========================================================================
#  Start today's Virtual Boxes. At 1:30 we just start the most important
#  VBoxes, which generate the nightly binaries for CentOS/RHEL 6 and
#  Ubuntu 12.04. Other VBoxes are started depending on DAYOFWEEK in the
#  script site_specific_finish.cmake.
# ===========================================================================
START_VBOXES("precise;oracle6")

