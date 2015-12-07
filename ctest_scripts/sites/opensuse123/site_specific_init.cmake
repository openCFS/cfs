MESSAGE(
"
=============================================================================
 Making a local copy of WC since openSUSE 12.3 uses different SVN...
=============================================================================
"
)

SET(OLD_WC_DIR "${SITE_BASE_DIR}/../../../CFS_TRUNK_NIGHTLY")
SET(DEV_DIR "${SITE_BASE_DIR}/../../../..")
SET(NEW_WC_DIR "${SITE_BASE_DIR}/../../../../CFS_TRUNK_NIGHTLY")

IF(EXISTS ${NEW_WC_DIR})
  FILE(REMOVE_RECURSE ${NEW_WC_DIR})
  FILE(MAKE_DIRECTORY ${NEW_WC_DIR})
ENDIF()

EXECUTE_PROCESS(
#  COMMAND cp -a ${OLD_WC_DIR} ${DEV_DIR}
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${OLD_WC_DIR} ${NEW_WC_DIR}
  RESULT_VARIABLE RETVAL
  )

MESSAGE(
"
=============================================================================
 Upgrading local SVN working copy...
=============================================================================
"
)
EXECUTE_PROCESS(
  COMMAND svn upgrade
  WORKING_DIRECTORY ${NEW_WC_DIR}
  RESULT_VARIABLE RETVAL
  )

