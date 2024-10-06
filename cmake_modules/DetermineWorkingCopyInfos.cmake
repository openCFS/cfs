#-----------------------------------------------------------------------------
# Determine informations in case of a Git working copy.
#-----------------------------------------------------------------------------
IF(EXISTS "${CFS_SOURCE_DIR}/.git")
  FIND_PACKAGE("Git" 1.4)

  SET(CFS_WC_TYPE "Git")

  SET(ENV{GIT_DIR} "${CFS_SOURCE_DIR}/.git")
  SET(ENV{GIT_WORK_TREE} "${CFS_SOURCE_DIR}")

  SET(_Subversion_SAVED_LC_ALL "$ENV{LC_ALL}")
  SET(ENV{LC_ALL} C)

  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} svn info
    OUTPUT_VARIABLE CFS_GIT_SVN_INFO
    ERROR_VARIABLE CFS_GIT_SVN_ERROR
    RESULT_VARIABLE RETVAL)

  IF(RETVAL)
      SET(CFS_WC_URL "N/A")
      SET(CFS_WC_ROOT "N/A")
      SET(CFS_WC_REVISION "N/A")
  ELSE()
    STRING(REGEX REPLACE "^(.*\n)?URL: ([^\n]+).*"
      "\\2" CFS_WC_URL "${CFS_GIT_SVN_INFO}")
    STRING(REGEX REPLACE "^(.*\n)?Repository Root: ([^\n]+).*"
      "\\2" CFS_WC_ROOT "${CFS_GIT_SVN_INFO}")
    STRING(REGEX REPLACE "^(.*\n)?Revision: ([^\n]+).*"
      "\\2" CFS_WC_REVISION "${CFS_GIT_SVN_INFO}")
  ENDIF()

#  MESSAGE("CFS_WC_URL ${CFS_WC_URL}")
#  MESSAGE("CFS_WC_REVISION ${CFS_WC_REVISION}")
#  MESSAGE("CFS_WC_ROOT ${CFS_WC_ROOT}")

  #---------------------------------------------------------------------------
  # Find out git commit identifier.
  #---------------------------------------------------------------------------
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} log -1
    OUTPUT_VARIABLE CFS_GIT_LOG
    RESULT_VARIABLE RETVAL)

#  MESSAGE("CFS_GIT_LOG ${CFS_GIT_LOG}")

  STRING(REGEX REPLACE "^(.*\n)?commit ([^\n]+).*"
    "\\2" CFS_GIT_COMMIT "${CFS_GIT_LOG}")

#  MESSAGE("CFS_GIT_COMMIT ${CFS_GIT_COMMIT}")


  #---------------------------------------------------------------------------
  # Check if WC has been modified.
  #---------------------------------------------------------------------------
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} diff
    OUTPUT_VARIABLE CFS_WC_MODIFIED
    RESULT_VARIABLE RETVAL)

#  MESSAGE("CFS_GIT_MODIFIED ${CFS_WC_MODIFIED}")
  
  IF(CFS_WC_MODIFIED)
    SET(CFS_WC_MODIFIED "(modified)")
  ENDIF(CFS_WC_MODIFIED)

#  MESSAGE("CFS_WC_MODIFIED ${CFS_WC_MODIFIED}")

  #---------------------------------------------------------------------------
  # Determine name of local git branch.
  #---------------------------------------------------------------------------
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} branch
    OUTPUT_VARIABLE CFS_GIT_BRANCH_LIST
    RESULT_VARIABLE RETVAL)

#  MESSAGE("CFS_GIT_BRANCH_LIST ${CFS_GIT_BRANCH_LIST}")
  
  STRING(REGEX REPLACE "^(.*\n)?\\* ([^\n]+).*"
    "\\2" CFS_GIT_BRANCH "${CFS_GIT_BRANCH_LIST}")

#  MESSAGE("CFS_GIT_BRANCH ${CFS_GIT_BRANCH}")

  # restore the previous LC_ALL
  SET(ENV{LC_ALL} ${_Subversion_SAVED_LC_ALL})

  # set 
  SET(CFS_WC_REVISION ${CFS_GIT_COMMIT})

ENDIF(EXISTS "${CFS_SOURCE_DIR}/.git")

IF(CFS_WC_MODIFIED)
  SET(CFS_WC_REVISION "${CFS_WC_REVISION} ${CFS_WC_MODIFIED}")
ENDIF(CFS_WC_MODIFIED)
