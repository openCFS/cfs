#-----------------------------------------------------------------------------
# Determine informations in case of a Subversion working copy.
#-----------------------------------------------------------------------------
IF(EXISTS "${CFS_SOURCE_DIR}/.svn")

  #---------------------------------------------------------------------------
  # For information on how to get SVN WC infos at build time
  # cf. http://stackoverflow.com/questions/3780667/
  #       use-cmake-to-get-build-time-svn-revision
  #---------------------------------------------------------------------------
  
  SET(CFS_WC_TYPE "Subversion")
  
  #---------------------------------------------------------------------------
  # If previously configured to use SVNKit we use our own FindSVNKit.cmake.
  #---------------------------------------------------------------------------
  IF(NOT Subversion_SVN_EXECUTABLE MATCHES "jsvn")
    FIND_PACKAGE(Subversion 1.4)
  ELSE(NOT Subversion_SVN_EXECUTABLE MATCHES "jsvn")
    FIND_PACKAGE(SVNKit 1.4)
  ENDIF(NOT Subversion_SVN_EXECUTABLE MATCHES "jsvn")

  # MESSAGE("Subversion_VERSION_SVN ${Subversion_VERSION_SVN}")

  #---------------------------------------------------------------------------
  # Determine version of working copy from first line in .svn/entries and
  # the Repo_VERSION_SVN accordingly in order to compare it to the version
  # of the SVN executable. If the WC version is higher than the version
  # supported by the Subversion exe, we download SVNKit instead.
  #---------------------------------------------------------------------------
  FILE(STRINGS "${CFS_SOURCE_DIR}/.svn/entries" SVNREPOVERSION LIMIT_COUNT 1)
  # MESSAGE("SVNREPOVERSION ${SVNREPOVERSION}")

  SET(Repo_VERSION_SVN "1.4")
  IF(SVNREPOVERSION EQUAL 9)
    SET(Repo_VERSION_SVN "1.5")
  ENDIF(SVNREPOVERSION EQUAL 9)
  IF(SVNREPOVERSION EQUAL 10)
    SET(Repo_VERSION_SVN "1.6")
  ENDIF(SVNREPOVERSION EQUAL 10)
  IF(SVNREPOVERSION EQUAL 11)
    SET(Repo_VERSION_SVN "1.7")
  ENDIF(SVNREPOVERSION EQUAL 11)

  # MESSAGE("Repo_VERSION_SVN ${Repo_VERSION_SVN}")

  #---------------------------------------------------------------------------
  # Check if Subversion has been found and if it has a high enough version,
  # otherwise download and configure SVNKit.
  #---------------------------------------------------------------------------
  IF(NOT Subversion_FOUND OR
     Subversion_VERSION_SVN VERSION_LESS Repo_VERSION_SVN)

   IF(NOT EXISTS "${CFS_BINARY_DIR}/svnkit-1.7.4")
     FILE(DOWNLOAD
       "http://www.svnkit.com/org.tmatesoft.svn_1.7.4.standalone.zip"
       "${CFS_BINARY_DIR}/tmp/org.tmatesoft.svn_1.7.4.standalone.zip")
     
     execute_process(COMMAND unzip "${CFS_BINARY_DIR}/tmp/org.tmatesoft.svn_1.7.4.standalone.zip"
       WORKING_DIRECTORY "${CFS_BINARY_DIR}")
   ENDIF(NOT EXISTS "${CFS_BINARY_DIR}/svnkit-1.7.4")

   SET(Subversion_SVN_EXECUTABLE "${CFS_BINARY_DIR}/svnkit-1.7.4/bin/jsvn" CACHE FILEPATH
     "subversion command line client" FORCE)

   FIND_PACKAGE("SVNKit" 1.4)
  ENDIF(NOT Subversion_FOUND OR
    Subversion_VERSION_SVN VERSION_LESS Repo_VERSION_SVN)

  
  #---------------------------------------------------------------------------
  # Determine informations about CFS WC.
  #---------------------------------------------------------------------------
  Subversion_WC_INFO("${CFS_SOURCE_DIR}" "CFS")

  # MESSAGE("CFS_WC_URL ${CFS_WC_URL}")
  # MESSAGE("CFS_WC_REVISION ${CFS_WC_REVISION}")

  #---------------------------------------------------------------------------
  # A new macro to determine if CFS WC has been modified.
  #---------------------------------------------------------------------------
  MACRO(Subversion_WC_MODIFIED dir prefix)
    # the subversion commands should be executed with the C locale, otherwise
    # the message (which are parsed) may be translated, Alex
    SET(_Subversion_SAVED_LC_ALL "$ENV{LC_ALL}")
    SET(ENV{LC_ALL} C)

    EXECUTE_PROCESS(COMMAND ${Subversion_SVN_EXECUTABLE} status -q ${dir}
      OUTPUT_VARIABLE ${prefix}_WC_MODIFIED
      ERROR_VARIABLE Subversion_svn_modified_error
      RESULT_VARIABLE Subversion_svn_modified_result
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    IF(NOT ${Subversion_svn_modified_result} EQUAL 0)
      MESSAGE(SEND_ERROR "Command \"${Subversion_SVN_EXECUTABLE} status -q ${dir}\" failed with output:\n${Subversion_svn_info_error}")
    ELSE(NOT ${Subversion_svn_modified_result} EQUAL 0)

      IF(${prefix}_WC_MODIFIED)
	SET(${prefix}_WC_MODIFIED "modified")
      ENDIF(${prefix}_WC_MODIFIED)

    ENDIF(NOT ${Subversion_svn_modified_result} EQUAL 0)

    # restore the previous LC_ALL
    SET(ENV{LC_ALL} ${_Subversion_SAVED_LC_ALL})

  ENDMACRO(Subversion_WC_MODIFIED)

  Subversion_WC_MODIFIED("${CFS_SOURCE_DIR}" "CFS")

  # MESSAGE("CFS_WC_MODIFIED ${CFS_WC_MODIFIED}")

ENDIF(EXISTS "${CFS_SOURCE_DIR}/.svn")


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
    SET(CFS_WC_MODIFIED "modified")
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

ENDIF(EXISTS "${CFS_SOURCE_DIR}/.git")

IF(CFS_WC_MODIFIED)
  SET(CFS_WC_REVISION "${CFS_WC_REVISION} ${CFS_WC_MODIFIED}")
ENDIF(CFS_WC_MODIFIED)
