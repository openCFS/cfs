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
  # Determine the svn format of the working copy and compare it to the svn
  # format used by the SVN executable. If the working copy svn format is 
  # higher than the format supported by the SVN executable, we download the
  # SVNKit instead.
  # http://stackoverflow.com/questions/1364618/how-do-i-determine-the-svn-working-copy-layout-version
  # http://stackoverflow.com/questions/19265363/find-out-svn-working-copy-version-1-7-or-1-8
  # http://svn.apache.org/repos/asf/subversion/trunk/subversion/libsvn_wc/wc.h
  #---------------------------------------------------------------------------
  SET(SVN_REPO_FORMAT 0)
  IF(EXISTS "${CFS_SOURCE_DIR}/.svn/format")
    FILE(STRINGS "${CFS_SOURCE_DIR}/.svn/format" SVN_REPO_FORMAT_TMP LIMIT_COUNT 1)
    IF(SVN_REPO_FORMAT_TMP GREATER SVN_REPO_FORMAT)
      SET(SVN_REPO_FORMAT ${SVN_REPO_FORMAT_TMP})
    ENDIF(SVN_REPO_FORMAT_TMP GREATER SVN_REPO_FORMAT)
  ENDIF(EXISTS "${CFS_SOURCE_DIR}/.svn/format")

  IF(EXISTS "${CFS_SOURCE_DIR}/.svn/entries")
    FILE(STRINGS "${CFS_SOURCE_DIR}/.svn/entries" SVN_REPO_FORMAT_TMP LIMIT_COUNT 1)
    IF(SVN_REPO_FORMAT_TMP GREATER SVN_REPO_FORMAT)
      SET(SVN_REPO_FORMAT ${SVN_REPO_FORMAT_TMP})
    ENDIF(SVN_REPO_FORMAT_TMP GREATER SVN_REPO_FORMAT)
  ENDIF(EXISTS "${CFS_SOURCE_DIR}/.svn/entries")

  IF(EXISTS "${CFS_SOURCE_DIR}/.svn/wc.db")
    EXECUTE_PROCESS(
      COMMAND od -An -j63 -N1 -t dC "${CFS_SOURCE_DIR}/.svn/wc.db"
      OUTPUT_VARIABLE SVN_REPO_FORMAT_TMP
      RESULT_VARIABLE RETVAL
    )
    IF(NOT RETVAL AND SVN_REPO_FORMAT_TMP GREATER SVN_REPO_FORMAT)
      STRING(STRIP ${SVN_REPO_FORMAT_TMP} SVN_REPO_FORMAT_TMP)
      SET(SVN_REPO_FORMAT ${SVN_REPO_FORMAT_TMP})
    ENDIF(NOT RETVAL AND SVN_REPO_FORMAT_TMP GREATER SVN_REPO_FORMAT)
  ENDIF(EXISTS "${CFS_SOURCE_DIR}/.svn/wc.db")

  IF(SVN_REPO_FORMAT EQUAL 0)
    MESSAGE(FATAL_ERROR "Could not read SVN repository format from ${CFS_SOURCE_DIR}/.svn/format, ${CFS_SOURCE_DIR}/.svn/entries or ${CFS_SOURCE_DIR}/.svn/wc.db")
  ENDIF(SVN_REPO_FORMAT EQUAL 0)

  SET(SVN_REPO_FORMAT_MATCHES FALSE)
  IF(Subversion_FOUND)
    SET(Subversion_USED_FORMAT 1)
    IF(Subversion_VERSION_SVN VERSION_EQUAL 1.0 OR Subversion_VERSION_SVN VERSION_GREATER 1.0)
      SET(Subversion_USED_FORMAT 4)
    ENDIF(Subversion_VERSION_SVN VERSION_EQUAL 1.0 OR Subversion_VERSION_SVN VERSION_GREATER 1.0)
    IF(Subversion_VERSION_SVN VERSION_EQUAL 1.4 OR Subversion_VERSION_SVN VERSION_GREATER 1.4)
      SET(Subversion_USED_FORMAT 8)
    ENDIF(Subversion_VERSION_SVN VERSION_EQUAL 1.4 OR Subversion_VERSION_SVN VERSION_GREATER 1.4)
    IF(Subversion_VERSION_SVN VERSION_EQUAL 1.5 OR Subversion_VERSION_SVN VERSION_GREATER 1.5)
      SET(Subversion_USED_FORMAT 9)
    ENDIF(Subversion_VERSION_SVN VERSION_EQUAL 1.5 OR Subversion_VERSION_SVN VERSION_GREATER 1.5)
    IF(Subversion_VERSION_SVN VERSION_EQUAL 1.6 OR Subversion_VERSION_SVN VERSION_GREATER 1.6)
      SET(Subversion_USED_FORMAT 10)
    ENDIF(Subversion_VERSION_SVN VERSION_EQUAL 1.6 OR Subversion_VERSION_SVN VERSION_GREATER 1.6)
    IF(Subversion_VERSION_SVN VERSION_EQUAL 1.7 OR Subversion_VERSION_SVN VERSION_GREATER 1.7)
      SET(Subversion_USED_FORMAT 29)
    ENDIF(Subversion_VERSION_SVN VERSION_EQUAL 1.7 OR Subversion_VERSION_SVN VERSION_GREATER 1.7)
    IF(Subversion_VERSION_SVN VERSION_EQUAL 1.8 OR Subversion_VERSION_SVN VERSION_GREATER 1.8)
      SET(Subversion_USED_FORMAT 31)
    ENDIF(Subversion_VERSION_SVN VERSION_EQUAL 1.8 OR Subversion_VERSION_SVN VERSION_GREATER 1.8)
    IF(Subversion_USED_FORMAT EQUAL SVN_REPO_FORMAT OR Subversion_USED_FORMAT GREATER SVN_REPO_FORMAT)
      SET(SVN_REPO_FORMAT_MATCHES TRUE)
    ENDIF(Subversion_USED_FORMAT EQUAL SVN_REPO_FORMAT OR Subversion_USED_FORMAT GREATER SVN_REPO_FORMAT)
  ENDIF(Subversion_FOUND)

  #--------------------------------------------------------------------------
  # Check if Subversion has been found and if it supports the repository svn
  # format download and configure SVNKit.
  # Also use SVNKit for smaller repo versions, needed for nightly builds on
  # virtual boxes.
  #---------------------------------------------------------------------------
  IF(NOT Subversion_FOUND OR NOT SVN_REPO_FORMAT_MATCHES)

   IF(NOT EXISTS "${CFS_BINARY_DIR}/svnkit-1.8.10")
     FILE(DOWNLOAD
       "http://www.svnkit.com/org.tmatesoft.svn_1.8.10.standalone.zip"
       "${CFS_BINARY_DIR}/tmp/org.tmatesoft.svn_1.8.10.standalone.zip")

     #-----------------------------------------------------------------------------
     # Workaround for cmake download bug, may appear related to
     # https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=716798
     # http://public.kitware.com/pipermail/cmake-developers/2013-July/007759.html
     #-----------------------------------------------------------------------------
     EXECUTE_PROCESS(COMMAND du -k "${CFS_BINARY_DIR}/tmp/org.tmatesoft.svn_1.8.10.standalone.zip" 
       OUTPUT_VARIABLE SVN_FILESIZE)
     STRING(REGEX REPLACE "${CFS_BINARY_DIR}/tmp/org.tmatesoft.svn_1.8.10.standalone.zip" "" SVN_FILESIZE ${SVN_FILESIZE})
     STRING(STRIP ${SVN_FILESIZE} SVN_FILESIZE)
     IF(SVN_FILESIZE EQUAL 0)
       EXECUTE_PROCESS(COMMAND wget -q "http://www.svnkit.com/org.tmatesoft.svn_1.8.10.standalone.zip" -O "${CFS_BINARY_DIR}/tmp/org.tmatesoft.svn_1.8.10.standalone.zip")
     ENDIF(SVN_FILESIZE EQUAL 0)
  
     
     EXECUTE_PROCESS(COMMAND unzip -qq
       "${CFS_BINARY_DIR}/tmp/org.tmatesoft.svn_1.8.10.standalone.zip"
       WORKING_DIRECTORY "${CFS_BINARY_DIR}")
   ENDIF(NOT EXISTS "${CFS_BINARY_DIR}/svnkit-1.8.10")

   SET(Subversion_SVN_EXECUTABLE "${CFS_BINARY_DIR}/svnkit-1.8.10/bin/jsvn"
     CACHE FILEPATH
     "subversion command line client" FORCE)

   FIND_PACKAGE("SVNKit" 1.4)
  ENDIF(NOT Subversion_FOUND OR NOT SVN_REPO_FORMAT_MATCHES)

  
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
      SET(MSG "Command \"${Subversion_SVN_EXECUTABLE} status -q ${dir}\"")
      SET(MSG "${MSG} failed with output:\n${Subversion_svn_info_error}")
      MESSAGE(SEND_ERROR "${MSG}")
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
    set(CFS_GIT_COMMIT "${CFS_GIT_COMMIT} (${CFS_WC_MODIFIED})")
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
