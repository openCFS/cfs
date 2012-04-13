# FIND_PACKAGE("Java" 1.4)

MARK_AS_ADVANCED(Subversion_SVN_EXECUTABLE)

IF(Subversion_SVN_EXECUTABLE)
  # the subversion commands should be executed with the C locale, otherwise
  # the message (which are parsed) may be translated, Alex
  SET(_Subversion_SAVED_LC_ALL "$ENV{LC_ALL}")
  SET(ENV{LC_ALL} C)

  EXECUTE_PROCESS(COMMAND ${Subversion_SVN_EXECUTABLE} --version
    OUTPUT_VARIABLE Subversion_VERSION_SVN
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # restore the previous LC_ALL
  SET(ENV{LC_ALL} ${_Subversion_SAVED_LC_ALL})

  STRING(REGEX REPLACE "jsvn, version ([.0-9]+).*"
    "\\1" Subversion_VERSION_SVN "${Subversion_VERSION_SVN}")

#  MESSAGE(FATAL_ERROR "Subversion_VERSION_SVN ${Subversion_VERSION_SVN}")

  MACRO(Subversion_WC_INFO dir prefix)
    # the subversion commands should be executed with the C locale, otherwise
    # the message (which are parsed) may be translated, Alex
    SET(_Subversion_SAVED_LC_ALL "$ENV{LC_ALL}")
    SET(ENV{LC_ALL} C)

    EXECUTE_PROCESS(COMMAND ${Subversion_SVN_EXECUTABLE} info ${dir}
      OUTPUT_VARIABLE ${prefix}_WC_INFO
      ERROR_VARIABLE Subversion_svn_info_error
      RESULT_VARIABLE Subversion_svn_info_result
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    IF(NOT ${Subversion_svn_info_result} EQUAL 0)
      MESSAGE(SEND_ERROR "Command \"${Subversion_SVN_EXECUTABLE} info ${dir}\" failed with output:\n${Subversion_svn_info_error}")
    ELSE(NOT ${Subversion_svn_info_result} EQUAL 0)

      STRING(REGEX REPLACE "^(.*\n)?URL: ([^\n]+).*"
        "\\2" ${prefix}_WC_URL "${${prefix}_WC_INFO}")
      STRING(REGEX REPLACE "^(.*\n)?Repository Root: ([^\n]+).*"
        "\\2" ${prefix}_WC_ROOT "${${prefix}_WC_INFO}")
      STRING(REGEX REPLACE "^(.*\n)?Revision: ([^\n]+).*"
        "\\2" ${prefix}_WC_REVISION "${${prefix}_WC_INFO}")
      STRING(REGEX REPLACE "^(.*\n)?Last Changed Author: ([^\n]+).*"
        "\\2" ${prefix}_WC_LAST_CHANGED_AUTHOR "${${prefix}_WC_INFO}")
      STRING(REGEX REPLACE "^(.*\n)?Last Changed Rev: ([^\n]+).*"
        "\\2" ${prefix}_WC_LAST_CHANGED_REV "${${prefix}_WC_INFO}")
      STRING(REGEX REPLACE "^(.*\n)?Last Changed Date: ([^\n]+).*"
        "\\2" ${prefix}_WC_LAST_CHANGED_DATE "${${prefix}_WC_INFO}")

    ENDIF(NOT ${Subversion_svn_info_result} EQUAL 0)

    # restore the previous LC_ALL
    SET(ENV{LC_ALL} ${_Subversion_SAVED_LC_ALL})

  ENDMACRO(Subversion_WC_INFO)

  MACRO(Subversion_WC_LOG dir prefix)
    # This macro can block if the certificate is not signed:
    # svn ask you to accept the certificate and wait for your answer
    # This macro requires a svn server network access (Internet most of the time)
    # and can also be slow since it access the svn server
    EXECUTE_PROCESS(COMMAND
      ${Subversion_SVN_EXECUTABLE} log -r BASE ${dir}
      OUTPUT_VARIABLE ${prefix}_LAST_CHANGED_LOG
      ERROR_VARIABLE Subversion_svn_log_error
      RESULT_VARIABLE Subversion_svn_log_result
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    IF(NOT ${Subversion_svn_log_result} EQUAL 0)
      MESSAGE(SEND_ERROR "Command \"${Subversion_SVN_EXECUTABLE} log -r BASE ${dir}\" failed with output:\n${Subversion_svn_log_error}")
    ENDIF(NOT ${Subversion_svn_log_result} EQUAL 0)
  ENDMACRO(Subversion_WC_LOG)

ENDIF(Subversion_SVN_EXECUTABLE)

INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Subversion REQUIRED_VARS Subversion_SVN_EXECUTABLE
                                             VERSION_VAR Subversion_VERSION_SVN )

# for compatibility
SET(Subversion_FOUND ${SUBVERSION_FOUND})
SET(Subversion_SVN_FOUND ${SUBVERSION_FOUND})
