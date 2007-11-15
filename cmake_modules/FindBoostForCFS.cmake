#-------------------------------------------------------------------------------
# Try to find Boost using CMake standard package.
#-------------------------------------------------------------------------------
FIND_PACKAGE(Boost)

#-------------------------------------------------------------------------------
# If standard algorithm could not find Boost use our own one.
#-------------------------------------------------------------------------------
IF(NOT Boost_FOUND)
  SET(Boost_FOUND 1)

  #-----------------------------------------------------------------------------
  # Search for Boost in ${CFSDEPS_INCLUDE_DIR}, C:/dev/boost, and the paths
  # specified in the registry by the Boost installer
  #-----------------------------------------------------------------------------
  SET (BOOST_POSSIBLE_INCLUDE_PATHS
    ${CFSDEPS_INCLUDE_DIR}
    )

  #-----------------------------------------------------------------------------
  # Look for boost/version.hpp in the possible Boost include paths.
  #-----------------------------------------------------------------------------
  FIND_FILE(Boost_INCLUDE_DIR
    NAMES boost/version.hpp
    PATHS ${BOOST_POSSIBLE_INCLUDE_PATHS}
    )

  IF(NOT Boost_INCLUDE_DIR)
    SET(Boost_FOUND 0)
  ELSE(NOT Boost_INCLUDE_DIR)
    #---------------------------------------------------------------------------
    # Remove boost/version.hpp from Boost_INCLUDE_DIR and set cache value.
    #---------------------------------------------------------------------------
    STRING(REPLACE "/boost/version.hpp"
      ""
      MY_BOOST_INCLUDE_DIR
      ${Boost_INCLUDE_DIR})
    SET(Boost_INCLUDE_DIR ${MY_BOOST_INCLUDE_DIR} CACHE STRING
      "The directory containing the Boost include files." FORCE)

  ENDIF(NOT Boost_INCLUDE_DIR)
ENDIF(NOT Boost_FOUND)

#-------------------------------------------------------------------------------
# On Windows we have to determine the correct suffix for 
# the Boost libraries according to our compiler and build type.
#-------------------------------------------------------------------------------
IF(Boost_FOUND)

  SET(Boost_FOUND 0)
  #---------------------------------------------------------------------------
  # Try to find the Boost library directory and set cache value.
  #---------------------------------------------------------------------------
  IF(EXISTS "${Boost_INCLUDE_DIR}/lib64")
    SET (Boost_LIBRARY_DIR ${Boost_INCLUDE_DIR}/lib64 CACHE STRING
      "The directory containing the Boost libraries." FORCE)
    SET(Boost_FOUND 1)
  ENDIF(EXISTS "${Boost_INCLUDE_DIR}/lib64")

  IF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/lib")
    SET (Boost_LIBRARY_DIR ${Boost_INCLUDE_DIR}/lib CACHE STRING
      "The directory containing the Boost libraries." FORCE)
    SET(Boost_FOUND 1)
  ENDIF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/lib")

  IF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/../lib64")
    SET (Boost_LIBRARY_DIR ${Boost_INCLUDE_DIR}/../lib64 CACHE STRING
      "The directory containing the Boost libraries." FORCE)
    SET(Boost_FOUND 1)
  ENDIF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/../lib64")

  IF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/../lib")
    SET (Boost_LIBRARY_DIR ${Boost_INCLUDE_DIR}/../lib CACHE STRING
      "The directory containing the Boost libraries." FORCE)
    SET(Boost_FOUND 1)
  ENDIF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/../lib")

  IF(Boost_FOUND)

    #-----------------------------------------------------------------------------
    # Specify path to header which has to be preprocessed to get some infos
    # about Boost.
    #-----------------------------------------------------------------------------
    SET(BOOST_TEST_FILE "${CFS_BINARY_DIR}/include/get_boost_infos.hh")
    
    FILE(WRITE ${BOOST_TEST_FILE}
      "#include <boost/config/auto_link.hpp>\n"
      "#include <boost/version.hpp>\n\n"
      "<CFS_MSC_VER>_MSC_VER</CFS_MSC_VER>\n"
      "<CFS_BOOST_VERSION>BOOST_LIB_VERSION</CFS_BOOST_VERSION>"
      )

    #-----------------------------------------------------------------------------
    # Build command line for preprocessing the Boost info header
    #-----------------------------------------------------------------------------
    SET(BOOST_INFO_CMD "${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS}")

    SET(BOOST_INFO_CMD "${BOOST_INFO_CMD} -I${Boost_INCLUDE_DIR} -DBOOST_LIB_NAME=boost_regex -E ${BOOST_TEST_FILE}")

    # MESSAGE("CFS_CXX_COMPILER ${CFS_CXX_COMPILER}")
    # MESSAGE("BOOST_INFO_CMD ${BOOST_INFO_CMD}")

    #-----------------------------------------------------------------------------
    # Preprocess the Boost info header
    #-----------------------------------------------------------------------------
    EXEC_PROGRAM(${BOOST_INFO_CMD}
      ARGS
      OUTPUT_VARIABLE BOOST_INFO_STR
      RETURN_VALUE RETVAL)

    FILE(REMOVE ${BOOST_TEST_FILE})
    #    MESSAGE("BOOST_INFO_STR ${BOOST_INFO_STR}")


    #-----------------------------------------------------------------------------
    # Obtain Boost version.
    #-----------------------------------------------------------------------------
    STRING(REGEX MATCH "<CFS_BOOST_VERSION>.*</CFS_BOOST_VERSION>"
      BOOST_VERSION "${BOOST_INFO_STR}")

    STRING(REGEX REPLACE "\"" ""
      BOOST_VERSION "${BOOST_VERSION}")
    STRING(REGEX MATCH "[0-9]+\\_[0-9]+"
      BOOST_VERSION
      "${BOOST_VERSION}")

    STRING(REGEX REPLACE "_" "."
      CFS_BOOST_VERSION "${BOOST_VERSION}")
#    MESSAGE("CFS_BOOST_VERSION ${CFS_BOOST_VERSION}")

    SET(BOOST_LIB_SUFFIX ".so")
    SET(BOOST_LIB_PREFIX "${Boost_LIBRARY_DIR}/lib")
  ENDIF(Boost_FOUND)

ENDIF(Boost_FOUND)

#-----------------------------------------------------------------------------
# Set Boost library paths in cache.
#-----------------------------------------------------------------------------
SET(BOOST_DATE_TIME_LIB
  "${BOOST_LIB_PREFIX}boost_date_time${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_DATE_TIME_LIB")
MARK_AS_ADVANCED(BOOST_DATE_TIME_LIB)
SET(BOOST_FILESYSTEM_LIB
  "${BOOST_LIB_PREFIX}boost_filesystem${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_FILESYSTEM_LIB")
MARK_AS_ADVANCED(BOOST_FILESYSTEM_LIB)
SET(BOOST_IOSTREAMS_LIB
  "${BOOST_LIB_PREFIX}boost_iostreams${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_IOSTREAMS_LIB")
MARK_AS_ADVANCED(BOOST_IOSTREAMS_LIB)
SET(BOOST_PRG_EXEC_MONITOR_LIB
  "${BOOST_LIB_PREFIX}boost_prg_exec_monitor${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_PRG_EXEC_MONITOR_LIB")
MARK_AS_ADVANCED(BOOST_PRG_EXEC_MONITOR_LIB)
SET(BOOST_PROGRAM_OPTIONS_LIB
  "${BOOST_LIB_PREFIX}boost_program_options${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_PROGRAM_OPTIONS_LIB")
MARK_AS_ADVANCED(BOOST_PROGRAM_OPTIONS_LIB)
SET(BOOST_PYTHON_LIB
  "${BOOST_LIB_PREFIX}boost_python${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_PYTHON_LIB")
MARK_AS_ADVANCED(BOOST_PYTHON_LIB)
SET(BOOST_REGEX_LIB
  "${BOOST_LIB_PREFIX}boost_regex${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_REGEX_LIB")
MARK_AS_ADVANCED(BOOST_REGEX_LIB)
SET(BOOST_SERIALIZATION_LIB
  "${BOOST_LIB_PREFIX}boost_serialization${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_SERIALIZATION_LIB")
MARK_AS_ADVANCED(BOOST_SERIALIZATION_LIB)
SET(BOOST_SIGNALS_LIB
  "${BOOST_LIB_PREFIX}boost_signals${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_SIGNALS_LIB")
MARK_AS_ADVANCED(BOOST_SIGNALS_LIB)
SET(BOOST_TEST_EXEC_MONITOR_LIB
  "${BOOST_LIB_PREFIX}boost_test_exec_monitor${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_TEST_EXEC_MONITOR_LIB")
MARK_AS_ADVANCED(BOOST_TEST_EXEC_MONITOR_LIB)
SET(BOOST_THREAD-MT_LIB
  "${BOOST_LIB_PREFIX}boost_thread-mt${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_THREAD-MT_LIB")
MARK_AS_ADVANCED(BOOST_THREAD-MT_LIB)
SET(BOOST_THREAD_LIB
  "${BOOST_LIB_PREFIX}boost_thread${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_THREAD_LIB")
MARK_AS_ADVANCED(BOOST_THREAD_LIB)
SET(BOOST_UNIT_TEST_FRAMEWORK_LIB
  "${BOOST_LIB_PREFIX}boost_unit_test_framework${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_UNIT_TEST_FRAMEWORK_LIB")
MARK_AS_ADVANCED(BOOST_UNIT_TEST_FRAMEWORK_LIB)
SET(BOOST_WSERIALIZATION_LIB
  "${BOOST_LIB_PREFIX}boost_wserialization${BOOST_LIB_SUFFIX}"
  CACHE STRING "BOOST_WSERIALIZATION_LIB")
MARK_AS_ADVANCED(BOOST_WSERIALIZATION_LIB)


MARK_AS_ADVANCED(Boost_INCLUDE_DIR)
MARK_AS_ADVANCED(Boost_LIBRARY_DIR)

IF(NOT Boost_FOUND)
  MESSAGE("Warning: Boost could not be found! Please specify proper paths.")
ENDIF(NOT Boost_FOUND)
