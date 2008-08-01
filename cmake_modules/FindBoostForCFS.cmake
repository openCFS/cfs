SET(Boost_FOUND 1)

#-------------------------------------------------------------------------------
# Determine paths of Boost libraries.
#-------------------------------------------------------------------------------
SET (Boost_POSSIBLE_LIB_PATHS
  ${CFSDEPS_LIBRARY_DIR}
  /usr/lib64
  /usr/lib
  /usr/local/lib64
  /usr/local/lib
)

FIND_LIBRARY(BOOST_DATE_TIME_LIB_TEST
  NAMES boost_date_time
  PATHS ${Boost_POSSIBLE_LIB_PATHS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  )

SET(BOOST_DATE_TIME_LIB_TEST "${BOOST_DATE_TIME_LIB_TEST}"
  CACHE INTERNAL "BOOST_DATE_TIME_LIB_TEST")

#-------------------------------------------------------------------------------
# Look for Boost header.
#-------------------------------------------------------------------------------
SET (Boost_POSSIBLE_INCLUDE_PATHS
  ${CFSDEPS_INCLUDE_DIR}
  /usr/include
  /usr/local/include
  )

FIND_PATH(Boost_INCLUDE_DIR
  NAMES boost/version.hpp 
  PATHS ${Boost_POSSIBLE_INCLUDE_PATHS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  )

MARK_AS_ADVANCED(Boost_INCLUDE_DIR)


IF(BOOST_DATE_TIME_LIB_TEST AND Boost_INCLUDE_DIR)
  SET(Boost_FOUND 1)

  STRING(REGEX REPLACE "/libboost.*" ""
    Boost_LIBRARY_DIR "${BOOST_DATE_TIME_LIB_TEST}")

  SET (Boost_LIBRARY_DIR "${Boost_LIBRARY_DIR}" CACHE STRING
    "The directory containing the Boost libraries." FORCE)

ELSE(BOOST_DATE_TIME_LIB_TEST AND Boost_INCLUDE_DIR)
  #-----------------------------------------------------------------------------
  # Try to find Boost using CMake standard package.
  #-----------------------------------------------------------------------------
  FIND_PACKAGE(Boost)

  #---------------------------------------------------------------------------
  # Try to find the Boost library directory and set cache value.
  #---------------------------------------------------------------------------
  IF(NOT Boost_LIBRARY_DIR)
    IF(CFS_ARCH STREQUAL "X86_64")
      IF(EXISTS "${Boost_INCLUDE_DIR}/lib64")
	SET (Boost_LIBRARY_DIR ${Boost_INCLUDE_DIR}/lib64 CACHE STRING
	  "The directory containing the Boost libraries." FORCE)
	SET(Boost_FOUND 1)
      ENDIF(EXISTS "${Boost_INCLUDE_DIR}/lib64")
    ENDIF(CFS_ARCH STREQUAL "X86_64")

    IF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/lib")
      SET (Boost_LIBRARY_DIR ${Boost_INCLUDE_DIR}/lib CACHE STRING
	"The directory containing the Boost libraries." FORCE)
      SET(Boost_FOUND 1)
    ENDIF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/lib")

    IF(CFS_ARCH STREQUAL "X86_64")
      IF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/../lib64")
	SET (Boost_LIBRARY_DIR ${Boost_INCLUDE_DIR}/../lib64 CACHE STRING
	  "The directory containing the Boost libraries." FORCE)
	SET(Boost_FOUND 1)
      ENDIF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/../lib64")
    ENDIF(CFS_ARCH STREQUAL "X86_64")

    IF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/../lib")
      SET (Boost_LIBRARY_DIR ${Boost_INCLUDE_DIR}/../lib CACHE STRING
	"The directory containing the Boost libraries." FORCE)
      SET(Boost_FOUND 1)
    ENDIF(NOT Boost_FOUND AND EXISTS "${Boost_INCLUDE_DIR}/../lib")
  ENDIF(NOT Boost_LIBRARY_DIR)

ENDIF(BOOST_DATE_TIME_LIB_TEST AND Boost_INCLUDE_DIR)



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

  SET(BOOST_LIB_SUFFIX ".a")
  SET(BOOST_LIB_PREFIX "${Boost_LIBRARY_DIR}/lib")
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
IF(CFS_BOOST_VERSION GREATER 1.34)
  SET(BOOST_SYSTEM_LIB
    "${BOOST_LIB_PREFIX}boost_system${BOOST_LIB_SUFFIX}"
    CACHE STRING "BOOST_SYSTEM_LIB")
  MARK_AS_ADVANCED(BOOST_SYSTEM_LIB)
ENDIF(CFS_BOOST_VERSION GREATER 1.34)
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
