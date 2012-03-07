#-------------------------------------------------------------------------------
# Build external CFSDEPS library
#-------------------------------------------------------------------------------
MACRO(BUILD_EXTLIB EXTLIB_NAME SEARCH_FOR_FILE BUILD_PERL_SCRIPT BUILD_LOG_FILE)
  IF(NOT CFS_DEPS_ROOT)
    MESSAGE(FATAL_ERROR "CFS_DEPS_ROOT must be set in order to build ${EXTLIB_NAME}!")
  ENDIF(NOT CFS_DEPS_ROOT)

  IF(NOT EXISTS ${SEARCH_FOR_FILE})
    MESSAGE(STATUS "Building ${EXTLIB_NAME}...")

    EXECUTE_PROCESS(
      COMMAND ${PERL_EXECUTABLE} "${CFS_SOURCE_DIR}/share/scripts/build_tee.pl" ${PERL_EXECUTABLE} ${BUILD_PERL_SCRIPT} "${CFS_TEMP_DIR}/build_vars.pl" "${CFS_TEMP_DIR}/build.log"
      OUTPUT_VARIABLE BUILD_OUT
      ERROR_VARIABLE BUILD_OUT
      RESULT_VARIABLE BUILD_RET)
      
    FILE(WRITE "${CFS_TEMP_DIR}/${BUILD_LOG_FILE}" "${BUILD_OUT}")
    
    IF(BUILD_RET)
      SET(NLINES "12")

      EXECUTE_PROCESS(
	COMMAND "tail" "-${NLINES}" "${CFS_TEMP_DIR}/${BUILD_LOG_FILE}"
	OUTPUT_VARIABLE BUILD_LOG_OUT
	ERROR_VARIABLE BUILD_LOG_OUT
	RESULT_VARIABLE BUILD_RET)

        MESSAGE(FATAL_ERROR "Build of ${EXTLIB_NAME} failed!\n"
	  "Here are the last ${NLINES} lines of the build logfile "
	  "'${CFS_TEMP_DIR}/${BUILD_LOG_FILE}':\n\n${BUILD_LOG_OUT}")
    ENDIF(BUILD_RET)
  ENDIF(NOT EXISTS ${SEARCH_FOR_FILE})
ENDMACRO(BUILD_EXTLIB SEARCH_FOR_FILE BUILD_PERL_SCRIPT BUILD_LOG_FILE)

#-------------------------------------------------------------------------------
# Check if C++ source code runs. Taken from CheckCXXSourceRuns.cmake.
#-------------------------------------------------------------------------------
# - Check if the C++ source code provided in the SOURCE argument compiles and runs.
# CHECK_CXX_SOURCE_RUNS(SOURCE VAR)
#
#  SOURCE - source code to try to compile
#  VAR    - variable to store the result, 1 for success, empty for failure
#
# The following variables may be set before calling this macro to
# modify the way the check is run:
#
#  CMAKE_REQUIRED_FLAGS = string of compile command line flags
#  CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#  CMAKE_REQUIRED_INCLUDES = list of include directories
#  CMAKE_REQUIRED_LIBRARIES = list of libraries to link
#  CFS_REQUIRED_LINK_DIRECTORIES = list of LINK_DIRECTORIES
MACRO(CFS_CHECK_CXX_SOURCE_RUNS SOURCE VAR LINK_DIRS)
  IF("${VAR}" MATCHES "")
    SET(MACRO_CHECK_FUNCTION_DEFINITIONS 
      "-D${VAR} ${CMAKE_REQUIRED_FLAGS}")
    IF(CMAKE_REQUIRED_LIBRARIES)
      SET(CHECK_CXX_SOURCE_COMPILES_ADD_LIBRARIES
        "-DLINK_LIBRARIES:STRING=${CMAKE_REQUIRED_LIBRARIES}")
    ELSE(CMAKE_REQUIRED_LIBRARIES)
      SET(CHECK_CXX_SOURCE_COMPILES_ADD_LIBRARIES)
    ENDIF(CMAKE_REQUIRED_LIBRARIES)
    IF(CMAKE_REQUIRED_INCLUDES)
      SET(CHECK_CXX_SOURCE_COMPILES_ADD_INCLUDES
        "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}")
    ELSE(CMAKE_REQUIRED_INCLUDES)
      SET(CHECK_CXX_SOURCE_COMPILES_ADD_INCLUDES)
    ENDIF(CMAKE_REQUIRED_INCLUDES)
    FILE(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx"
      "${SOURCE}\n")

    MESSAGE(STATUS "Performing Test ${VAR}")
    TRY_RUN(${VAR}_EXITCODE ${VAR}_COMPILED
      ${CMAKE_BINARY_DIR}
      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src.cxx
      COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
      CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS}
      -DCMAKE_SKIP_RPATH:BOOL=${CMAKE_SKIP_RPATH}
      "${CHECK_CXX_SOURCE_COMPILES_ADD_LIBRARIES}"
      "${CHECK_CXX_SOURCE_COMPILES_ADD_INCLUDES}"
      "-DLINK_DIRECTORIES:STRING=${LINK_DIRS}"
      COMPILE_OUTPUT_VARIABLE OUTPUT)
      
    # if it did not compile make the return value fail code of 1
    IF(NOT ${VAR}_COMPILED)
      SET(${VAR}_EXITCODE 1)
    ENDIF(NOT ${VAR}_COMPILED)
    # if the return value was 0 then it worked
    IF("${${VAR}_EXITCODE}" EQUAL 0)
      SET(${VAR} 1 CACHE INTERNAL "Test ${VAR}")
      MESSAGE(STATUS "Performing Test ${VAR} - Success")
      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log 
        "Performing C++ SOURCE FILE Test ${VAR} succeded with the following output:\n"
        "${OUTPUT}\n" 
        "Return value: ${${VAR}}\n"
        "Source file was:\n${SOURCE}\n")
    ELSE("${${VAR}_EXITCODE}" EQUAL 0)
      IF(CMAKE_CROSSCOMPILING AND "${${VAR}_EXITCODE}" MATCHES  "FAILED_TO_RUN")
        SET(${VAR} "${${VAR}_EXITCODE}")
      ELSE(CMAKE_CROSSCOMPILING AND "${${VAR}_EXITCODE}" MATCHES  "FAILED_TO_RUN")
        SET(${VAR} "" CACHE INTERNAL "Test ${VAR}")
      ENDIF(CMAKE_CROSSCOMPILING AND "${${VAR}_EXITCODE}" MATCHES  "FAILED_TO_RUN")

      MESSAGE(STATUS "Performing Test ${VAR} - Failed")
      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log 
        "Performing C++ SOURCE FILE Test ${VAR} failed with the following output:\n"
        "${OUTPUT}\n"  
        "Return value: ${${VAR}_EXITCODE}\n"
        "Source file was:\n${SOURCE}\n")
    ENDIF("${${VAR}_EXITCODE}" EQUAL 0)
  ENDIF("${VAR}" MATCHES "")
ENDMACRO(CFS_CHECK_CXX_SOURCE_RUNS)


#-------------------------------------------------------------------------------
# Import a library that does itself depend on other libraries
#-------------------------------------------------------------------------------
# ADD_DEPENDENT_LIBRARY(NAME LIB DEPENDS)
#
#  NAME    - target name for the added library (targets are added name-1 to name-n) 
#  LIB     - PATH to the library file (.a), may be a list, in that case several targets are added, every target depending on all later targets
#  DEPENDS - general DEPENDENCIES, are added to IMPORTED_LINK_INTERFACE_LIBRARIES
#
FUNCTION(ADD_DEPENDENT_LIBRARY NAME LIB DEPENDS)
  IF(LIB)
    SET(I 0)
    FOREACH(L ${LIB})
      MATH(EXPR I "${I}+1")
      ADD_LIBRARY(${NAME}${I} STATIC IMPORTED)
      SET_PROPERTY(TARGET "${NAME}${I}" PROPERTY IMPORTED_LOCATION "${L}")
    ENDFOREACH(L)

    FOREACH(J RANGE 1 ${I})
      SET(TARGET_LL_DEPENDS "${DEPENDS}")
      FOREACH(K RANGE ${J} ${I})
        SET(TARGET_LL_DEPENDS "${TARGET_LL_DEPENDS}" "${NAME}${K}")
      ENDFOREACH(K)
      SET_PROPERTY(TARGET "${NAME}${J}" PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES "${TARGET_LL_DEPENDS}")
    ENDFOREACH(J)
  ELSE(LIB)
    MESSAGE(SEND_ERROR "ADD_DEPENDENT_LIBRARY called with no LIB")
  ENDIF(LIB)  
ENDFUNCTION(ADD_DEPENDENT_LIBRARY)

