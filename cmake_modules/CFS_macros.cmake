#
# A big shout out to the cmake gurus @ compiz
#
function (colormsg)
  string (ASCII 27 _escape)
  set(WHITE "29")
  set(GRAY "30")
  set(RED "31")
  set(GREEN "32")
  set(YELLOW "33")
  set(BLUE "34")
  set(MAG "35")
  set(CYAN "36")

  foreach (color WHITE GRAY RED GREEN YELLOW BLUE MAG CYAN)
    set(HI${color} "1\;${${color}}")
    set(LO${color} "2\;${${color}}")
    set(_${color}_ "4\;${${color}}")
    set(_HI${color}_ "1\;4\;${${color}}")
    set(_LO${color}_ "2\;4\;${${color}}")
  endforeach()

  execute_process(COMMAND echo "echo $PPID"
                  COMMAND sh
                  COMMAND xargs ps -p RESULT_VARIABLE RES OUTPUT_VARIABLE OUT ERROR_VARIABLE ERR)
#  MESSAGE("RES ${RES}")
#  MESSAGE("OUT ${OUT}")
#  MESSAGE("ERR ${ERR}")

  SET(DISABLE_COLOR 0)
  # Do not display color in ccmake or cmake-gui
  if(OUT MATCHES "ccmake")
    SET(DISABLE_COLOR 1)
  endif()
  if(OUT MATCHES "cmake-gui")
    SET(DISABLE_COLOR 1)
  endif()

  set(str "")
  set(coloron FALSE)
  foreach(arg ${ARGV})
    if (NOT ${${arg}} STREQUAL "")
      if (CMAKE_COLOR_MAKEFILE AND NOT DISABLE_COLOR)
        set(str "${str}${_escape}[${${arg}}m")
        set(coloron TRUE)
      endif()
    else()
      set(str "${str}${arg}")
      if (coloron AND NOT DISABLE_COLOR)
        set(str "${str}${_escape}[0m")
        set(coloron FALSE)
      endif()
      set(str "${str} ")
    endif()
  endforeach()
  message(STATUS ${str})
endfunction()

# colormsg("Colors:"  
#   WHITE "white" GRAY "gray" GREEN "green" 
#   RED "red" YELLOW "yellow" BLUE "blue" MAG "mag" CYAN "cyan" 
#   _WHITE_ "white" _GRAY_ "gray" _GREEN_ "green" 
#   _RED_ "red" _YELLOW_ "yellow" _BLUE_ "blue" _MAG_ "mag" _CYAN_ "cyan" 
#   _HIWHITE_ "white" _HIGRAY_ "gray" _HIGREEN_ "green" 
#   _HIRED_ "red" _HIYELLOW_ "yellow" _HIBLUE_ "blue" _HIMAG_ "mag" _HICYAN_ "cyan" 
#   HIWHITE "white" HIGRAY "gray" HIGREEN "green" 
#   HIRED "red" HIYELLOW "yellow" HIBLUE "blue" HIMAG "mag" HICYAN "cyan" 
#   "right?")

#-------------------------------------------------------------------------------
# Build external CFSDEPS library
#-------------------------------------------------------------------------------
MACRO(BUILD_EXTLIB EXTLIB_NAME SEARCH_FOR_FILE BUILD_PERL_SCRIPT BUILD_LOG_FILE)
  IF(NOT CFS_DEPS_ROOT)
    MESSAGE(FATAL_ERROR "CFS_DEPS_ROOT must be set in order to build ${EXTLIB_NAME}!")
  ENDIF(NOT CFS_DEPS_ROOT)

  IF(NOT EXISTS ${SEARCH_FOR_FILE})
    colormsg(HIBLUE "Building ${EXTLIB_NAME}...")
    
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
# Get the current date
#-------------------------------------------------------------------------------
MACRO (TODAY RESULT)
    IF (WIN32)
        EXECUTE_PROCESS(COMMAND "date" "/T" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1"
${RESULT} ${${RESULT}})
    ELSEIF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1"
${RESULT} ${${RESULT}})
    ELSE (WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${RESULT} 000000)
    ENDIF (WIN32)
ENDMACRO (TODAY)

#-------------------------------------------------------------------------------
# Convert new line style of files
#-------------------------------------------------------------------------------
MACRO (CHANGE_NEWLINE_STYLE FILES TEMPFILE STYLE)
  foreach(file ${FILES})
    MESSAGE("Converting newlines to ${STYLE} for ${file}")
    CONFIGURE_FILE(${file} ${TEMPFILE} @ONLY NEWLINE_STYLE ${STYLE})
    CONFIGURE_FILE(${TEMPFILE} ${file} @ONLY NEWLINE_STYLE ${STYLE})
  endforeach(file)
ENDMACRO (CHANGE_NEWLINE_STYLE)

#-------------------------------------------------------------------------------
# Apply patches
#-------------------------------------------------------------------------------
MACRO (APPLY_PATCHES PATCHES INPUTDIR)
  foreach(patch ${PATCHES})
    MESSAGE("Applying patch ${patch}")
    EXECUTE_PROCESS(
      COMMAND patch -p0 -i "${INPUTDIR}/${patch}" OUTPUT_QUIET
      )
  endforeach(patch)
ENDMACRO (APPLY_PATCHES)
