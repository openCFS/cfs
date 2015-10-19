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

  IF(MINGW OR WIN32 OR CYGWIN)
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
      IF(MINGW)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE OUT)
#        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1"
#         ${RESULT} ${${RESULT}})
        string(STRIP "${OUT}" OUT)
        SET(${RESULT} ${OUT})
      ELSE()
        EXECUTE_PROCESS(COMMAND cmd /E:ON /C "${CFS_SOURCE_DIR}/share/scripts/getdate.bat" OUTPUT_VARIABLE OUT)
#        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1"
#	  ${RESULT} ${${RESULT}})
        string(STRIP "${OUT}" OUT)
	SET(${RESULT} ${OUT})
      ENDIF()
    ELSEIF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE OUT)
#        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\3\\2\\1"
#	  ${RESULT} ${${RESULT}})
        string(STRIP "${OUT}" OUT)
        SET(${RESULT} ${OUT})
    ELSE (WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${RESULT} 000000)
    ENDIF (WIN32)
ENDMACRO (TODAY)

#-------------------------------------------------------------------------------
# Convert new line style of files
#-------------------------------------------------------------------------------
MACRO(CHANGE_NEWLINE_STYLE FILES TEMPFILE STYLE)
  foreach(file ${FILES})
    MESSAGE(STATUS "Converting newlines to ${STYLE} for ${file}")
    CONFIGURE_FILE(${file} ${TEMPFILE} @ONLY NEWLINE_STYLE ${STYLE})
    CONFIGURE_FILE(${TEMPFILE} ${file} @ONLY NEWLINE_STYLE ${STYLE})
  endforeach(file)
ENDMACRO(CHANGE_NEWLINE_STYLE)

#-------------------------------------------------------------------------------
# Apply patches
#-------------------------------------------------------------------------------
MACRO (APPLY_PATCHES PATCHES INPUTDIR)
  foreach(patch ${PATCHES})
    MESSAGE(STATUS "Applying patch ${patch}")
#    EXECUTE_PROCESS(
#      COMMAND pwd
#      OUTPUT_QUIET
#      RESULT_VARIABLE RES
#      )
    EXECUTE_PROCESS(
      COMMAND patch -p0 --binary -i "${INPUTDIR}/${patch}"
#      OUTPUT_QUIET
      RESULT_VARIABLE RES
      )

    IF(NOT RES EQUAL 0)
      MESSAGE("A problem occurred while trying to apply patch '${patch}'.")
    ENDIF()
  endforeach(patch)
ENDMACRO (APPLY_PATCHES)


MACRO(DOWNLOAD_CFSDEPS LOCAL_FILE MD5_SUM MIRROR_LIST)
  SET(PERFORM_DOWNLOAD 0)                             
  SET(DOWNLOAD_OKAY 0)                             
  STRING(STRIP ${MD5_SUM} MD5_SUM)                    
  STRING(TOLOWER ${MD5_SUM} MD5_SUM)                  
  SET(TIMEOUT 60)                         

  IF(EXISTS ${LOCAL_FILE})
#    MESSAGE("'${LOCAL_FILE}' already exists!\nComparing MD5 sums...")

    FILE(MD5 ${LOCAL_FILE} ACTUAL_MD5)
    STRING(TOLOWER ${ACTUAL_MD5} ACTUAL_MD5)

#    MESSAGE("ACTUAL_MD5 ${ACTUAL_MD5}")
#    MESSAGE("EXPECTED_MD5 ${MD5_SUM}") 

    STRING(COMPARE EQUAL ${MD5_SUM} ${ACTUAL_MD5} MD5_EQUAL)

    IF(NOT MD5_EQUAL)
      FILE(REMOVE ${LOCAL_FILE})
      SET(PERFORM_DOWNLOAD 1)   
#      MESSAGE("MD5 sums do not match!\nDeleting '${LOCAL_FILE}'...")
    ELSE()                                                          
#      MESSAGE("MD5 sums match! Very fine!")
      SET(DOWNLOAD_OKAY 1)
    ENDIF()
  ELSE()
    SET(PERFORM_DOWNLOAD 1)
  ENDIF()

  IF(PERFORM_DOWNLOAD)
    FOREACH(URL IN ITEMS ${MIRROR_LIST})
      MESSAGE("downloading...
     src='${URL}'
     dst='${LOCAL_FILE}'
     timeout=${TIMEOUT}")

      FILE(DOWNLOAD
        ${URL}
        ${LOCAL_FILE}
        INACTIVITY_TIMEOUT ${TIMEOUT}
        STATUS DL_STATUS
        LOG DL_LOG
        SHOW_PROGRESS)

      LIST(GET DL_STATUS 0 DL_FAIL)
      LIST(GET DL_STATUS 1 DL_MSG)

      IF(NOT DL_FAIL)
        FILE(MD5 ${LOCAL_FILE} ACTUAL_MD5)
        STRING(TOLOWER ${ACTUAL_MD5} ACTUAL_MD5)

#        MESSAGE("ACTUAL_MD5 ${ACTUAL_MD5}")
#        MESSAGE("EXPECTED_MD5 ${MD5_SUM}")

        STRING(COMPARE EQUAL ${MD5_SUM} ${ACTUAL_MD5} MD5_EQUAL)

        IF(MD5_EQUAL)
          SET(DOWNLOAD_OKAY 1)
          BREAK()
        ELSE()
          FILE(REMOVE ${LOCAL_FILE})
        ENDIF()

      ELSE()
        MESSAGE("Download failed: ${DL_MSG}")
        FILE(REMOVE ${LOCAL_FILE})
      ENDIF()

    ENDFOREACH()
  ENDIF()

  IF(NOT DOWNLOAD_OKAY)
    MESSAGE(FATAL_ERROR "Download failed!")
  ENDIF()
ENDMACRO()

#-------------------------------------------------------------------------------
# Extract ZIP_FILE to TARGET_DIR
#-------------------------------------------------------------------------------
MACRO(ZIP_FROM_CACHE ZIP_FILE TARGET_DIR)
  IF(EXISTS "${ZIP_FILE}")
    MESSAGE("Found precompiled version ${ZIP_FILE}.")
    EXECUTE_PROCESS(
      COMMAND "${CMAKE_COMMAND}" -E tar xzf "${ZIP_FILE}"
      WORKING_DIRECTORY "${TARGET_DIR}"
      OUTPUT_QUIET
      RESULT_VARIABLE rv
      )
    IF(NOT "${rv}" STREQUAL "0")
      MESSAGE(SEND_ERROR "Could not extract ${ZIP_FILE}.")
    ENDIF()
  ELSE()
    MESSAGE(SEND_ERROR "Could not find precompiled ${ZIP_FILE}.")
  ENDIF()
  
ENDMACRO()

#-------------------------------------------------------------------------------
# Create ZIP_FILE
# If a TMP_DIR/src/*-build/install_manifest.txt exists, we zip all files
# listed in there else we zip TMP_INSTALL_DIR.
#-------------------------------------------------------------------------------
MACRO(ZIP_TO_CACHE ZIP_FILE TMP_DIR)
  STRING(REGEX REPLACE "^.+[/\\]" "" ZIP_NAME ${ZIP_FILE})
  STRING(REGEX REPLACE "${ZIP_NAME}$" "" TARGET_DIR ${ZIP_FILE})
  
  IF(NOT EXISTS "${TARGET_DIR}")
    FILE(MAKE_DIRECTORY "${TARGET_DIR}")
  ENDIF()
  
  FILE(GLOB MANIFESTS "${TMP_DIR}/src/*-build/install_manifest.txt")
  IF("${MANIFESTS}" STREQUAL "")
    # No manifests exists -> zip TMP_DIR
    EXECUTE_PROCESS(
      COMMAND zip -g -r ${ZIP_FILE} "."
      WORKING_DIRECTORY "${TMP_DIR}"
      RESULT_VARIABLE rv
    )
  ELSE()
    # Manifests exists -> zip files listed therein
    FOREACH(manifest ${MANIFESTS})
      EXECUTE_PROCESS(
        COMMAND sed "s@${CMAKE_CURRENT_BINARY_DIR}/@@g" ${manifest} 
        COMMAND zip -@ -g ${ZIP_FILE}
        OUTPUT_QUIET
        RESULT_VARIABLE rv
      )
    ENDFOREACH()
  ENDIF()
  IF(NOT "${rv}" STREQUAL "0")
    MESSAGE(WARNING "Could not create ${ZIP_NAME} at ${TARGET_DIR}.")
  ENDIF()
ENDMACRO()

# ------------------------------------------------------------------------------
# Generate a package name for the precompiled zip file. 
# Names the compiler. If C/C++ and Fortan are different, both are named. Otherwise it would
# lead to issues with clang
# ------------------------------------------------------------------------------
MACRO(PRECOMPILED_ZIP RETVAL IN_PACKAGE_NAME IN_PACKAGE_VER)
  # in the legacy cfs there was for WIN32 ${CMAKE_BUILD_TYPE} instead of the compiler stuff
  IF(${CMAKE_CXX_COMPILER_VERSION} STREQUAL ${CMAKE_Fortran_COMPILER_VERSION})
    SET(${RETVAL} "${CFS_DEPS_CACHE_DIR}/precompiled/${IN_PACKAGE_NAME}_${IN_PACKAGE_VER}_${CFS_ARCH_STR}_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}_${CMAKE_BUILD_TYPE}.zip")
  ELSE()
    # in the intel case:
    # FC_VERSION=16.0 20150815
    # CMAKE_Fortran_COMPILER_VERSION=16.0.0.20150815
    SET(${RETVAL} "${CFS_DEPS_CACHE_DIR}/precompiled/${IN_PACKAGE_NAME}_${IN_PACKAGE_VER}_${CFS_ARCH_STR}_C_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}_F_${CMAKE_Fortran_COMPILER_ID}_${CMAKE_Fortran_COMPILER_VERSION}_${CMAKE_BUILD_TYPE}.zip")  
  ENDIF()    
ENDMACRO()

# don't add Release or Debug when the package is built independently
MACRO(PRECOMPILED_ZIP_NOBUILD RETVAL IN_PACKAGE_NAME IN_PACKAGE_VER)
  IF(${CMAKE_CXX_COMPILER_VERSION} STREQUAL ${CMAKE_Fortran_COMPILER_VERSION})
    SET(${RETVAL} "${CFS_DEPS_CACHE_DIR}/precompiled/${IN_PACKAGE_NAME}_${IN_PACKAGE_VER}_${CFS_ARCH_STR}_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}.zip")
  ELSE()
    SET(${RETVAL} "${CFS_DEPS_CACHE_DIR}/precompiled/${IN_PACKAGE_NAME}_${IN_PACKAGE_VER}_${CFS_ARCH_STR}_C_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}_F_${CMAKE_Fortran_COMPILER_ID}_${CMAKE_Fortran_COMPILER_VERSION}.zip")  
  ENDIF()    
ENDMACRO()


#------------------------------------------------------
# Display all available variables
#------------------------------------------------------
MACRO(DISPLAY_ALL_VARIABLES)
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    message("${_variableName}=${${_variableName}}")
  endforeach()
ENDMACRO()