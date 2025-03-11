include(CMakePrintHelpers)
# easy check if already loaded
set(CFS_macros "1")
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

  IF(WIN32)
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

    #MESSAGE(STATUS "Performing Test ${VAR}")
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
      #MESSAGE(STATUS "Performing Test ${VAR} - Success")
      FILE(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log 
        "Performing C++ SOURCE FILE Test ${VAR} succeded with the following output:\n"
        "${OUTPUT}\n" 
        "Return value: ${${VAR}}\n"
        "Source file was:\n${SOURCE}\n")
    ELSE("${${VAR}_EXITCODE}" EQUAL 0)
      SET(${VAR} "" CACHE INTERNAL "Test ${VAR}")
 
       #MESSAGE(STATUS "Performing Test ${VAR} - Failed")
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
MACRO(TODAY RESULT)
  IF(WIN32)
    EXECUTE_PROCESS(COMMAND cmd /E:ON /C "${CFS_SOURCE_DIR}/share/scripts/getdate.bat" OUTPUT_VARIABLE OUT)
    string(STRIP "${OUT}" OUT)
    SET(${RESULT} ${OUT})
  ELSEIF(UNIX)
    EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE OUT)
    string(STRIP "${OUT}" OUT)
    SET(${RESULT} ${OUT})
  ELSE()
    MESSAGE(SEND_ERROR "date not implemented")
    SET(${RESULT} 000000)
  ENDIF()
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
#      RESULT_VARIABLE RES)

    EXECUTE_PROCESS(
      COMMAND patch -p0 --binary -i "${INPUTDIR}/${patch}"
#      OUTPUT_QUIET
      RESULT_VARIABLE RES)


    IF(NOT RES EQUAL 0)
      MESSAGE("A problem occurred while trying to apply patch '${patch}'.")
    ENDIF()
  endforeach(patch)
ENDMACRO (APPLY_PATCHES)


# dowloads a file if the local file does not already exist
MACRO(DOWNLOAD_CFSDEPS LOCAL_FILE MD5_SUM MIRROR_LIST)
  SET(PERFORM_DOWNLOAD 0)                             
  SET(DOWNLOAD_OKAY 0)                             
  STRING(STRIP ${MD5_SUM} MD5_SUM)                    
  STRING(TOLOWER ${MD5_SUM} MD5_SUM)                  
  SET(TIMEOUT 120)                         

  IF(EXISTS ${LOCAL_FILE})
    MESSAGE("'${LOCAL_FILE}' already exists!\nComparing MD5 sums...")

    FILE(MD5 ${LOCAL_FILE} ACTUAL_MD5)
    STRING(TOLOWER ${ACTUAL_MD5} ACTUAL_MD5)

    # MESSAGE("ACTUAL_MD5 ${ACTUAL_MD5}")
    # MESSAGE("EXPECTED_MD5 ${MD5_SUM}") 

    STRING(COMPARE EQUAL ${MD5_SUM} ${ACTUAL_MD5} MD5_EQUAL)

    IF(NOT MD5_EQUAL)
      FILE(REMOVE ${LOCAL_FILE})
      SET(PERFORM_DOWNLOAD 1)   
      MESSAGE("MD5 sums do not match!\nDeleting '${LOCAL_FILE}'...")
    ELSE()                                                          
      MESSAGE("MD5 sums match! Very fine!")
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

      IF(EXISTS ${URL})
        GET_FILENAME_COMPONENT(FOLDER ${LOCAL_FILE} DIRECTORY)
        MESSAGE("src is an existing local file: copy ${URL} -> ${LOCAL_FILE}")
        CONFIGURE_FILE("${URL}" "${LOCAL_FILE}" COPYONLY)
      ELSE()
        FILE(DOWNLOAD
          ${URL}
          ${LOCAL_FILE}
          INACTIVITY_TIMEOUT ${TIMEOUT}
          STATUS DL_STATUS
          LOG DL_LOG
          SHOW_PROGRESS
          TLS_VERIFY OFF)

        LIST(GET DL_STATUS 0 DL_FAIL)
        LIST(GET DL_STATUS 1 DL_MSG)
      ENDIF()

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
          message("md5sum ${ACTUAL_MD5} != expected ${MD5_SUM} -> remove local file and continue")
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
    MESSAGE("Found precompiled version ${ZIP_FILE}")
    EXECUTE_PROCESS(
      COMMAND "${CMAKE_COMMAND}" -E tar xvfz "${ZIP_FILE}"
      WORKING_DIRECTORY "${TARGET_DIR}"
      OUTPUT_VARIABLE ov
      ERROR_VARIABLE ev
      RESULT_VARIABLE rv
      )
    IF(NOT "${rv}" STREQUAL "0")
      MESSAGE(SEND_ERROR "Could not extract ${ZIP_FILE}:\nov=${out}\nev=${ev}")
    ENDIF()
  ELSE()
    MESSAGE(SEND_ERROR "Could not find precompiled ${ZIP_FILE}")
  ENDIF()
ENDMACRO()

#-------------------------------------------------------------------------------
# Create ZIP_FILE
# If a TMP_DIR/src/*-build/install_manifest.txt exists (when built with cmake), we zip all files listed in there.
# Else we try to be smart on TMP_DIR (move libs from lib64 to lib), make a zip out of it and copy the content to CMAKE_CURRENT_BINARY_DIR
#-------------------------------------------------------------------------------
MACRO(ZIP_TO_CACHE ZIP_FILE TMP_DIR)
  STRING(REGEX REPLACE "^.+[/\\]" "" ZIP_NAME ${ZIP_FILE})
  get_filename_component(TARGET_DIR ${ZIP_FILE} DIRECTORY)
  IF(NOT EXISTS "${TARGET_DIR}")
    message(STATUS "creating directroy: ZIP_FILE=${ZIP_FILE} -> ZIP_NAME=${ZIP_NAME}, TARGET_DIR=${TARGET_DIR}")
    FILE(MAKE_DIRECTORY "${TARGET_DIR}")
  ENDIF()
  
  # Fabian: was src/*-build/install_manifest.tx before but that fails for openblas
  FILE(GLOB MANIFESTS "${TMP_DIR}/src/*/install_manifest.txt")

  message("ZIP_TO_CACHE ZF=${ZIP_FILE} TD=${TMP_DIR} M=${MANIFESTS}")

  IF("${MANIFESTS}" STREQUAL "")
    message("ZIP_TO_CACHE: zip ${TMP_DIR} as manifest ${TMP_DIR}/src/*-build/install_manifest.txt was not found")
    # No manifests exists -> zip TMP_DIR

    # standard make or configure does not necessarily know about lib
    IF(NOT EXISTS "${TMP_DIR}/lib")
      FILE(MAKE_DIRECTORY "${TMP_DIR}/lib")
    ENDIF()

    # move any lib64 to lib. Extend to lib/files if necessary
    IF(EXISTS "${TMP_DIR}/lib64")
      message(STATUS "move files from ${TMP_DIR}/lib64 to lib")
	    FILE(COPY "${TMP_DIR}/lib64/" DESTINATION "${TMP_DIR}/lib")
      FILE(REMOVE_RECURSE "${TMP_DIR}/lib64")
    ENDIF()

    file(GLOB_RECURSE files_to_zip RELATIVE "${TMP_DIR}" "${TMP_DIR}/*")
    # message(STATUS "found ${files_to_zip}")
    string(REPLACE ";" "\n" files_to_zip "${files_to_zip}")
    # message(STATUS "packing ${files_to_zip}")
    file(WRITE "${TMP_DIR}/files_to_zip.txt" ${files_to_zip})
    # create the zip
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar cvf ${ZIP_FILE} --files-from=${TMP_DIR}/files_to_zip.txt --format=zip
      WORKING_DIRECTORY "${TMP_DIR}"
      RESULT_VARIABLE rv
      OUTPUT_VARIABLE ov
      ERROR_VARIABLE ev
    )
    # remove the file again
    file(REMOVE "${TMP_DIR}/files_to_zip.txt")
      
    # copy data from TMP_DIR install dir to the target where cfs needs it and where the zip would be extracted.
    # Install only the contents of TMP_DIR since trying to copy . too leads to permission error
    # file INSTALL cannot set permissions on CMAKE_CURRENT_BINARY_DIR
    # if the build user does not own the target directory (but is allowed to write into it)
    file(GLOB files LIST_DIRECTORIES true RELATIVE "${TMP_DIR}" "${TMP_DIR}/*")
    FOREACH(file ${files})
      # message(STATUS "installing contents of ${file}")
      FILE(INSTALL "${TMP_DIR}/${file}" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
    ENDFOREACH()
  ELSE()
    # Manifests exists -> zip files listed therein
    file(WRITE "${TMP_DIR}/files_to_zip.txt" "")
    FOREACH(manifest ${MANIFESTS})
      file(READ ${manifest} manifest_string)
      # remove the leading path, i.e. make relative
      string(REPLACE "${CMAKE_CURRENT_BINARY_DIR}/" "" manifest_string "${manifest_string}")
      # change to semicolon-seperated list
      string(REPLACE "\n" ";" manifest_string "${manifest_string}")
      foreach(file ${manifest_string})
        #message("file=${file}")
        file(APPEND "${TMP_DIR}/files_to_zip.txt" "${file}\n")
        #message("files_to_zip=${files_to_zip}")
      endforeach()
    ENDFOREACH()
    #message("files_to_zip=${files_to_zip}")
    EXECUTE_PROCESS(
      COMMAND ${CMAKE_COMMAND} -E tar "cf" ${ZIP_FILE} --files-from=${TMP_DIR}/files_to_zip.txt --format=zip
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      # OUTPUT_QUIET
      RESULT_VARIABLE rv
      OUTPUT_VARIABLE ov
      ERROR_VARIABLE ev
    )
  ENDIF()
  # ckeck if it went well
  if(NOT "${rv}" STREQUAL "0")
    MESSAGE(WARNING "Could not create ${ZIP_NAME} at ${TARGET_DIR}:\nout=${ov}\nerr=${ev}")
  endif()
 
ENDMACRO()

# ------------------------------------------------------------------------------
# Generate a package name for the precompiled zip file. 
# Names the compiler. If C/C++ and Fortan are different, both are named. Otherwise it would
# lead to issues with clang which needs either gfortran or ifort as companion
# Also check PRECOMPILED_ZIP_NOBUILD, it might be more appropriate
# ------------------------------------------------------------------------------
MACRO(PRECOMPILED_ZIP RETVAL IN_PACKAGE_NAME IN_PACKAGE_VER)
  # there is complex issue with ifort. On Ubuntu there was the case CMAKE_Fortran_COMPILER_VERSION is not set?!

  # use our own variable
  IF(DEFINED CMAKE_Fortran_COMPILER_VERSION AND NOT "${CMAKE_Fortran_COMPILER_VERSION}" STREQUAL "")
    SET(Fortran_COMPILER_VERSION "${CMAKE_Fortran_COMPILER_VERSION}")
  ELSE()
    message(WARNING "Fortran compiler version not recognised (by CMake)")
    SET(Fortran_COMPILER_VERSION "?")
  ENDIF()

  IF(${CMAKE_CXX_COMPILER_VERSION} STREQUAL ${Fortran_COMPILER_VERSION})
    SET(${RETVAL} "${CFS_DEPS_CACHE_DIR}/precompiled/${IN_PACKAGE_NAME}_${IN_PACKAGE_VER}_${CFS_ARCH_STR}_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}_${CMAKE_BUILD_TYPE}-lib.zip")
  ELSE()
    SET(${RETVAL} "${CFS_DEPS_CACHE_DIR}/precompiled/${IN_PACKAGE_NAME}_${IN_PACKAGE_VER}_${CFS_ARCH_STR}_C_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}_F_${CMAKE_Fortran_COMPILER_ID}_${Fortran_COMPILER_VERSION}_${CMAKE_BUILD_TYPE}-lib.zip")
  ENDIF()    
  
ENDMACRO()

# don't add Release or Debug when the package is built independently
MACRO(PRECOMPILED_ZIP_NOBUILD RETVAL IN_PACKAGE_NAME IN_PACKAGE_VER)
  # see PRECOMPILED_ZIP
  IF(DEFINED CMAKE_Fortran_COMPILER_VERSION AND NOT "${CMAKE_Fortran_COMPILER_VERSION}" STREQUAL "")
    SET(Fortran_COMPILER_VERSION "${CMAKE_Fortran_COMPILER_VERSION}")
  ELSE()
    SET(Fortran_COMPILER_VERSION "?")
  ENDIF()

  IF(${CMAKE_CXX_COMPILER_VERSION} STREQUAL ${Fortran_COMPILER_VERSION})
    SET(${RETVAL} "${CFS_DEPS_CACHE_DIR}/precompiled/${IN_PACKAGE_NAME}_${IN_PACKAGE_VER}_${CFS_ARCH_STR}_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}-lib.zip")
  ELSE()
    SET(${RETVAL} "${CFS_DEPS_CACHE_DIR}/precompiled/${IN_PACKAGE_NAME}_${IN_PACKAGE_VER}_${CFS_ARCH_STR}_C_${CMAKE_CXX_COMPILER_ID}_${CMAKE_CXX_COMPILER_VERSION}_F_${CMAKE_Fortran_COMPILER_ID}_${Fortran_COMPILER_VERSION}-lib.zip")
  ENDIF()    

ENDMACRO()

# overwrite a cmake variable by the content of a environment variable with same name
macro(set_from_env VAR)
  if(DEFINED ENV{${VAR}})
    set(${VAR} $ENV{${VAR}})
  endif()
endmacro()    

# display all available variables. 
# You can add optionally add a match key. E.g. dump_variables("OpenMP")
# from  https://stackoverflow.com/questions/9298278/cmake-print-out-all-accessible-variables-in-a-script
function(dump_variables)
  get_cmake_property(_variableNames VARIABLES)
  list (SORT _variableNames)
  foreach (_variableName ${_variableNames})
    if ((NOT DEFINED ARGV0) OR _variableName MATCHES ${ARGV0})
      message(STATUS "${_variableName}=${${_variableName}}")
    endif()
  endforeach()
endfunction()


# dump the content of the given directory
macro(DUMP_DIR DIR)
  message("DUMP_DIR ${DIR}:")
  file(GLOB_RECURSE _DIR_VAR FOLLOW_SYMLINKS ${DIR})
  foreach (_NAME ${_DIR_VAR})
    message("${_NAME}")
  endforeach()  
endmacro()


# This is a simple assert to check if a cmake variable is set.
# use like assert_set(CTEST_SCRIPTS_DIR) or assert_set("CTEST_SCRIPTS_DIR")
macro(assert_set test)
  if(NOT ${test})
     message(FATAL_ERROR "assertion failed, the variable ${test} is not set")
  endif()
endmacro()

# This simple assert checks if a string or similar is not "")
macro(assert_not_empty test)
  if(${test} STREQUAL "")
    message(FATAL_ERROR "assert failed, string is empty")
  endif()
endmacro()

# This is a simple assert to check if a cmake variable is unset.
# see assert_set() 
macro(assert_unset test)
  if(${test})
     message(FATAL_ERROR  "assertion failed: variable ${test} is already set")
  endif()
endmacro()

# Assert a value is different from a required one
macro(assert_not test compare)
  if(${test} STREQUAL ${compare})
    message(FATAL_ERROR "assert failed, string is ${test} but shall be different")
  endif()
endmacro()

# assert that a file exists
macro(assert_file_exists test)
  if(NOT EXISTS ${test})
     message(FATAL_ERROR "file ${test} not found")
  endif()
endmacro()      

# assert that a directorx exists
macro(assert_dir_exists test)
  if(NOT IS_DIRECTORY ${test})
     message(FATAL_ERROR "directory ${test} not found")
  endif()
endmacro()      


#-------------------------------------------------------------------------------
# Run a program with optional logging and throw error if return code is non-zero
#-------------------------------------------------------------------------------
function(run_program workdir logfile cmd)
  if(NOT IS_DIRECTORY "${workdir}")
    MESSAGE(FATAL_ERROR "Working directory '${workdir}' does not exist or is not a directory.")
  endif(NOT IS_DIRECTORY "${workdir}")

  if("${logfile}" STREQUAL "")
  
    execute_process(COMMAND ${cmd} ${ARGN}
                    WORKING_DIRECTORY ${workdir}
                    RESULT_VARIABLE rv)
    if(NOT ${rv} EQUAL 0)
      message(FATAL_ERROR "Command '${cmd} ${ARGN}' failed with return code ${rv}.")
    endif(NOT ${rv} EQUAL 0)
                    
  else("${logfile}" STREQUAL "")
  
    execute_process(COMMAND ${cmd} ${ARGN}
                    WORKING_DIRECTORY ${workdir}
                    RESULT_VARIABLE rv
                    OUTPUT_VARIABLE outerr
                    ERROR_VARIABLE outerr)
    file(WRITE "${logfile}" "${outerr}")
    if(NOT ${rv} EQUAL 0)
      message(FATAL_ERROR "Command ${cmd} ${ARGN} failed with return code ${rv}.
        See ${logfile} for details.")
    endif(NOT ${rv} EQUAL 0)
    
  endif("${logfile}" STREQUAL "")

endfunction(run_program)
