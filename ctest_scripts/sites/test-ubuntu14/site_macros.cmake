# define some common macros for this site
# this reduces excessive copy & paste

macro(INIT_TEST)
  # this sets the variables
  #   NIGHTLY_TESTSUITE_DIR
  #   NIGHTLY_CACHE_DIR
  #   NIGHTLY_SOURCE_DIR
  #   NIGHTLY_BUILD_DIR
  include("${CMAKE_CURRENT_LIST_DIR}/site_specific_vars.cmake")

  # Set source & binary directories on the machine
  # We take the name of this file and append it to NIGHTLY_BUILD_DIR
  SET(CTEST_SOURCE_DIRECTORY "${NIGHTLY_SOURCE_DIR}")
  get_filename_component(FNAME "${CMAKE_CURRENT_LIST_FILE}" NAME_WE)
  set(CTEST_BINARY_DIRECTORY "${NIGHTLY_BUILD_DIR}/${FNAME}")

  # Include macros for CFS++ specific testing (e.g. default initial cache, etc)
  INCLUDE("${CTEST_SOURCE_DIRECTORY}/ctest_scripts/shared/test_macros.cmake")
  # Identify distro.
  IDENTIFY_DISTRO()
  # empty binary directory
  #ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY}) # this fails if there is no "CMakeLists.txt"
  file(REMOVE_RECURSE "${CTEST_BINARY_DIRECTORY}")
  file(MAKE_DIRECTORY "${CTEST_BINARY_DIRECTORY}")

  # write ctest config
  WRITE_CTEST_CONFIG()

  # set ctest & cmake commands to the executables corresponding to the one used to run this script.
  SET(CTEST_COMMAND  "\"${CTEST_EXECUTABLE_NAME}\"")
  SET(CTEST_CMAKE_COMMAND  "\"${CMAKE_EXECUTABLE_NAME}\"")
endmacro()

macro(DO_TESTING)
  #-----------------------------------------------------------------------------
  # Start the Update, Configure, Build and Testing
  #-----------------------------------------------------------------------------
  message("Start dashboard...")
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  # see http://www.vtk.org/Wiki/CTest:Nightly,_Experimental,_Continuous
  ctest_start(Nightly)

  message("  Update")
  ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)

  message("  Configure")
  file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" ${CTEST_INITIAL_CACHE})
  # some configurations require to repeat the same ctest_configure() command
  ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
  ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

  message("  Build")
  ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

  message("  Test")
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

  message("  Submit")
  ctest_submit(RETURN_VALUE res)
  message("  All done")	
endmacro()


macro(SET_COMPILER_ENV COMPILER_TYPE)
  message("Setting environment for ${COMPILER_TYPE}")
  # =========================================================================
  #
  #  Set evironment for different compilers (gcc and intel)
  #
  # =========================================================================
  if(${COMPILER_TYPE} STREQUAL "GCC")
    #-----------------------------------------------------------------------------
    # Set the following environment variables for the test run. This can be used
    # to specifiy the compilers and that all messages should be output in English
    # language, so that CTest may properly parse them.
    #-----------------------------------------------------------------------------
    SET(ENV{CC} "/usr/bin/gcc")
    SET(ENV{CXX} "/usr/bin/g++")
    SET(ENV{FC} "/usr/bin/gfortran")
    SET(ENV{LC_MESSAGES} "C")
    SET(ENV{LC_ALL} "C")
    SET(ENV{LANG} "C")
    SET(ENV{LANGUAGE} "C")
    # determine compiler info
    SET(IDCOMP_TEMPL "${CTEST_SOURCE_DIRECTORY}/share/scripts/identify_compiler.cmake.in")
    SET(COMPILER_ID_FILE "${CTEST_BINARY_DIRECTORY}/CMakeFiles/out.cmake")
    SET(IDENTIFY_COMPILER_SRC "${CTEST_SOURCE_DIRECTORY}/share/scripts/IdentifyCXXCompiler.cpp")
    SET(COMPILER "$ENV{CXX}")
    SET(ID_CXX "${CTEST_BINARY_DIRECTORY}/share/scripts/identify_cxx.cmake")
    CONFIGURE_FILE("${IDCOMP_TEMPL}" "${ID_CXX}" @ONLY)
    # create build/tmp directory
    EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -E make_directory "${CTEST_BINARY_DIRECTORY}/tmp" WORKING_DIRECTORY "${CTEST_BINARY_DIRECTORY}" RESULT_VARIABLE RETVAL)
    # read ID_CXX
    EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -P "${ID_CXX}" WORKING_DIRECTORY "${CTEST_BINARY_DIRECTORY}/tmp" RESULT_VARIABLE RETVAL)
    INCLUDE("${COMPILER_ID_FILE}")

  elseif(${COMPILER_TYPE} STREQUAL "ICC")

    if(NOT INTEL_COMPILER_PATH)
      set(INTEL_COMPILER_PATH "/share/programs/intel/composer_xe_2015.2.164")
      message("INTEL_COMPILER_PATH not defined, guessing ${INTEL_COMPILER_PATH}")
    endif(NOT INTEL_COMPILER_PATH)

    SET(INTEL_COMPVARS_SH "${CTEST_BINARY_DIRECTORY}/CMakeFiles/out.sh")
    SET(INTEL_COMPVARS_CMAKE "${CTEST_BINARY_DIRECTORY}/CMakeFiles/out.cmake")

    FILE(WRITE "${INTEL_COMPVARS_SH}"
    "
    source ${INTEL_COMPILER_PATH}/bin/compilervars.sh intel64
    ${CTEST_CMAKE_COMMAND} -E environment
    ")

    SET(SHCMD
      bash
      "${INTEL_COMPVARS_SH}"
    )
    EXECUTE_PROCESS(
      COMMAND ${SHCMD}
      OUTPUT_VARIABLE INTEL_ENV_SH
    )
    STRING(REPLACE "\n" ";" INTEL_ENV_SH ${INTEL_ENV_SH})

    SET(INTEL_ENV_CMAKE "")
    FOREACH(LINE IN ITEMS ${INTEL_ENV_SH})
      STRING(REPLACE "=" ";" LINE_TOKENS ${LINE})
      LIST(GET LINE_TOKENS 0 VAR_NAME)
      LIST(GET LINE_TOKENS 1 VAR_VALUE)
      IF(NOT VAR_NAME STREQUAL "INTEL_LICENSE_FILE" AND
         NOT VAR_NAME STREQUAL "PWD" AND
         NOT VAR_NAME STREQUAL "SHLVL" AND
         NOT VAR_NAME STREQUAL "_")
        SET(INTEL_ENV_CMAKE "${INTEL_ENV_CMAKE}\nSET(ENV{${VAR_NAME}} \"${VAR_VALUE}\")")
      ENDIF()
    ENDFOREACH()

    #message("INTEL_ENV_CMAKE XXX ${INTEL_ENV_CMAKE} XXX INTEL_ENV_CMAKE")
    FILE(WRITE "${INTEL_COMPVARS_CMAKE}" "${INTEL_ENV_CMAKE}")
    INCLUDE("${INTEL_COMPVARS_CMAKE}")

    SET(ENV{LM_LICENSE_FILE} "/share/programs/intel/licenses/pmklicserv.lic")
    SET(ENV{INTEL_LICENSE_FILE} "$ENV{LM_LICENSE_FILE}")
    SET(ENV{CC} "icc")
    SET(ENV{CXX} "icpc")
    SET(ENV{FC} "ifort")
    SET(ENV{LC_MESSAGES} "C")
    SET(ENV{LC_ALL} "C")
    SET(ENV{LANG} "C")
    SET(ENV{LANGUAGE} "C")

    SET(IDCOMP_TEMPL "${CTEST_SOURCE_DIRECTORY}/share/scripts/identify_compiler.cmake.in")
    SET(COMPILER_ID_FILE "${CTEST_BINARY_DIRECTORY}/CMakeFiles/out.cmake")
    SET(IDENTIFY_COMPILER_SRC "${CTEST_SOURCE_DIRECTORY}/share/scripts/IdentifyCXXCompiler.cpp")

    # determine compiler info
    SET(COMPILER "$ENV{CXX}")
    SET(ID_CXX "${CTEST_BINARY_DIRECTORY}/share/scripts/identify_cxx.cmake")
    CONFIGURE_FILE("${IDCOMP_TEMPL}" "${ID_CXX}" @ONLY)
    # create build/tmp directory
    EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -E make_directory "${CTEST_BINARY_DIRECTORY}/tmp" WORKING_DIRECTORY "${CTEST_BINARY_DIRECTORY}" RESULT_VARIABLE RETVAL)
    # read ID_CXX
    EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -P "${ID_CXX}" WORKING_DIRECTORY "${CTEST_BINARY_DIRECTORY}/tmp" RESULT_VARIABLE RETVAL)
    # include compiler info
    INCLUDE("${COMPILER_ID_FILE}")

  else(${COMPILER_TYPE} STREQUAL "GCC")

    message("can only set compiler environment for GCC or ICC, not for ${COMPILER}!")

  endif(${COMPILER_TYPE} STREQUAL "GCC")
endmacro()


