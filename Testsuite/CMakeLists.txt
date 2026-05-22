#-------------------------------------------------------------------------------
#              Main CMakeLists.txt file for CFS++ Testsuite.
# 
# To use testsuite please specify the base path to this file in cmake for CFS++.
#-------------------------------------------------------------------------------
# 3.16 is minimum for advaned options of FindPython.cmake
cmake_minimum_required(VERSION 3.16)

# check if we run from within cfs (i.e. testsuite is included)
if(NOT PROJECT_NAME) # we are in standalone mode
  message(STATUS "Configuring Testsuite in standalone mode")
  set(CMAKE_MAKE_PROGRAM "none") # we do not need any
  project(TESTSUITE NONE)
  set(CFS_INSTALL_DIR "${CMAKE_BINARY_DIR}" CACHE PATH "location of cfs installation to test")
  # find the binary
  if(UNIX) # linux and mac
    set(CFS_EXECUTABLE "${CFS_INSTALL_DIR}/bin/cfs")
  else() # windows
    set(CFS_EXECUTABLE "${CFS_INSTALL_DIR}/bin/cfs.exe")
  endif()
  # message("CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}")
  if(EXISTS "${CFS_EXECUTABLE}")
    # run cfs and prinf version info
    execute_process(COMMAND ${CFS_INSTALL_DIR}/bin/cfs --version --noColor
                    RESULT_VARIABLE res
                    OUTPUT_VARIABLE cfs_version_output)
    message(STATUS "running: cfs --version\n${cfs_version_output}")
    # run cfs to get build config information
    execute_process(COMMAND "${CFS_EXECUTABLE}" "--dependencies" "${CMAKE_CURRENT_BINARY_DIR}/config.cmake"
	  RESULT_VARIABLE res 
	  OUTPUT_VARIABLE STDOUT)
    if(NOT ${res})
      set(CFS_BUILD_CONFIG_FILE "${CMAKE_BINARY_DIR}/config.cmake")
      set(TESTSUITE_DIR ${CMAKE_SOURCE_DIR})
    else()
      message(ERROR "res=${res}, fail: ${out}")
    endif()
  else() # could not find executable
    message(STATUS "cannot find cfs binary at ${CFS_EXECUTABLE}: is CFS_INSTALL_DIR correct?")
    if(TESTSUITE_CFS) # we want to test CFS but cannot find it
	    message(FATAL_ERROR "  with TESTSUITE_CFS=ON we want to test CFS but did not find it!")
    else() # we do not need to test CFS
      message(STATUS "  since TESTSUITE_CFS=OFF we do not need a CFS binary, what else is in ${CFS_INSTALL_DIR} ?")
      if(EXISTS "${CFS_INSTALL_DIR}/share")
        message(STATUS "    testing the scripts in ${CFS_INSTALL_DIR}/share")
      else()
        message(STATUS "    testing python scripts from the source directory ${CMAKE_SOURCE_DIR}/../share")
        set(CFS_INSTALL_DIR "${CMAKE_SOURCE_DIR}/..")
      endif()
    endif()
    set(TESTSUITE_DIR ${CMAKE_SOURCE_DIR})
  endif()
  # include the common config file
  include("ctest_settings.cmake")
  # ctest must be included from the top level CMakeLists
  include(CTest)
  enable_testing()
  mark_as_advanced(BUILD_TESTING) # hide this CTests internal option
else() # included from within CFS
  set(CFS_INSTALL_DIR "${CMAKE_BINARY_DIR}")
endif(NOT PROJECT_NAME)

# replace possible backslashes in path names by forward slashes (when passed from windows)
if(TESTSUITE_DIR MATCHES "\\\\")
  message(STATUS "Replacing \\ in TESTSUITE_DIR=${TESTSUITE_DIR} please use forward slashes in path names")
  STRING(REGEX REPLACE "\\\\" "/" TESTSUITE_DIR ${TESTSUITE_DIR})
endif()
if(CFS_BUILD_CONFIG_FILE MATCHES "\\\\")
  message(STATUS "Replacing \\ in CFS_BUILD_CONFIG_FILE=${CFS_BUILD_CONFIG_FILE} please use forward slashes in path names")
  STRING(REGEX REPLACE "\\\\" "/" CFS_BUILD_CONFIG_FILE ${CFS_BUILD_CONFIG_FILE})
endif()

# set the directory of the CFS installation (legacy name)
set(CFS_BINARY_DIR ${CFS_INSTALL_DIR})

# set the binary dir of the testsuite
# can we use CMAKE_BINARY_DIR for both options (testsuite incuded from cfs or standalone)? 
set(TESTSUITE_BIN_DIR "${CMAKE_BINARY_DIR}/testsuite")

# CFS_BUILD_CONFIG_FILE can be set
#  a) above if Testsuite is used standalone
#  b) via -D if one wants to override defaults for testing.
if(CFS_BUILD_CONFIG_FILE)
  message(STATUS "Testsuite is loading CFS build config from file '${CFS_BUILD_CONFIG_FILE}' to define tests")
  include("${CFS_BUILD_CONFIG_FILE}")
endif()

# configure the testsuite

# set the number of thread to use for the testing runs. We set to > 1 to test parallelization. 
# The golden rule for fast tests ist minimal number of threads and testing parallel with 
# ctest -j N where N*TESTSUITE_THREADS=system threads (including hyperthreading)
# fastest would be TESTSUITE_THREADS=1, e.g. 2 min with ctest -j 8. TESTSUITE_THREADS=2 with ctest -j 4 is 3 minutes, then it becomes even slower.
set(TESTSUITE_THREADS "2" CACHE STRING "Number of threads to use for the cfs testing runs, i.e. use argument '-t TESTSUITE_THREADS' in cfs runs.")
# check if integer is set
if(NOT TESTSUITE_THREADS MATCHES "^[1-9][0-9]*$")
  message(FATAL_ERROR "TESTSUITE_THREADS must be integer valued.")
endif()

# set the mode of cfstool
set(TESTSUITE_CFSTOOL_MODE "relL2diff" CACHE STRING "The mode used to compare hdf5 files by cfstool.")
set_property(CACHE TESTSUITE_CFSTOOL_MODE PROPERTY STRINGS scalardiff relL2diff absL2diff)

# set the test timeout
set(TESTSUITE_TIMEOUT 300 CACHE STRING "timeout in s for each test in the testsuite (sets DART_TESTING_TIMEOUT)")
set(DART_TESTING_TIMEOUT "${TESTSUITE_TIMEOUT}" CACHE INTERNAL "=TESTSUITE_TIMEOUT") # see https://cmake.org/cmake/help/v3.3/manual/ctest.1.html#ctest-test-step
mark_as_advanced(TESTSUITE_TIMEOUT)

# check options
# check if we have dependent options set correctly
# In order to perform tests we need CFSTOOL!
if(CFS_TESTSUITE AND NOT BUILD_CFSTOOL)
    message(FATAL_ERROR "CFS_TESTSUITE=ON requires BUILD_CFSTOOL=ON.")
endif()

# enbale some custom ctest settings
configure_file("CTestCustom.cmake.in" "${CMAKE_BINARY_DIR}/CTestCustom.cmake" @ONLY)

#-------------------------------------------------------------------------------
# Include some macros. Especially the one for generating test names out of
# directory names.
#-------------------------------------------------------------------------------
include("TestMacros.cmake")

# Set file name for configured version of CFS++ standard test macro (there is a .in version)
# * cfstool for comparisson of h5 files (default)
# * compare_info_xml.py for comparison of special info.xml stuff (eigenfrequencies) via parameter "info.xml"
set(CFS_STANDARD_TEST "${CFS_BINARY_DIR}/testsuite/StandardTest.cmake")
set(PYTHON_TEST "${CFS_BINARY_DIR}/testsuite/PythonTest.cmake")

# Define path to info.xml compare script (for eigenvalue and other stuff)
set(COMPARE_INFO_XML "${CMAKE_CURRENT_SOURCE_DIR}/compare_info_xml.py")

# Define path to the compare_hist_file.py to compare text outputs of cfs
set(COMPARE_HIST_FILE "${CMAKE_CURRENT_SOURCE_DIR}/compare_hist_file.py")

# python tests including matviz
set(MATVIZ_PY "${CFS_BINARY_DIR}/share/python/matviz.py")
set(BASECELL_PY "${CFS_BINARY_DIR}/share/python/basecell.py")

# Checks matviz.py results
set(COMPARE_IMAGES "${CMAKE_CURRENT_SOURCE_DIR}/compare_images.py")

include(CMakeDependentOption)
# Enables tests of Python scripts e.g. via doctest
cmake_dependent_option(TESTSUITE_PYTHON "testing of python modules" ON "CFS_TESTSUITE OR ${PROJECT_NAME} MATCHES TESTSUITE" OFF)
mark_as_advanced(TESTSUITE_PYTHON)

# option to use tests which need 'compare_info_xml.py'
#set(TESTSUITE_INFO_XML ON CACHE BOOL "Include tests which need 'compare_info_xml.py'")
cmake_dependent_option(TESTSUITE_INFO_XML "Include tests which need 'compare_info_xml.py'" ON "CFS_TESTSUITE OR ${PROJECT_NAME} MATCHES TESTSUITE" OFF)
mark_as_advanced(TESTSUITE_INFO_XML)

# Enables tests of Python scripts e.g. via doctest
cmake_dependent_option(TESTSUITE_CFS "testing of CFS (the usual tests in TESTSUITE)" ON "CFS_TESTSUITE OR ${PROJECT_NAME} MATCHES TESTSUITE" OFF)
mark_as_advanced(TESTSUITE_CFS)

# check for python
if((TESTSUITE_PYTHON OR TESTSUITE_INFO_XML) AND NOT "$ENV{CI_JOB_STAGE}" MATCHES "build" )
  # in case of USE_EMBEDDED_PYTHON, cfs also searched for Python -> see FindCFSDEPS.cmake
  # use -DPython_ROOT_DIR=... to help in case or check cmake's FindPython for other hints
  set(Python_FIND_VIRTUALENV "FIRST")
  find_package(Python MODULE COMPONENTS Interpreter)
  if(NOT Python_Interpreter_FOUND OR NOT Python_EXECUTABLE)
    message(WARNING "cannot find Python for testsuite - set Python_ROOT_DIR or Python_EXECUTABLE")
  else()
    message(STATUS "using ${Python_EXECUTABLE} for python related tests (see options TESTSUITE_PYTHON and TESTSUITE_INFO_XML)")
  endif()

  # check for required python modules
  set(PythonModules_FIND_COMPONENTS "")
  # use the OPTIONAL_COMPONENTS to define our custom error message
  set(OPTIONAL_COMPONENTS lxml argparse os numpy h5py vtk scipy matplotlib)
  set(PythonModules_FIND_OPTIONAL_COMPONENTS ${OPTIONAL_COMPONENTS})
  # this will only work on the first execution with a set Python_EXECUTABLE
  include("${CMAKE_CURRENT_SOURCE_DIR}/FindPythonModules.cmake") # from https://github.com/cmake-basis/find-modules
  #cmake_print_variables(PythonModules_FOUND PythonModules_OPTIONAL_FOUND)
  if(NOT PythonModules_OPTIONAL_FOUND) # report which ones are missong
    #message("${OPTIONAL_COMPONENTS}")
    foreach(module IN ITEMS ${OPTIONAL_COMPONENTS})
      #message(STATUS "${module} ${PythonModules_${module}_FOUND}")
      if(NOT ${PythonModules_${module}_FOUND})
	      message(WARNING "Python module '${module}' missing!")
      endif()
    endforeach()
  endif()
  mark_as_advanced(PythonModules_DIR)

  # Find h5diff from hdf5-tools — used by Python-script tests to compare
  # .h5 result against .h5ref reference.  Independent of cfstool so the CI
  # pipeline (which does not build cfstool) can still run these tests.
  find_program(H5DIFF_EXECUTABLE NAMES h5diff)
  if(TESTSUITE_PYTHON AND NOT H5DIFF_EXECUTABLE)
    message(WARNING "h5diff not found — Python HDF5 comparison tests will fail. "
                    "Install hdf5-tools (apt: hdf5-tools, conda: hdf5).")
  else()
    message(STATUS "using ${H5DIFF_EXECUTABLE} for Python HDF5 comparison tests")
  endif()
endif()

#-------------------------------------------------------------------------------
# Set default tolerance for scalar difference tests using cfstool.
#-------------------------------------------------------------------------------
set(CFS_DEFAULT_EPSILON "1e-6")

# add packing and build test
#  - only used in CI pipeline
if("$ENV{CI_JOB_STAGE}" MATCHES "build")
  message(STATUS "running in CI pipeline 'build' stage detected -> creating build and pack tests")
  # Create distributable binary via cpack running ctest, alternatively use create_binary_archive.py
  if(UNIX) # pack to a tar.gz archive
    add_test(NAME Pipeline_Build_cpack-TGZ COMMAND cpack -G TGZ --verbose WORKING_DIRECTORY "${CFS_BINARY_DIR}")
  else() # probably windows, we do zip
    add_test(NAME Pipeline_Build_cpack-ZIP COMMAND cpack -G ZIP --verbose WORKING_DIRECTORY "${CFS_BINARY_DIR}")
  endif()
  # Use an additional test during the build stage to grep the compiler errors
  add_test(NAME Pipeline_Build_check-compiler-errors COMMAND grep "0 Compiler errors" "${CFS_BINARY_DIR}/build.log")
  message(STATUS "submitting tests to: ${SUBMIT_URL}")
endif()

# Subdirectory for normal test suite
if(TESTSUITE_CFS)
  message(STATUS "Testsuite will be testing: ${CFS_EXECUTABLE}")
  add_subdirectory("TESTSUIT")
endif()

# Subdirectory for pypthon doctest, matviz visualization test suite and basecell
if(TESTSUITE_PYTHON)
  add_subdirectory("PYTHON")
endif()

# Configure CFS++ standard test macro.
configure_file("StandardTest.cmake.in" "${CFS_STANDARD_TEST}" @ONLY)

# python tests including matviz.py is used in the optimization group
configure_file("PythonTest.cmake.in" "${PYTHON_TEST}" @ONLY)
