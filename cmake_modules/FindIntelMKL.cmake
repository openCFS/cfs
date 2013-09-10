IF(MINGW OR MSVC)

  #-----------------------------------------------------------------------------
  # If not specified by the user, try to determine proper MKL root directory.
  #-----------------------------------------------------------------------------
  IF(NOT MKL_ROOT_DIR)
    IF(NOT CFS_BUILD_OS STREQUAL CFS_TARGET_OS)
      SET(MKL_POSSIBLE_ROOT_DIRS
	"/opt/pckg/mkl_win/composer_xe_2013"
	"/opt/pckg/mkl_win/10.0.5.025"
	)
    ELSE()
      SET(MKL_POSSIBLE_ROOT_DIRS
	"e:/dev/intel/MKL/composer_xe_2013"
	"e:/dev/intel/MKL/10.0.5.025"
	)
    ENDIF()

    find_file(MKL_H
      "include/mkl.h"
      PATHS ${MKL_POSSIBLE_ROOT_DIRS}
      DOC "MKL root dir"
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH
      )

    STRING(REPLACE "/include/mkl.h" "" MKL_ROOT_DIR "${MKL_H}")
    MESSAGE(FATAL_ERROR "MKL_ROOT_DIR ${MKL_ROOT_DIR}")
  ENDIF()


  #-----------------------------------------------------------------------------
  # Read mkl.h and try to determine MKL version.
  #-----------------------------------------------------------------------------
  FILE(STRINGS "${MKL_ROOT_DIR}/include/mkl.h" MKL_HEADER)

  foreach(line IN LISTS MKL_HEADER)
    IF(line MATCHES "#define __INTEL_MKL")
      STRING(REPLACE "_" ";" MKL_LIST "${line}")
      LIST(GET MKL_LIST 4 ITEM_FOUR)

      IF(ITEM_FOUR STREQUAL "")
	LIST(GET MKL_LIST 5 MKL_MAJOR_VERSION)
        STRING(STRIP "${MKL_MAJOR_VERSION}" MKL_MAJOR_VERSION)
      ELSEIF(ITEM_FOUR STREQUAL "MINOR")
	LIST(GET MKL_LIST 6 MKL_MINOR_VERSION)
        STRING(STRIP "${MKL_MINOR_VERSION}" MKL_MINOR_VERSION)
      ELSEIF(ITEM_FOUR STREQUAL "UPDATE")
	LIST(GET MKL_LIST 6 MKL_UPDATE)
        STRING(STRIP "${MKL_UPDATE}" MKL_UPDATE)
      ENDIF()
    ENDIF()
  endforeach()

  SET(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include")

  # MESSAGE("MKL_VERSION: ${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE}")

  #-----------------------------------------------------------------------------
  # If MKL version is 11 bundle together linker libraries.
  #-----------------------------------------------------------------------------
  IF("${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE}" MATCHES "^11")

    IF(CFS_ARCH STREQUAL "I386")
      SET(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib/ia32")
      
      SET(MKL_BLAS_LIB
        # This works with MinGW GCC 4.7.1 on Windows XP
        -Wl,--start-group
	${MKL_LIB_DIR}/mkl_intel_c_dll.lib
	# ${MKL_LIB_DIR}/mkl_sequential_dll.lib
	${MKL_LIB_DIR}/mkl_intel_thread_dll.lib
	${MKL_LIB_DIR}/mkl_core_dll.lib
        -Wl,--end-group
        # libiomp5md is not needed for mkl_sequential lib
	${MKL_ROOT_DIR}/compiler/lib/ia32/libiomp5md.lib
	${CFS_SOURCE_DIR}/cfsdeps/mkl/msvc90/ia32/runtmchk.lib
	)
      SET(MKL_LAPACK_LIB ${MKL_BLAS_LIB})
      
      SET(MKL_PARDISO_LIB "")
    ELSE(CFS_ARCH STREQUAL "I386")
      
      SET(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib/intel64")
      
      SET(MKL_BLAS_LIB
        # This works with MinGW GCC 4.5.3 on CentOS/Oracle 6
        -Wl,--start-group
	${MKL_LIB_DIR}/mkl_intel_lp64.lib
	# ${MKL_LIB_DIR}/mkl_intel_thread.lib
	${MKL_LIB_DIR}/mkl_sequential.lib
	${MKL_LIB_DIR}/mkl_core.lib
        -Wl,--end-group
	# ${MKL_ROOT_DIR}/compiler/lib/intel64/libiomp5md.lib
	${CFS_SOURCE_DIR}/cfsdeps/mkl/msvc100/amd64/runtmchk.lib
	)
      SET(MKL_LAPACK_LIB ${MKL_BLAS_LIB})
      
      SET(MKL_PARDISO_LIB "")
    ENDIF(CFS_ARCH STREQUAL "I386")
  ELSE()  

    IF(CFS_ARCH STREQUAL "I386")
      SET(MKL_LIB_DIR "${MKL_ROOT_DIR}/ia32/lib")

      SET(MKL_BLAS_LIB
        -Wl,--start-group
        ${MKL_LIB_DIR}/mkl_intel_c.lib
        ${MKL_LIB_DIR}/mkl_intel_thread.lib
        ${MKL_LIB_DIR}/mkl_core.lib
        -Wl,--end-group
        ${MKL_LIB_DIR}/libiomp5md.lib
        )
      IF(MINGW)
        LIST(APPEND MKL_BLAS_LIB 
          ${CFS_SOURCE_DIR}/cfsdeps/mkl/msvc90/ia32/runtmchk.lib
        )
      ENDIF()

      SET(MKL_LAPACK_LIB ${MKL_BLAS_LIB})
      
      SET(MKL_PARDISO_LIB
        ${MKL_LIB_DIR}/mkl_solver.lib
        )
    ELSE(CFS_ARCH STREQUAL "I386")
      SET(MKL_LIB_DIR "${MKL_ROOT_DIR}/em64t/lib")

      SET(MKL_BLAS_LIB
        # This works with MinGW GCC 4.5.3 on CentOS/Oracle 6
        -Wl,--start-group
        ${MKL_LIB_DIR}/mkl_intel_lp64.lib
        ${MKL_LIB_DIR}/mkl_intel_thread.lib
        ${MKL_LIB_DIR}/mkl_core.lib
        -Wl,--end-group
        ${MKL_LIB_DIR}/libiomp5mt.lib
        #    wrap-chkstk
	# ${MKL_LIB_DIR}/libguide.lib
        )
      IF(MINGW)
        LIST(APPEND MKL_BLAS_LIB 
          ${CFS_SOURCE_DIR}/cfsdeps/mkl/msvc90/amd64/runtmchk.lib
        )
      ENDIF()
      SET(MKL_LAPACK_LIB ${MKL_BLAS_LIB})
      
      SET(MKL_PARDISO_LIB
        ${MKL_LIB_DIR}/mkl_solver_lp64.lib
        )
      
      SET(DEPS_SEQUENTIAL
	${MKL_LIB_DIR}/mkl_solver_lp64_sequential.lib
	${MKL_LIB_DIR}/mkl_intel_lp64.lib
	${MKL_LIB_DIR}/mkl_intel_thread.lib
	${MKL_LIB_DIR}/mkl_core.lib
	${MKL_LIB_DIR}/mkl_intel_lp64.lib
	${MKL_LIB_DIR}/mkl_sequential.lib
	${MKL_LIB_DIR}/libiomp5mt.lib
	${CFS_SOURCE_DIR}/cfsdeps/mkl/msvc90/amd64/runtmchk.lib
	#    ${MKL_LIB_DIR}/libguide40.lib
	#    /home/strieben/Documents/MKL/test/runtmchk.lib
	)

    ENDIF(CFS_ARCH STREQUAL "I386")
    
  ENDIF()

  SET(MKL_ROOT_DIR ${MKL_ROOT_DIR} CACHE PATH "Directory of MKL." FORCE)

ELSEIF(UNIX)

#-------------------------------------------------------------------------------
# The idea behind the algorithm implemented for finding MKL, is to let the
# make files in the MKL example directories tell us which linker flags we need.
# We therefore tell the MKLROOT (which contains a include/mkl.h) to an external
# shell script (compile_mkl_test.sh) which will then copy one of the example 
# make files to a temporary directory, make some necessary modifications to it
# and compile a simple test program which will afterwards print out the version
# of MKL. All outputs of the shell script and test executable are in CMake 
# syntax and get written to CMakeFiles/mkl.cmake which will in turn be included
# from the current script.
#
# This somewhat complicated procedure has been chosen, since Intel seems to 
# change the file and directory layout of their compiler and MKL suite even 
# with every minor version. Therefore Intel should know best, what linker 
# switches to use for what version and we happily copy that information.
#
# ATTENTION: The path to MKL should not contain any spaces!
#-------------------------------------------------------------------------------

SET (MKL_POSSIBLE_PATHS
  # Path the user may have specified via -DMKL_ROOT_DIR:PATH=...
  ${MKL_ROOT_DIR}
  # Path the user may have specified as environment variable
  $ENV{MKL_ROOT_DIR}
  # Path which may have been given in a platform_defaults_*.cmake file 
  ${MKL_ROOT_DIR_DEFAULT}
  # Path set by ifortvars.sh resp. iccvars.sh
  $ENV{MKLROOT}
  # Local paths
  /opt/intel/composerxe-2011.4.191
  /opt/intel/Compiler/11.1/069/mkl
  /opt/intel/Compiler/11.0/059/Frameworks/mkl
  /opt/intel/Compiler/11.0/081/mkl
  /opt/intel/Compiler/11.0/074/mkl
  /opt/intel/mkl/10.0.5.025
  /opt/intel/mkl/9.1.023
  /opt/intel/mkl/9.1.021
  # Paths on Woody
  /apps/intel/ComposerXE/mkl
  /apps/intel/ComposerXE/composerxe-2011.4.191/mkl
  /apps/intel/Compiler/11.0/069/mkl
  /apps/intel/mkl/10.0.011
  /apps/intel/ict/3.0/cmkl/9.0
 )

#-------------------------------------------------------------------------------
# Find the first directory containing include/mkl.h
#-------------------------------------------------------------------------------
FIND_FILE(MKL_ROOT_DIR
  NAMES include/mkl.h
  PATHS ${MKL_POSSIBLE_PATHS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  )

#-------------------------------------------------------------------------------
# Replace include/mkl.h with nothing to get the MKL root dir.
#-------------------------------------------------------------------------------
STRING(REPLACE "/include/mkl.h" "" MKL_ROOT_DIR "${MKL_ROOT_DIR}")

#-------------------------------------------------------------------------------
# Make sure that there are no spaces in MKLROOT.
#-------------------------------------------------------------------------------
STRING(REGEX MATCH " " SPACES "${MKL_ROOT_DIR}")
STRING(LENGTH "${SPACES}" LEN)
IF(LEN GREATER 0)
  MESSAGE(SEND_ERROR "No spaces are allowed in MKL root directory: ${MKL_ROOT_DIR}")
ENDIF(LEN GREATER 0)

IF(NOT MKL_ROOT_DIR)
  MESSAGE(FATAL_ERROR "Could not find Intel MKL library. This script just "
    "checks explicitly for MKL installations in certain locactions. "
    "If you have installed MKL into a different directory pleas do one of the following things: "
    " 1. Set an environment variable MKL_ROOT_DIR or MKLROOT " 
    " 2. Set '-DMKL_ROOT_DIR:PATH=...' on the CMake command line "
    " 3. Or adapt the MKL_POSSIBLE_PATHS at the beginning of cmake_modules/FindIntelMKL.cmake. "
    "You may also just copy MKL from LSE /home/shareAll/linux_bin/intel/mkl. "
    "Another option is to change CFS_BLAS_LAPACK and CFS_PARDISO to different "
    "implementations. GotoBLAS and ACML should provide similar performance like MKL.")
ELSE(NOT MKL_ROOT_DIR)
  SET(MKL_FOUND 1)
ENDIF(NOT MKL_ROOT_DIR)

SET(MKL_ROOT_DIR ${MKL_ROOT_DIR} CACHE PATH "Directory of MKL." FORCE)
SET(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include")

#-------------------------------------------------------------------------------
# Determine version of MKL and linker flags by compiling a little test program
# and   parsing  the   command   line   generated  by   a   makefile  from 
# MKLROOT/examples/solver
#-------------------------------------------------------------------------------
If(WIN32)
  STRING(REGEX REPLACE ".*MKL/"
    "" MKL_VERSION
    ${MKL_ROOT_DIR})
ELSE(WIN32)
  #-----------------------------------------------------------------------------
  # Specify desired architecture and output type for make file.
  #-----------------------------------------------------------------------------
  IF(CFS_ARCH STREQUAL "I386")
    SET(MKL_ARCH_ID "lib32")
  ENDIF(CFS_ARCH STREQUAL "I386")

  IF(CFS_ARCH STREQUAL "X86_64")
    SET(MKL_ARCH_ID "libem64t")
  ENDIF(CFS_ARCH STREQUAL "X86_64")

  IF(CFS_ARCH STREQUAL "IA64")
    SET(MKL_ARCH_ID "lib64")
  ENDIF(CFS_ARCH STREQUAL "IA64")

  #-----------------------------------------------------------------------------
  # Specify compiler type for make file.
  #-----------------------------------------------------------------------------
  IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
    SET(MKL_COMPILER_ID "intel")
  ELSE(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
    SET(MKL_COMPILER_ID "gnu")
  ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC") 

  #-----------------------------------------------------------------------------
  # We always want to link against the parallel version of MKL.
  #-----------------------------------------------------------------------------
#  IF(USE_OPENMP)
    SET(MKL_THREADING_ID "parallel")
#  ELSE(USE_OPENMP)
#    SET(MKL_THREADING_ID "sequential")
#  ENDIF(USE_OPENMP)  

  #-----------------------------------------------------------------------------
  # Configure our little build script using the previously defined variables.
  #-----------------------------------------------------------------------------
  FILE(MAKE_DIRECTORY "${CFS_BINARY_DIR}/tmp/mkl_test")
  CONFIGURE_FILE(
    "${CFS_SOURCE_DIR}/share/scripts/compile_mkl_test.sh.in"
    "${CFS_BINARY_DIR}/tmp/mkl_test/compile_mkl_test.sh" @ONLY)

  #-----------------------------------------------------------------------------
  # Execute the build script and put the generated CMake code into mkl.cmake.
  #-----------------------------------------------------------------------------
  EXECUTE_PROCESS(
    COMMAND sh "${CFS_BINARY_DIR}/tmp/mkl_test/compile_mkl_test.sh"
    OUTPUT_VARIABLE MKL_INFO
    RESULT_VARIABLE RETVAL)

  IF(RETVAL EQUAL 0)
    FILE(WRITE "${CFS_BINARY_DIR}/CMakeFiles/mkl.cmake" ${MKL_INFO})
    #---------------------------------------------------------------------------
    # Finally just include the linker flags and version info.
    #---------------------------------------------------------------------------
    INCLUDE(${CFS_BINARY_DIR}/CMakeFiles/mkl.cmake)
  ELSE(RETVAL EQUAL 0)
    MESSAGE(SEND_ERROR "A problem occurred during determination of MKL linker
                        flags. Please run ${CFS_BINARY_DIR}/tmp/compile_mkl_test.sh
                        by hand to investigate the problem.")
  ENDIF(RETVAL EQUAL 0)

ENDIF(WIN32)
# MESSAGE("MKL_ROOT_DIR ${MKL_ROOT_DIR}")
# MESSAGE("MKL_MAJOR_VERSION ${MKL_MAJOR_VERSION}")
# MESSAGE("MKL_MINOR_VERSION ${MKL_MINOR_VERSION}")

SET(CFS_MKL_VERSION ${MKL_MAJOR_VERSION})

#-------------------------------------------------------------------------------
# Mark libraries as advanced variables.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(MKL_BLAS_LIB)
MARK_AS_ADVANCED(MKL_LAPACK_LIB)
MARK_AS_ADVANCED(MKL_PARDISO_LIB)
MARK_AS_ADVANCED(MKL_ROOT_DIR)

# TODO: libguide und libiomp durch die von MKL ersetzen.
# TODO: für libmkl_intel_lp64.a libmkl_gf_lp64.a für gnu einsetzen.
ENDIF() # MINGW OR MSVC

#-------------------------------------------------------------------------------
# Set BLAS, LAPACK and PARDISO libraries depending on the MKL version.
#-------------------------------------------------------------------------------
IF(CFS_BLAS_LAPACK STREQUAL "MKL")
  SET(BLAS_LIBRARY "${MKL_BLAS_LIB}")
  IF(MKL_MAJOR_VERSION LESS 10)
    SET(LAPACK_LIBRARY "${MKL_LAPACK_LIB}")
  ENDIF(MKL_MAJOR_VERSION LESS 10)  
ENDIF(CFS_BLAS_LAPACK STREQUAL "MKL")
SET(PARDISO_LIBRARY "${MKL_PARDISO_LIB}")

