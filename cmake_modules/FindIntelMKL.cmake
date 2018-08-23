#-------------------------------------------------------------------------------
# This script  is responsible for the  determination of the correct  include and
# linker  parameters for  the Intel  Math Kernel  Library (MKL).  We distinguish
# between Windows/MinGW builds and Unix (Linux  and MacOS X) builds. On Linux we
# just support the  MinGW (cross-)compilers at the moment,  which are officially
# incompatible  with  MKL.   This  is  due  to  the   different  OpenMP  runtime
# environments for  GCC and  Intel/MS compilers.  Anyway, we  found a  number of
# working combinations  of GCC and  MKL versions. These combinations  are pretty
# much hardcoded at the moment.
#
# On Unix  we determine the  required linker flags by  compiling one of  the MKL
# examples and reading the flags from the Makefile output (cf. below).
# another option: https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#-------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Function to read MKL version from header file mkl.h or mkl_version.h
#-----------------------------------------------------------------------------
function(MKL_VERSION_FROM_HEADER)

  IF(EXISTS "${MKL_ROOT_DIR}/include/mkl_version.h")
    FILE(STRINGS "${MKL_ROOT_DIR}/include/mkl_version.h" MKL_HEADER)
  ELSEIF(EXISTS "${MKL_ROOT_DIR}/include/mkl.h")
    FILE(STRINGS "${MKL_ROOT_DIR}/include/mkl.h" MKL_HEADER)
  ELSE(EXISTS "${MKL_ROOT_DIR}/include/mkl_version.h")
    SET(MSG "Please download the file ")
    SET(MSG "${MSG}'${CFS_DS_WEBDAV}/cfsdeps/sources/mkl/mkl_win.zip'")
    SET(MSG "${MSG}, unpack it and set a proper MKL_ROOT_DIR.")

    colormsg(HIRED "${MSG}")
    MESSAGE(FATAL_ERROR "MKL for Windows could not be found!")
  ENDIF(EXISTS "${MKL_ROOT_DIR}/include/mkl_version.h")

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

  SET(MKL_MAJOR_VERSION ${MKL_MAJOR_VERSION} PARENT_SCOPE)
  SET(MKL_MINOR_VERSION ${MKL_MINOR_VERSION} PARENT_SCOPE)
  SET(MKL_UPDATE ${MKL_UPDATE} PARENT_SCOPE)
  # MESSAGE("MKL_VERSION: ${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE}")

endfunction(MKL_VERSION_FROM_HEADER)


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


  #---------------------------------------------------------------------------
  # Read mkl version from header files
  #---------------------------------------------------------------------------
  MKL_VERSION_FROM_HEADER()
  SET(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include")
  
  SET(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib/intel64_win")

  # mkl_solver* libs are deprecated as of v10
  # see: https://web.archive.org/web/20091027212420/https://software.intel.com/en-us/articles/mkl_solver_libraries_are_deprecated_libraries_since_version_10_2_Update_2/

  SET(MKL_BLAS_LIB
    # This works with MinGW GCC 4.5.3 on CentOS/Oracle 6
    -Wl,--start-group
    ${MKL_LIB_DIR}/mkl_intel_lp64.lib
    ${MKL_LIB_DIR}/mkl_intel_thread.lib
    ${MKL_LIB_DIR}/mkl_core.lib
    ${MKL_LIB_DIR}/mkl_blas95_lp64.lib
    -Wl,--end-group
    #    wrap-chkstk
# ${MKL_LIB_DIR}/libguide.lib
    )

  SET(MKL_LAPACK_LIB ${MKL_BLAS_LIB})

  SET(DEPS_SEQUENTIAL
  ${MKL_LIB_DIR}/mkl_intel_lp64.lib
  ${MKL_LIB_DIR}/mkl_intel_thread.lib
  ${MKL_LIB_DIR}/mkl_core.lib
  ${MKL_LIB_DIR}/mkl_intel_lp64.lib
  ${MKL_LIB_DIR}/mkl_sequential.lib
  ${MKL_ROOT_DIR}/../compiler/lib/intel64_win/libiomp5md.lib # md or mt? I have no idea
  -openmp
  #${MKL_ROOT_DIR}/../msvcrt/msvc90/amd64/runtmchk.lib
  #    ${MKL_LIB_DIR}/libguide40.lib
  )

  SET(MKL_ROOT_DIR ${MKL_ROOT_DIR} CACHE PATH "Directory of MKL." FORCE)

ELSEIF(UNIX AND NOT CMAKE_CROSSCOMPILING) # end MINGW OR MSVC 
  #-------------------------------------------------------------------------------
  # The idea behind the algorithm implemented for	 finding MKL, is to let the make
  # files in the MKL  example directories tell us which linker  flags we need.  We
  # therefore tell  the MKLROOT  (which contains a  include/mkl.h) to  an external
  # shell script	(compile_mkl_test.sh) which  will then copy  one of  the example
  # make files to	 a temporary directory, make some necessary  modifications to it
  # and compile a simple test program  which will afterwards print out the version
  # of MKL.  All outputs	of the	shell script  and test	executable are	in CMake
  # syntax and get written to CMakeFiles/mkl.cmake  which will in turn be included
  # from the current script.
  # 
  # This	somewhat complicated  procedure has  been chosen,  since Intel	seems to
  # change the file and directory layout of their compiler and MKL suite even with
  # every minor version. Therefore Intel should know best, what linker switches to
  # use for what version and we happily copy that information.
  # 
  # ATTENTION: The path to MKL should not contain any spaces!
  #-------------------------------------------------------------------------------
  
  SET (MKL_POSSIBLE_PATHS
    # Path the user may have specified via -DMKL_ROOT_DIR:PATH=...
    ${MKL_ROOT_DIR}
    # Path the user may have specified as environment variable
    $ENV{MKL_ROOT_DIR}
    $ENV{MKL_BASE}
    # Path which may have been given in a platform_defaults_*.cmake file 
    ${MKL_ROOT_DIR_DEFAULT}
    # Path set by ifortvars.sh resp. iccvars.sh
    # This happens automatically when loading mkl 
    # modules at Erlangens RRZE HPC cluster
    $ENV{MKLROOT}
    # Local paths
    /opt/intel/composerxe-2011.4.191
    /opt/intel/Compiler/11.1/069/mkl
    /opt/intel/Compiler/11.0/059/Frameworks/mkl
    /opt/intel/Compiler/11.0/081/mkl
    /opt/intel/Compiler/11.0/074/mkl
    /opt/intel/mkl/10.0.5.025
    /opt/intel/composer_xe_2015.2.164/mkl
    # path for mdmt
    /share/programs/intel-mkl/latest/mkl
    # Paths on Lima
    /apps/intel/ComposerXE2013/composer_xe_2013.5.192/mkl
    # intel 2016 tools
    /opt/intel/compilers_and_libraries/linux/mkl
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
      "implementations. OpenBLAS and should provide similar performance to MKL.")
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
  
    #---------------------------------------------------------------------------
    # Read mkl version from header files
    #---------------------------------------------------------------------------
    MKL_VERSION_FROM_HEADER()
  
    #-----------------------------------------------------------------------------
    # Specify desired architecture and output type for make file.
    #-----------------------------------------------------------------------------
  
    IF(CFS_ARCH STREQUAL "X86_64")
      SET(MKL_ARCH_ID "intel64")
    ENDIF()
  
    #-----------------------------------------------------------------------------
    # Specify compiler type for make file.
    #-----------------------------------------------------------------------------
    IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
      SET(MKL_COMPILER_ID "intel")
    ELSE()
      SET(MKL_COMPILER_ID "gnu")
    ENDIF() 
  
    #-----------------------------------------------------------------------------
    # We always want to link against the parallel version of MKL.
    #-----------------------------------------------------------------------------
    IF(${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE} VERSION_LESS 11.3.3)
      SET(MKL_THREADING_ID "parallel")
    ELSE()
      SET(MKL_THREADING_ID "omp")
    ENDIF()
  
    #-----------------------------------------------------------------------------
    # Make a directory for the test compile
    #-----------------------------------------------------------------------------
    SET(COMPILE_MKL_TEST_DIR "${CFS_BINARY_DIR}/tmp/mkl_test")
    FILE(MAKE_DIRECTORY "${COMPILE_MKL_TEST_DIR}")
  
    #-----------------------------------------------------------------------------
    # Search for the proper $MKLROOT/examples/solver[c]/makefile.
    #-----------------------------------------------------------------------------
    foreach(fname "solver" "solverc") # search for example directly 
      file(GLOB_RECURSE MKL_SOLVERMAKEFILE FOLLOW_SYMLINKS "${MKL_ROOT_DIR}/*/examples/${fname}/[Mm]akefile")
        if(MKL_SOLVERMAKEFILE)
          break()
        endif(MKL_SOLVERMAKEFILE)
    endforeach(fname)
    
    if(NOT MKL_SOLVERMAKEFILE) # serach for archive and extract it if found
      foreach(fname "examples_core.tgz" "examples_core_c.tgz")
        file(GLOB_RECURSE MKL_SOLVERMAKEFILE_TGZ FOLLOW_SYMLINKS "${MKL_ROOT_DIR}/*/${fname}")
        #message("${fname}:${MKL_SOLVERMAKEFILE_TGZ}")
        if(MKL_SOLVERMAKEFILE_TGZ)
          execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${MKL_SOLVERMAKEFILE_TGZ}
            WORKING_DIRECTORY "${COMPILE_MKL_TEST_DIR}"
          )
          file(GLOB_RECURSE MKL_SOLVERMAKEFILE FOLLOW_SYMLINKS "${COMPILE_MKL_TEST_DIR}/solverc/makefile")
          break()
        endif(MKL_SOLVERMAKEFILE_TGZ)
      endforeach(fname)
    endif(NOT MKL_SOLVERMAKEFILE)
  
    if(NOT MKL_SOLVERMAKEFILE)
      message(FATAL_ERROR "Did not find example solver makefile. Have you sourced MKL in bashrc?")
    endif(NOT MKL_SOLVERMAKEFILE)
  
    #-----------------------------------------------------------------------------
    # Configure our little build script using the previously defined variables.
    #-----------------------------------------------------------------------------
    CONFIGURE_FILE(
      "${CFS_SOURCE_DIR}/share/scripts/compile_mkl_test.sh.in"
      "${COMPILE_MKL_TEST_DIR}/compile_mkl_test.sh" @ONLY)
   
    # Make sure we can execute the script
    EXECUTE_PROCESS(
      COMMAND chmod 755 "${COMPILE_MKL_TEST_DIR}/compile_mkl_test.sh"
      WORKING_DIRECTORY "${COMPILE_MKL_TEST_DIR}")
    
    #-----------------------------------------------------------------------------
    # Execute the build script and put the generated CMake code into mkl.cmake.
    #-----------------------------------------------------------------------------
    EXECUTE_PROCESS(
      COMMAND sh "${COMPILE_MKL_TEST_DIR}/compile_mkl_test.sh"
      OUTPUT_VARIABLE MKL_INFO
      WORKING_DIRECTORY "${COMPILE_MKL_TEST_DIR}"
      RESULT_VARIABLE RETVAL)
  
    IF(RETVAL EQUAL 0)
      FILE(WRITE "${CFS_BINARY_DIR}/CMakeFiles/mkl.cmake" ${MKL_INFO})
      FILE(REMOVE_RECURSE "${COMPILE_MKL_TEST_DIR}")
  
      #---------------------------------------------------------------------------
      # Finally just include the linker flags and version info.
      #---------------------------------------------------------------------------
      INCLUDE(${CFS_BINARY_DIR}/CMakeFiles/mkl.cmake)
      #---------------------------------------------------------------------------
      # Create a custom target to copy over libiomp5.so
      # here one could include the option USE_OMP=OFF which does not need to not copy libomp5
      #---------------------------------------------------------------------------
      ADD_CUSTOM_TARGET(mkl_libomp5 ALL
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${MKL_OMP_DIR}/libiomp5.so ${LIBRARY_OUTPUT_PATH}/libiomp5.so
        BYPRODUCTS ${LIBRARY_OUTPUT_PATH}/libiomp5.so
        COMMENT "Copying libiomp5.so to ${LIBRARY_OUTPUT_PATH} folder..."
        USES_TERMINAL)
  
    ELSE()
      MESSAGE("executing ${COMPILE_MKL_TEST_DIR}/compile_mkl_test.sh in ${COMPILE_MKL_TEST_DIR} failed")
      MESSAGE("RETVAL=${RETVAL} and MKL_INFO=${MKL_INFO}")
      MESSAGE(SEND_ERROR "A problem occurred during determination of MKL linker
                          flags. Please run ${CFS_BINARY_DIR}/tmp/mkl_test/compile_mkl_test.sh
                          by hand to investigate the problem. Script output was\n ${MKL_INFO}")
    ENDIF()
  
  SET(CFS_MKL_VERSION ${MKL_MAJOR_VERSION})
  
  #-------------------------------------------------------------------------------
  # Mark libraries as advanced variables.
  #-------------------------------------------------------------------------------
  MARK_AS_ADVANCED(MKL_BLAS_LIB)
  MARK_AS_ADVANCED(MKL_LAPACK_LIB)
  MARK_AS_ADVANCED(MKL_PARDISO_LIB)
  MARK_AS_ADVANCED(MKL_ROOT_DIR)
  MARK_AS_ADVANCED(MKL_OMP_DIR)
  
  # TODO: libguide und libiomp durch die von MKL ersetzen.
  # TODO: fuer libmkl_intel_lp64.a libmkl_gf_lp64.a fuer gnu einsetzen.

ENDIF() # end of UNIX AND NOT CMAKE_CROSSCOMPILING

#-------------------------------------------------------------------------------
# Set BLAS, LAPACK and PARDISO libraries depending on the MKL version.
#-------------------------------------------------------------------------------

# see also External_OpenBLAS and External_LAPACK, where the setting is quite different:
# e.g. BLAS_LIBRARY=-Wl,--start-group;/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_intel_lp64.a;/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_gnu_thread.a;/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64/libmkl_core.a;-Wl,--end-group;-L/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/../compiler/lib/intel64;-liomp5;-lpthread;-lm;-ldl
# e.g. LAPACK_LIBRARY=
IF(CFS_BLAS_LAPACK STREQUAL "MKL")
  SET(BLAS_LIBRARY "${MKL_BLAS_LIB}")
  IF(MKL_MAJOR_VERSION LESS 10)
    SET(LAPACK_LIBRARY "${MKL_LAPACK_LIB}")
  ENDIF(MKL_MAJOR_VERSION LESS 10)  
ENDIF(CFS_BLAS_LAPACK STREQUAL "MKL")
SET(PARDISO_LIBRARY "${MKL_PARDISO_LIB}")

#-------------------------------------------------------------------------------
# Status message of found MKL
#-------------------------------------------------------------------------------
IF(DEFINED MKL_UPDATE)
  MESSAGE(STATUS "Found Intel MKL version ${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.${MKL_UPDATE}.")
ELSE()
  MESSAGE(STATUS "Found Intel MKL version ${MKL_MAJOR_VERSION}.${MKL_MINOR_VERSION}.")
ENDIF()



