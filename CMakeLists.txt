#-----------------------------------------------------------------------------
# Set the version number of CMake which is required to build CFS++.
#-----------------------------------------------------------------------------

SET(REQ_MAJOR_VER 3)
SET(REQ_MINOR_VER 5)
SET(REQ_PATCH_VER 2)

#-----------------------------------------------------------------------------
# Include global definitions for CFS++ development server (host name, etc.)
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/DevelopmentServer.cmake")

#-----------------------------------------------------------------------------
# Give some hints to the user, if CMake is too old.
#-----------------------------------------------------------------------------
IF(CMAKE_MAJOR_VERSION LESS REQ_MAJOR_VER AND
   CMAKE_MINOR_VERSION LESS REQ_MINOR_VER AND
   CMAKE_PATCH_VERSION LESS REQ_PATCH_VER)

MESSAGE(
"
==============================================================================
! ATTENTION: The CMake binaries you are using are too old for building CFS++ !
! But do not worry. There are a number of ways to remedy this situation:     !
==============================================================================

1. On Unix download and run the latest bootstrapper script for CFS++

    $ CFSDSFTP='${CFS_DS_FTP}'
    $ wget $CFSDSFTP/scripts/bootstrap_devel_machine.sh
    $ sudo bash
    $ sh bootstrap_devel_machine.sh

   This will not only install CMake but also all other packages
   required to build CFS++.

2. Run the bootstrapper script from the CFS++ source directory (Unix)

    $ cd CFS_SOURCE_DIR
    $ sudo bash
    $ bash share/scripts/bootstrap_devel_machine.sh

3. Download just the CMake binaries for your platform from Kitware at
   http://cmake.org/cmake/resources/software.html (works on Windows/Linux/Mac)
   You may have to install additional packages by hand! (cf. Developer's
   Manual and share/scripts/bootstrap_devel_machine.sh)

==============================================================================
")

ENDIF(CMAKE_MAJOR_VERSION LESS REQ_MAJOR_VER AND
      CMAKE_MINOR_VERSION LESS REQ_MINOR_VER AND
      CMAKE_PATCH_VERSION LESS REQ_PATCH_VER)

#-----------------------------------------------------------------------------
# Don't build anything unless the version of CMake is high enough.
#-----------------------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(
  VERSION "${REQ_MAJOR_VER}.${REQ_MINOR_VER}.${REQ_PATCH_VER}"
  FATAL_ERROR
)

#-----------------------------------------------------------------------------
# Set CMake policies for newer versions of CMake.
# Enter 'cmake --help-policy CMP000X' to get information about the meaning
# of a certain policy.
#-----------------------------------------------------------------------------
IF(COMMAND CMAKE_POLICY)
  IF(POLICY CMP0038)
    CMAKE_POLICY(SET CMP0038 OLD)
  ENDIF(POLICY CMP0038)
ENDIF(COMMAND CMAKE_POLICY)

#-----------------------------------------------------------------------------
# The name of our project is "CFS".  CMakeLists files in this
# project can refer to the root source directory of the project as
# ${CFS_SOURCE_DIR} and to the root binary directory of the project
# as ${CFS_BINARY_DIR}. We use C, C++ and Fortran as languages.
#-----------------------------------------------------------------------------
PROJECT(CFS C CXX Fortran)

#-----------------------------------------------------------------------------
# Make sure that output of started programs is in English, so that we can
# grep for defined messages.
#-----------------------------------------------------------------------------
SET(ENV{LC_MESSAGES} "en_EN")
SET(ENV{LANG} C)

# Include some standard macros
include(CMakeDependentOption)

#-----------------------------------------------------------------------------
# Include some utility macros (for building ext. libs, etc...)
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/CFS_macros.cmake")

#-----------------------------------------------------------------------------
# CFS version number.  Simply the last modification date
#-----------------------------------------------------------------------------
SET(CFS_VERSION_YEAR "17")
SET(CFS_VERSION_MONTH "12")
SET(CFS_VERSION "${CFS_VERSION_YEAR}.${CFS_VERSION_MONTH}")

#-----------------------------------------------------------------------------
# Version Name. Every major release gets a label/name
#-----------------------------------------------------------------------------
SET(CFS_NAME "AMGme")

colormsg(HIBLUE "Configuration for CFS++ ${CFS_VERSION}: ${CFS_NAME} starting")

#-----------------------------------------------------------------------------
# Set hostname of build host
#-----------------------------------------------------------------------------
SITE_NAME(CFS_BUILD_HOST)
SET(CFS_BUILD_HOST "${CFS_BUILD_HOST}" CACHE INTERNAL "CFS_BUILD_HOST")

#-----------------------------------------------------------------------------
# Set user who builds CFS++, in order to be able to track back who built a
# given binary.
#-----------------------------------------------------------------------------
IF(WIN32 AND NOT CMAKE_CROSSCOMPILING)
  SET(REG_KEY "HKEY_CURRENT_USER\\Software\\Microsoft\\Windows")
  SET(REG_KEY "${REG_KEY}\\CurrentVersion\\Explorer")
  EXECUTE_PROCESS(
    COMMAND reg query "${REG_KEY}" /v "Logon User Name"
    WORKING_DIRECTORY "."
    OUTPUT_VARIABLE OUTPUT
    )

  STRING(REPLACE "\n" ";" OUTPUT "${OUTPUT}")

  FOREACH(line IN LISTS OUTPUT)
    IF(line MATCHES "Logon User Name")
      STRING(REPLACE "Logon User Name" "" line "${line}")
      STRING(REPLACE "REG_SZ" "" line "${line}")
      STRING(STRIP "${line}" CFS_BUILD_USER)
    ENDIF()
  ENDFOREACH()

ELSEIF(UNIX OR MINGW)
  EXECUTE_PROCESS(
    COMMAND id -unr # get user name (possibly DOMAIN\user)
    COMMAND sed s/[A-Z]\\+\\\\// # remove 'DOMAIN\' from 'DOMAIN\user' ('\' must be escaped in cmake)
    WORKING_DIRECTORY "."
    OUTPUT_VARIABLE CFS_BUILD_USER RESULT_VARIABLE res OUTPUT_VARIABLE out
    )
  if(res)
    message("could not determine user: ${res}, ${out}")
  endif(res)
  STRING(STRIP "${CFS_BUILD_USER}" CFS_BUILD_USER)
ENDIF()

SET(CMAKE_MODULE_PATH ${CFS_SOURCE_DIR}/cmake_modules)

#-----------------------------------------------------------------------------
# Set default values according to environment or platform. The user can
# specify his own defaults in $HOME/.cfs_platform_defaults.cmake
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/platform_defaults.cmake")

#-----------------------------------------------------------------------------
# Switch on/off usage of cfsdepscache
#-----------------------------------------------------------------------------
OPTION(CFS_DEPS_PRECOMPILED
  "Use precompiled cfsdeps. If no precompiled version exists, compile, zip and move package into cache."
  ${CFS_DEPS_PRECOMPILED_DEFAULT})

#-----------------------------------------------------------------------------
# Switch on/off usage of custom-built Python
#-----------------------------------------------------------------------------
# external anaconda
OPTION(BUILD_ANACONDA3
  "build anaconda3 from continuum analytics, not necessary for CFS++ itself"
  OFF)

#-----------------------------------------------------------------------------
# Switch on/off external libraries and other options for compile
#-----------------------------------------------------------------------------
OPTION(USE_MESH "Use the mesh file reader for input from ANSYS/GiD."  ${USE_MESH_DEFAULT})
OPTION(USE_COMSOL "Use Comsol .mphtxt meshes for simulation input." ${USE_COMSOL_DEFAULT})
OPTION(USE_GIDPOST "Use the GiD library for simulation output."  ${USE_GIDPOST_DEFAULT})
OPTION(USE_GMSH "Use Gmsh meshes for simulation input." ${USE_GMSH_DEFAULT})
OPTION(USE_GMV "Use General Mesh Viewer (GMV) for simulation in/output."  ${USE_GMV_DEFAULT})
OPTION(USE_UNV "Use Universal File for simulation output." ${USE_UNV_DEFAULT})
# USE_ENSIGHT will we mirrored to USE_VTK which is only used for Ensight
OPTION(USE_ENSIGHT "Support for Ensight Gold file format for input.." ${USE_ENSIGHT_DEFAULT})
OPTION(USE_SGPP    "Enable support for sparse grid interpolation."   ${USE_SGPP_DEFAULT})

IF(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  OPTION(USE_ANSYSRST
    "Use ANSYS .RST format for simulation output."
    ${USE_ANSYSRST_DEFAULT})
ENDIF()

#-----------------------------------------------------------------------------
# There are two options for an xml reader: xerces or libxml2. Depending on
# the switch XmlReader::Parse*() use one of them. When a schema is given
# validation and default values are applied based on DOM, otherwise a faster
# interface is used. Xerces is the established but tends to make compile time
# problems on some systems. XML_READER must be either xerces or libxml2. This
# toggles the hidden USE_XERCES and USE_LIBXML
#-----------------------------------------------------------------------------
SET(XML_READER ${XML_READER_DEFAULT}
  CACHE STRING
  "Choose either 'xerces' or 'libxml2' as xml reader.")

#-----------------------------------------------------------------------------
# Options concerning linear algebra and  direct solver libraries The following
# combinations are implemented
#
#         OPENBLAS MKL NETLIB
# SCHENK    x       x    x
#    MKL    -       x    -
#
# Since ILUPACK uses  MUMPS and no longer Pardiso for  its matching operations
# it should now compile with all  other combinations of BLAS/LAPACK If we want
# to try out the Pardiso matching this would have to be implemented separately
#-----------------------------------------------------------------------------
SET(CFS_BLAS_LAPACK ${CFS_BLAS_LAPACK_DEFAULT} CACHE STRING
  "BLAS/LAPACK implementation. May be one of OPENBLAS, MKL, NETLIB")
SET(CFS_PARDISO "MKL" CACHE STRING 
    "Pardiso implementation for CFS. May be one of MKL or SCHENK. Irrelevant if disabled via USE_PARDISO")
OPTION(USE_PARDISO "Use Pardiso, see also advanced option CFS_PARDISO (MKL or SCHENK)" ${USE_PARDISO_DEFAULT})

OPTION(USE_ARPACK "Use the ARPACK eigensolver library." ${USE_ARPACK_DEFAULT})
OPTION(USE_ILUPACK "Iterative solver package with own algebraic multilevel prconditioners." ${USE_ILUPACK_DEFAULT})
OPTION(USE_SUITESPARSE "CholMod and UMFPACK direct solvers from SuiteSparse." ${USE_SUITESPARSE_DEFAULT})
OPTION(USE_LIS "Suite of OpenMP paralell solvers.." ${USE_LIS_DEFAULT})
OPTION(USE_SUPERLU "Use SuperLU direct solver and ILUs." ${USE_SUPERLU_DEFAULT})
OPTION(USE_PETSC "Suite of MPI parallel linear solver, requires local mpi installation" ${USE_PETSC_DEFAULT})
OPTION(USE_PHIST "phist eigensolver library from sppexa." ${USE_PHIST_DEFAULT})
OPTION(USE_GHOST "ghost kernel library (optionally cuda based) from sppexa/FAU. Required for phist" ${USE_GHOST_DEFAULT})

#-----------------------------------------------------------------------------
# Since MKL Pardiso only works with MKL BLAS/LAPACK, we force that here.
#-----------------------------------------------------------------------------
if(USE_PARDISO AND CFS_PARDISO STREQUAL "MKL")
  set(CFS_BLAS_LAPACK "MKL" CACHE STRING
    "BLAS/LAPACK implementation. May be one of OPENBLAS, MKL, NETLIB, Disable USE_PARDISO or set CFS_PARDSIO to SCHENK if not MKL"
    FORCE)
  # MESSAGE("With USE_PARDISO and CFS_PARDISO=MKL we enforce CFS_BLAS_LAPACK=MKL")
endif(USE_PARDISO AND CFS_PARDISO STREQUAL "MKL")

# Do we have mpi? There are multiple mpi installations, one should stick with one
# but we don't build mpi by cfsdeps to allow usage of tuned hpc mpi instances.
# phist requires ghost, ghost requires hwloc. One can set hwloc also for petsc
if(USE_PETSC)
  set(USE_MPI "ON")
  INCLUDE("cmake_modules/mpi.cmake")
else()
  set(USE_MPI "OFF") 
endif()

# phist requires ghost
if(USE_PHIST)
  set(USE_GHOST ON CACHE BOOL "ghost kernel library (optionally cuda based) from sppexa/FAU. Required for phist" FORCE)
endif()

if(USE_GHOST)
  set(BUILD_HWLOC ON CACHE BOOL "toggled by USE_PHIST and USE_GHOST" FORCE)  
else() 
  set(BUILD_HWLOC OFF CACHE BOOL "toggled by USE_PHIST and USE_GHOST" FORCE)  
endif()

# Disables Optimization in CFS++. Optimization is (and shall) still compiled but
# optimization problems don't run. Furthermore CFS_Optimiztion.xsd will be replaced
# in the build copy of the schema. Meant for deploying CFS++ to the public.
# Note that this option is advanced, so you need to toggle in ccmake to see it.
OPTION(DISABLE_OPTIMIZATION "Disables to run optimization problmes (still compiled)" ${DISABLE_OPTIMIZATION_DEFAULT})

#-----------------------------------------------------------------------------
# set USE_XERCES and USE_LIBXML2 (disjunct!) according to XML_READER choice
IF(NOT XML_READER STREQUAL "xerces" AND NOT XML_READER STREQUAL "libxml2")
  MESSAGE(FATAL_ERROR "Set XML_READER to either 'xerces' or 'libxml2'!")
ENDIF()

IF(XML_READER STREQUAL "xerces")
  SET(USE_XERCES ON CACHE BOOL "toggled by XML_READER" FORCE)
  SET(USE_LIBXML2 OFF CACHE BOOL "toggled by XML_READER" FORCE)
ELSE()
  SET(USE_XERCES OFF CACHE BOOL "toggled by XML_READER" FORCE)
  SET(USE_LIBXML2 ON CACHE BOOL "toggled by XML_READER" FORCE)
ENDIF()  


#-----------------------------------------------------------------------------
# Options which have become so important that they are not optional anymore.
#-----------------------------------------------------------------------------
SET(USE_HDF5 ${USE_HDF5_DEFAULT}
  CACHE INTERNAL "Use HDF5 for simulation in/output.")
SET(USE_METIS ${USE_METIS_DEFAULT}
  CACHE INTERNAL "Use Metis Graph partitioning.")

#-----------------------------------------------------------------------------
# FFTW - Fastest Fourier Transform in the West (C subroutine library) 
#-----------------------------------------------------------------------------
OPTION(USE_FFTW
  "Activate DFT and filtering capabilities of cfstool." 
  ${USE_FFTW_DEFAULT})
 
#-----------------------------------------------------------------------------
# CGNS - CFD General Notation System
#-----------------------------------------------------------------------------
OPTION(USE_CGNS "Enable support for the CFD General Notation System."
  ${USE_CGNS_DEFAULT})

# Since the CGNS library depends on the HDF5 library, we need to force
# USE_HDF5 to ON.
IF(USE_CGNS)
    SET(USE_HDF5 ON CACHE BOOL "Use HDF5 for simulation in/output." FORCE)
ENDIF()

#-----------------------------------------------------------------------------
# CCMIO - STAR-CCM+ I/O library
#-----------------------------------------------------------------------------
OPTION(USE_CCMIO "Enable support for the STAR-CCM+ I/O library."
  ${USE_CCMIO_DEFAULT})

# Since the CCMIO library depends on the ADF library from CGNS, we need to
# force USE_CGNS to ON.
IF(USE_CCMIO AND
    NOT USE_CGNS)
    SET(USE_CGNS ON CACHE BOOL
      "Enable support for the CFD General Notation System." FORCE)
ENDIF()

#-----------------------------------------------------------------------------
# CFXIO - Support for reading ANSYS CFX result files
#-----------------------------------------------------------------------------
IF(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  OPTION(USE_CFXIO "Enable support for reading ANSYS CFX result files."
    ${USE_CFXIO_DEFAULT})
ENDIF()

#-----------------------------------------------------------------------------
# SCPIP optimizer library.
#-----------------------------------------------------------------------------
OPTION(USE_SCPIP "Use the external academic SCPIP optimizer." ${USE_SCPIP_DEFAULT})

#-----------------------------------------------------------------------------
# commercial SNOPT optimizer library.
#-----------------------------------------------------------------------------
OPTION(USE_SNOPT "Use the external commercial SNOPT optimizer." ${USE_SNOPT_DEFAULT})

#-----------------------------------------------------------------------------
# open source IPOPT optimizer library.
#-----------------------------------------------------------------------------
OPTION(USE_IPOPT "Use the external open source IPOPT optimizer." ${USE_IPOPT_DEFAULT})


#-----------------------------------------------------------------------------
# Use cgal toggles the inclusion of the CGAL library.
#-----------------------------------------------------------------------------
OPTION(USE_CGAL
  "Activate fast element mapping using CGAL library."
  ${USE_CGAL_DEFAULT})

#-----------------------------------------------------------------------------
# Use libfbi toggles the inclusion of the Fast Box Intersection library.
# Please note, that either CGAL or libfbi can be used at a time, but not both!
#-----------------------------------------------------------------------------
OPTION(USE_LIBFBI
  "Use fast box intersection library."
  ${USE_LIBFBI_DEFAULT})

IF(USE_CGAL AND USE_LIBFBI)
  MESSAGE(FATAL_ERROR
    "Either USE_CGAL or USE_LIBFBI can be actived, but not both!")
ENDIF()

#-----------------------------------------------------------------------------
# Use FLANN toggles the inclusion of the FLANN library.
#-----------------------------------------------------------------------------
OPTION(USE_FLANN
  "Enable support for FLANN - Fast Library for Approximate Nearest Neighbors."
  ${USE_FLANN_DEFAULT})

OPTION(DEBUG "Build debug version of the code." ${DEBUG_DEFAULT})
OPTION(PROFILING "Set profiling compiler options for Release." ${PROFILING_DEFAULT})
OPTION(COVERAGE "Turn on coverage compier options."  ${COVERAGE_DEFAULT})
OPTION(CFSTOOL "Build CFS++ converter tool (cfstool)."  ${CFSTOOL_DEFAULT})
OPTION(CFSDAT "Build CFS++ data analysis tool (cfsdat)."  ${CFSDAT_DEFAULT})
OPTION(UNIT_TESTS "Build CFS++ unit tests (cfstest)."  ${UNIT_TESTS_DEFAULT})

#-----------------------------------------------------------------------------
# Shall we build HDFView?
#-----------------------------------------------------------------------------
OPTION(BUILD_HDFVIEW
  "Build the HDF 5 file viewer?"
   ${BUILD_HDFVIEW_DEFAULT})

#-----------------------------------------------------------------------------
# Do we want to use explicit or implicit template instantiation
#-----------------------------------------------------------------------------
OPTION(EXPLICIT_TEMPLATE_INSTANTIATION
  "Use explicit template instantiation."
  ${EXPLICIT_TEMPLATE_INSTANTIATION_DEFAULT})

#-----------------------------------------------------------------------------
# Set a sane default value for the CFS++ schema path
#-----------------------------------------------------------------------------
SET(XMLSCHEMA "${CFS_SOURCE_DIR}/share/xml" CACHE STRING
  "Path where the XML schema for the input files is located.")

#-----------------------------------------------------------------------------
# Specify path to CFS++ testsuite. If this path exists, CTest compatible
# test cases will be generated in the current binary directory.
#-----------------------------------------------------------------------------
IF(NOT TESTSUITE_DIR_DEFAULT)
  SET(TESTSUITE_DIR_DEFAULT "TESTSUITE_DIR_NOT_FOUND")
ENDIF(NOT TESTSUITE_DIR_DEFAULT) 

SET(TESTSUITE_DIR "${TESTSUITE_DIR_DEFAULT}" CACHE PATH "Path where the CFS++ testsuite is located.")

IF(NOT TESTSUITE_DIR STREQUAL "TESTSUITE_DIR_NOT_FOUND")
  IF(EXISTS ${TESTSUITE_DIR}/TESTSUIT)
    # In order to perform tests we need the CFSTOOL and HDF5!
    SET(CFSTOOL ON CACHE BOOL
      "Build CFS++ converter tool (cfstool)." FORCE)
    
    SET(USE_HDF5 ON CACHE INTERNAL
      "Use HDF5 for simulation in/output." FORCE)
    
    ADD_SUBDIRECTORY(${TESTSUITE_DIR}
      "${CMAKE_CURRENT_BINARY_DIR}/testsuite")

  ELSE(EXISTS ${TESTSUITE_DIR}/TESTSUIT)
    SET(MSG "Testsuite directory '${TESTSUITE_DIR}' does not contain")
    SET(MSG "${MSG} TESTSUIT.\nPlease specify a proper location!")

    MESSAGE("${MSG}")
  ENDIF(EXISTS ${TESTSUITE_DIR}/TESTSUIT)
ENDIF(NOT TESTSUITE_DIR STREQUAL "TESTSUITE_DIR_NOT_FOUND")

#-----------------------------------------------------------------------------
# Mark some of the options as advanced to clear up the standard selection
#-----------------------------------------------------------------------------
MARK_AS_ADVANCED(CFS_PARDISO)
MARK_AS_ADVANCED(CFS_DEPS_TOCACHE)
MARK_AS_ADVANCED(USE_FFTW)
MARK_AS_ADVANCED(USE_UNV)
MARK_AS_ADVANCED(USE_GMV)
MARK_AS_ADVANCED(USE_MESH)
MARK_AS_ADVANCED(INTEL_OPT_REPORT)
MARK_AS_ADVANCED(CMAKE_BUILD_TYPE)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(EXPLICIT_TEMPLATE_INSTANTIATION)
MARK_AS_ADVANCED(BUILD_HDFVIEW)
MARK_AS_ADVANCED(COVERAGE)
MARK_AS_ADVANCED(PROFILING)
MARK_AS_ADVANCED(XMLSCHEMA)
MARK_AS_ADVANCED(USE_XERCES)
MARK_AS_ADVANCED(USE_LIBXML2)
#MARK_AS_ADVANCED(USE_GHOST)
MARK_AS_ADVANCED(BUILD_HWLOC)
MARK_AS_ADVANCED(DISABLE_OPTIMIZATION)

#-----------------------------------------------------------------------------
# Set some variables according to build type DEBUG / RELEASE
#-----------------------------------------------------------------------------
IF(DEBUG)
  SET(CMAKE_BUILD_TYPE "DEBUG" CACHE STRING
    "Only DEBUG and RELEASE types are valid." FORCE)
  
  #---------------------------------------------------------------------------
  # Check for certain error conditions (in DEBUG mode)
  #---------------------------------------------------------------------------
  SET(CHECK_INDEX 1)
  SET(CHECK_INITIALIZED 1)
  SET(CHECK_MEM_ALLOC 1)
  SET(CHECK_MEMORY 1)

ELSE(DEBUG)
  SET(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING
    "Only DEBUG and RELEASE types are valid." FORCE)
  
  #---------------------------------------------------------------------------
  # Activate expression templates (in RELEASE mode)
  #---------------------------------------------------------------------------
  SET(EXPR_TEMPLATES 1)
ENDIF(DEBUG)

#-----------------------------------------------------------------------------
# Find executables of a few required programs
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/FindPrograms.cmake")

#-----------------------------------------------------------------------------
# Set compiler flags according to output target and debug setting
# ( if necessesary )
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/compiler.cmake")

#-----------------------------------------------------------------------------
# Determine on what UNIX/Linux/Windows distribution we are building on
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/distro.cmake")

#-----------------------------------------------------------------------------
# Determine if some Fortran system libs are installed (blas, lapack, g2c...)
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/CheckFortranRuntime.cmake")

#-----------------------------------------------------------------------------
# Perform some configure checks
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/ConfigureChecks.cmake")

#-----------------------------------------------------------------------------
# Configure Dart testing support.  This should be done before any
# MESSAGE(FATAL_ERROR ...) commands are invoked (???) and after the languages
# for the project have been defined since CMake needs to know about the 
# language-specific build tools.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Switch support for testing off by default. It will be switched on by the
# CTest scripts. These scripts set BUILD_TESTING to ON in the initial CMake
# cache. Cf. ctest_scripts/*.cmake
#-----------------------------------------------------------------------------
OPTION(BUILD_TESTING "Build the testing tree." ${BUILD_TESTING_DEFAULT})

IF(BUILD_TESTING)
  SET(BUILDNAME
    "${CFS_ARCH_STR} ${CFS_CXX_COMPILER_NAME} ${CFS_CXX_COMPILER_VER} ${CMAKE_BUILD_TYPE}"
    CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)

  #---------------------------------------------------------------------------
  # Include CMake module for testing support
  # This adds the option BUILD_TESTING again! -> BUILD_TESTING overwritten?
  #---------------------------------------------------------------------------
  INCLUDE(CTest)
ENDIF(BUILD_TESTING)

#-----------------------------------------------------------------------------
# Enables python testsuite (contains matviz.py tests) used be the optimization group
#-----------------------------------------------------------------------------
OPTION(PYTHON_TESTS "Add testing of python modules (currently stuff from the optimization group)")
  
#-----------------------------------------------------------------------------
# Switch on/off usage of intel optimization report
#-----------------------------------------------------------------------------
OPTION(INTEL_OPT_REPORT "Tells Intel compiler to print out additional information on code optimization during compiler time" ${INTEL_OPT_REPORT_DEFAULT})

# The performance tests are meant for the nightly builds. For some strange reasons
# PERFORMANCE_TESTS doesn't work as it statys always OFF in the test suite?? 
OPTION(PERF_TESTS "Add performance tests to the test-cases")


#-----------------------------------------------------------------------------
# Enables tests of Python scripts via doctest
#-----------------------------------------------------------------------------
CMAKE_DEPENDENT_OPTION(PYTHON_TESTS "Add python tests" OFF
                       "BUILD_TESTING" OFF) # set to off all the time

#-----------------------------------------------------------------------------
# Set proper paths to development servers in testing readme.
#-----------------------------------------------------------------------------
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/ctest_scripts/README.md
  ${CMAKE_CURRENT_BINARY_DIR}/ctest_scripts/README.md
  @ONLY)

#-----------------------------------------------------------------------------
# Set path where shared/static libraries should go
#-----------------------------------------------------------------------------
IF(CFS_ARCH STREQUAL "X86_64")
  SET(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib64/${CFS_ARCH_STR}
    CACHE INTERNAL "Single output directory for building all libraries.")
ELSE(CFS_ARCH STREQUAL "X86_64")
  SET(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}
    CACHE INTERNAL "Single output directory for building all libraries.")
ENDIF(CFS_ARCH STREQUAL "X86_64")

#-----------------------------------------------------------------------------
# Set path where executables should go
#-----------------------------------------------------------------------------
SET(EXECUTABLE_OUTPUT_PATH ${CFS_BINARY_DIR}/bin/${CFS_ARCH_STR}
  CACHE INTERNAL "Single output directory for building all executables.")

#-----------------------------------------------------------------------------
# Create symlinks to lib & bin dirs allowing the build to run on (newer) linux
#-----------------------------------------------------------------------------
if(UNIX)
  # lib dir
  get_filename_component(wd ${LIBRARY_OUTPUT_PATH} DIRECTORY)
  add_custom_target(symlink_lib_dir ALL
     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CFS_ARCH_STR} "LINUX_${CFS_ARCH}"
     WORKING_DIRECTORY ${wd})
  # bin dir
  add_custom_target(symlink_bin_dir ALL
     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CFS_ARCH_STR} "LINUX_${CFS_ARCH}"
     WORKING_DIRECTORY "${CFS_BINARY_DIR}/bin")
endif(UNIX)

#-----------------------------------------------------------------------------
# Create and set a path where temporaries should go
#-----------------------------------------------------------------------------
SET(CFS_TEMP_DIR "${CFS_BINARY_DIR}/tmp")
FILE(MAKE_DIRECTORY "${CFS_TEMP_DIR}")

#-----------------------------------------------------------------------------
# Find CFSDEPS, i.e. gidpost, HDF5, pardiso, arpack...
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/FindCFSDEPS.cmake")

#-----------------------------------------------------------------------------
# Set path of XML schema. If env variable XMLSCHEMA is set, use it!
#-----------------------------------------------------------------------------
SET(MYXMLSCHEMA $ENV{XMLSCHEMA})
IF(NOT ${MYXMLSCHEMA} STREQUAL "")
  SET(XMLSCHEMA ${MYXMLSCHEMA} CACHE STRING
    "Path where the XML schema for the input files is located."
    FORCE)
ENDIF(NOT ${MYXMLSCHEMA} STREQUAL "")

#-----------------------------------------------------------------------------
# Set date of configure run. TODO: make platform independent!!!
#-----------------------------------------------------------------------------
TODAY(CFS_CONF_DATE)

#-----------------------------------------------------------------------------
# Determine URL of subversion repository, global revision of working copy
# and other informations from working copy
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/DetermineWorkingCopyInfos.cmake")

#-----------------------------------------------------------------------------
# Make sure the compiler can find include files from CFS.
#-----------------------------------------------------------------------------
INCLUDE_DIRECTORIES(
  ${Boost_INCLUDE_DIR} 
  ${CFS_BINARY_DIR}/include
  ${CFS_SOURCE_DIR}/source)

IF(USE_CGAL)
  INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIR} ${MPFR_INCLUDE_DIR})
ENDIF(USE_CGAL)

IF(USE_VTK)
  INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIR})
ENDIF(USE_VTK)

#-----------------------------------------------------------------------------
# Make sure the linker can find the CFS base libraries once it is built.
#-----------------------------------------------------------------------------
LINK_DIRECTORIES(${CFSDEPS_LIBRARY_DIR} ${LIBRARY_OUTPUT_PATH})

#-----------------------------------------------------------------------------
# Recurse into the "source" subdirectory.  This does not actually
# cause another cmake executable to run.  The same process will
# walk through the project's entire directory structure.
#-----------------------------------------------------------------------------
ADD_SUBDIRECTORY(source/OLAS/algsys)
ADD_SUBDIRECTORY(source/OLAS/graph)
ADD_SUBDIRECTORY(source/OLAS/precond)
ADD_SUBDIRECTORY(source/OLAS/solver)
ADD_SUBDIRECTORY(source/OLAS/utils)

IF(USE_ARPACK)
  ADD_SUBDIRECTORY(source/OLAS/external/arpack)
ENDIF(USE_ARPACK)

IF(USE_PHIST)
  ADD_SUBDIRECTORY(source/OLAS/external/phist)
ENDIF(USE_PHIST)

IF(USE_ILUPACK)
  ADD_SUBDIRECTORY(source/OLAS/external/ilupack)
ENDIF(USE_ILUPACK)

ADD_SUBDIRECTORY(source/OLAS/external/lapack)

IF(USE_PARDISO)
  ADD_SUBDIRECTORY(source/OLAS/external/pardiso)
ENDIF(USE_PARDISO)

IF(USE_SUITESPARSE)
  ADD_SUBDIRECTORY(source/OLAS/external/umfpack)
  ADD_SUBDIRECTORY(source/OLAS/external/cholmod)
ENDIF(USE_SUITESPARSE)

IF(USE_SUPERLU)
  ADD_SUBDIRECTORY(source/OLAS/external/superlu)
ENDIF(USE_SUPERLU)

IF(USE_LIS)
  ADD_SUBDIRECTORY(source/OLAS/external/lis)
ENDIF(USE_LIS)

IF(USE_PETSC)
  ADD_SUBDIRECTORY(source/OLAS/external/petsc)
ENDIF(USE_PETSC)

ADD_SUBDIRECTORY(source/CoupledPDE)
ADD_SUBDIRECTORY(source/DataInOut)
ADD_SUBDIRECTORY(source/DataInOut/Logging)
ADD_SUBDIRECTORY(source/DataInOut/ParamHandling)

ADD_SUBDIRECTORY(source/DataInOut/SimInOut/AnsysCDB)
IF(USE_MESH)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/AnsysFile)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/internalMesh)
ENDIF(USE_MESH)

IF(USE_COMSOL)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/COMSOL)
ENDIF(USE_COMSOL)

IF(USE_GMV)
 SET(USE_GMV_INPUT ON 
    CACHE INTERNAL
    "Use GMV files as input for simulation. Debugging may not work!")

  SET(USE_GMV_OUTPUT ON 
    CACHE INTERNAL "Use GMV files as output for simulation.")
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/gmv)
ENDIF(USE_GMV)

IF(USE_GMSH)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/gmsh)
ENDIF(USE_GMSH)

ADD_SUBDIRECTORY(source/DataInOut/SimInOut/RefElems)

IF(USE_HDF5)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/hdf5)
ENDIF(USE_HDF5)

IF(USE_GIDPOST)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/GiD)
ENDIF(USE_GIDPOST)

IF(USE_UNV)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/Unverg)
ENDIF(USE_UNV)
 
IF(USE_ENSIGHT)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/VTKBased/Ensight)
  SET(USE_VTK ON  CACHE INTERNAL "Surplus comment for cached variable") 
ELSE()
  SET(USE_VTK OFF CACHE INTERNAL "Surplus comment for cached variable")      
ENDIF()

IF(USE_ANSYSRST)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/AnsysRST)
ENDIF(USE_ANSYSRST)

IF(USE_CGNS)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/CGNS)
ENDIF(USE_CGNS)

ADD_SUBDIRECTORY(source/DataInOut/SimInOut/TextOutput)
ADD_SUBDIRECTORY(source/DataInOut/SimInOut/InfoResultOutput)
IF(NOT MINGW)
IF(NOT CFS_CXX_COMPILER_NAME STREQUAL "ICC")
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/Streaming)
ENDIF()
ENDIF()
ADD_SUBDIRECTORY(source/DataInOut/ScatteredDataInOut)

ADD_SUBDIRECTORY(source/Domain)
ADD_SUBDIRECTORY(source/Domain/CoefFunction)
ADD_SUBDIRECTORY(source/Domain/CoordinateSystems)
ADD_SUBDIRECTORY(source/Domain/ElemMapping)
ADD_SUBDIRECTORY(source/Domain/Mesh)
ADD_SUBDIRECTORY(source/Domain/Results)

ADD_SUBDIRECTORY(source/Driver)
ADD_SUBDIRECTORY(source/Driver/TimeSchemes)

ADD_SUBDIRECTORY(source/FeBasis)

ADD_SUBDIRECTORY(source/Forms)
#here we have to discussed
#as the forms are almost completely templates
#we have just headers right now. but perhaps this will change...
#ADD_SUBDIRECTORY(source/Forms/Operators)
ADD_SUBDIRECTORY(source/Forms/LinForms)
ADD_SUBDIRECTORY(source/Forms/BiLinForms)

ADD_SUBDIRECTORY(source/General)
ADD_SUBDIRECTORY(source/Materials)
ADD_SUBDIRECTORY(source/MatVec)
ADD_SUBDIRECTORY(source/ODESolve)
ADD_SUBDIRECTORY(source/PDE)
ADD_SUBDIRECTORY(source/PDE/LatticeBoltzmannSolver)
ADD_SUBDIRECTORY(source/Utils)
ADD_SUBDIRECTORY(source/Utils/mathParser)
ADD_SUBDIRECTORY(source/Utils/EvalIntegrals)
ADD_SUBDIRECTORY(source/Optimization)
ADD_SUBDIRECTORY(source/Optimization/Design)
ADD_SUBDIRECTORY(source/Optimization/Optimizer)
ADD_SUBDIRECTORY(source/main)

IF(CFSTOOL)
  ADD_SUBDIRECTORY(source/cfstool)
  #  SET(USE_FFTW ON CACHE BOOL 
  #  "Activate DFT and filtering capabilities of cfstool." FORCE)
ENDIF(CFSTOOL)

IF(CFSDAT)
  ADD_SUBDIRECTORY(source/cfsdat)
ENDIF(CFSDAT)

if(UNIT_TESTS)
  add_subdirectory(source/unittests)
endif(UNIT_TESTS)

#-----------------------------------------------------------------------------
# Build the actual headers with defines from templates.
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/configure_def_headers.cmake")

#-----------------------------------------------------------------------------
# Add targets for building documentation (see 'make help')
#-----------------------------------------------------------------------------
ADD_SUBDIRECTORY(share/doc)

#-----------------------------------------------------------------------------
# Add targets for scripts
#-----------------------------------------------------------------------------
ADD_SUBDIRECTORY(share/scripts)

#-----------------------------------------------------------------------------
# Add targets for mkmesh
#-----------------------------------------------------------------------------
ADD_SUBDIRECTORY(share/mkmesh)

#-----------------------------------------------------------------------------
# Add targets for XML schemata and XSL style sheets
#-----------------------------------------------------------------------------
ADD_SUBDIRECTORY(share/xml)
ADD_SUBDIRECTORY(share/xsl)
