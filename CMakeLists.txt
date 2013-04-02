#-----------------------------------------------------------------------------
# Set the version number of CMake which is required to build CFS++.
#-----------------------------------------------------------------------------
SET(REQ_MAJOR_VER 2)
SET(REQ_MINOR_VER 8)
SET(REQ_PATCH_VER 8)

#-----------------------------------------------------------------------------
# Give some hints to the user, if CMake is too old.
#-----------------------------------------------------------------------------
IF(CMAKE_MAJOR_VERSION LESS REQ_MAJOR_VER OR
   CMAKE_MINOR_VERSION LESS REQ_MINOR_VER OR
   CMAKE_PATCH_VERSION LESS REQ_PATCH_VER)

MESSAGE(
"
================================================================================
! ATTENTION: The CMake binaries you are using are too old for building CFS++!  !
! But do not worry. There are a number of ways to remedy this situation:       !
================================================================================

1. On Unix download and run the latest bootstrapper script for CFS++

    $ LSE17FTP='ftp://lse17.e-technik.uni-erlangen.de:40065'
    $ wget $LSE17FTP/scripts/bootstrap_devel_machine.sh
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
   You may have to install additional packages by hand! (cf. Developer's Manual)

================================================================================
")

ENDIF(CMAKE_MAJOR_VERSION LESS REQ_MAJOR_VER OR
      CMAKE_MINOR_VERSION LESS REQ_MINOR_VER OR
      CMAKE_PATCH_VERSION LESS REQ_PATCH_VER)

#-----------------------------------------------------------------------------
# Don't build anything unless the version of CMake is high enough.
#-----------------------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(
  VERSION "${REQ_MAJOR_VER}.${REQ_MINOR_VER}.${REQ_PATCH_VER}"
  FATAL_ERROR
)

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

#-----------------------------------------------------------------------------
# Include some utility macros (for building ext. libs, etc...) 
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/CFS_macros.cmake")

#-----------------------------------------------------------------------------
# CFS version number.  Simply the last modification date
#-----------------------------------------------------------------------------
SET(CFS_VERSION_YEAR "13")
SET(CFS_VERSION_MONTH "02")
SET(CFS_VERSION
  "${CFS_VERSION_YEAR}.${CFS_VERSION_MONTH}")

#-----------------------------------------------------------------------------
# Version Name. Every major release gets a label/name
#-----------------------------------------------------------------------------
SET(CFS_NAME "Cocky Coriolis")

colormsg(HIBLUE "Configuration for CFS++ ${CFS_VERSION}: ${CFS_NAME} starting")

#-----------------------------------------------------------------------------
# Set hostname of build host
#-----------------------------------------------------------------------------
SITE_NAME(CFS_BUILD_HOST)
SET(CFS_BUILD_HOST "${CFS_BUILD_HOST}"
  CACHE INTERNAL "CFS_BUILD_HOST")

#-----------------------------------------------------------------------------
# Set user who builds CFS++, in order to track back who built a given binary.
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
    COMMAND id -u -nr 
    WORKING_DIRECTORY "."
    OUTPUT_VARIABLE CFS_BUILD_USER
    )
  STRING(STRIP "${CFS_BUILD_USER}" CFS_BUILD_USER)
ENDIF()

SET(CMAKE_MODULE_PATH ${CFS_SOURCE_DIR}/cmake_modules)

#-----------------------------------------------------------------------------
# Set default values according to environment or platform
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/platform_defaults.cmake")

#-----------------------------------------------------------------------------
# Switch on/off usage of custom-built Python
#-----------------------------------------------------------------------------
OPTION(USE_PYTHON
  "Build our own Python for CFS++." 
  ${USE_PYTHON_DEFAULT})

#-----------------------------------------------------------------------------
# Switch on/off external libraries and other options for compile
#-----------------------------------------------------------------------------
OPTION(USE_MESH
  "Use the mesh file reader for input from ANSYS/GiD."
  ${USE_MESH_DEFAULT})
OPTION(USE_COMSOL
  "Use Comsol .mphtxt meshes for simulation input."
  ${USE_COMSOL_DEFAULT})
OPTION(USE_GIDPOST
  "Use the GiD library for simulation output."
  ${USE_GIDPOST_DEFAULT})
OPTION(USE_GMSH
  "Use Gmsh meshes for simulation input."
  ${USE_GMSH_DEFAULT})
OPTION(USE_GMV
  "Use General Mesh Viewer (GMV) for simulation in/output."
  ${USE_GMV_DEFAULT})
OPTION(USE_UNV
  "Use Universal File for simulation output."
  ${USE_UNV_DEFAULT})

IF(NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  OPTION(USE_ANSYSRST
    "Use ANSYS .RST format for simulation output."
    ${USE_ANSYSRST_DEFAULT})
ENDIF()

#-----------------------------------------------------------------------------
# Options concerning linear algebra and direct solver libraries
# The following combinations are implemented
#
#         ACML GOTO MKL NETLIB
# SCHENK   x    x    x    x
#    MKL   -    -    x    -
#
# Since ILUPACK uses MUMPS and no longer Pardiso for its matching operations
# it should now compile with all other combinations of BLAS/LAPACK
# If we want to try out the Pardiso matching this would have to be 
# implemented separately
#-----------------------------------------------------------------------------
SET(CFS_BLAS_LAPACK ${CFS_BLAS_LAPACK_DEFAULT}
  CACHE STRING
  "BLAS/LAPACK implementation for CFS. May be one of ACML, APPLE, GOTO, MKL, NETLIB")

SET(CFS_PARDISO  ${CFS_PARDISO_DEFAULT}
  CACHE STRING
  "Pardiso implementation for CFS. May be one of MKL or SCHENK.")

OPTION(USE_PARDISO "Use Pardiso" ${USE_PARDISO_DEFAULT})
OPTION(USE_ARPACK "Use the ARPACK eigensolver library." ${USE_ARPACK_DEFAULT})
OPTION(USE_ILUPACK "Use Preconditioners from ILUPACK." ${USE_ILUPACK_DEFAULT})
OPTION(USE_SUITESPARSE "Use CholMod and UMFPACK solvers from SuiteSparse." ${USE_SUITESPARSE_DEFAULT})
OPTION(USE_LIS "Use (L)ibrary of (I)terative (S)olvers." ${USE_LIS_DEFAULT})
OPTION(USE_SUPERLU "Use SuperLU direct solver and ILUs." ${USE_SUPERLU_DEFAULT})

IF(USE_PARDISO OR
   USE_ARPACK OR
   USE_ILUPACK OR
   USE_SUITESPARSE)
  SET(USE_BLAS ON CACHE BOOL
    "Use the Basic Linear Algebra Subroutines library." FORCE)
  SET(USE_LAPACK ON CACHE BOOL
    "Use the LAPACK linear algebra library." FORCE)
ELSE(USE_PARDISO OR
     USE_ARPACK OR
     USE_ILUPACK OR
     USE_SUITESPARSE)
  OPTION(USE_LAPACK
    "Use the LAPACK linear algebra library."
    ${USE_LAPACK_DEFAULT})
  OPTION(USE_BLAS
    "Use the Basic Linear Algebra Subroutines library."
    ${USE_BLAS_DEFAULT})
ENDIF(USE_PARDISO OR
      USE_ARPACK OR
      USE_ILUPACK OR
      USE_SUITESPARSE)

#-----------------------------------------------------------------------------
# Since MKL Pardiso only works with MKL BLAS/LAPACK, we force that here.
#-----------------------------------------------------------------------------
IF(USE_PARDISO AND CFS_PARDISO STREQUAL "MKL")
  SET(CFS_BLAS_LAPACK "MKL" CACHE STRING
    "BLAS/LAPACK implementation for CFS. May be one of ACML, GOTO, MKL, NETLIB"
    FORCE)
ENDIF(USE_PARDISO AND CFS_PARDISO STREQUAL "MKL")

#-----------------------------------------------------------------------------
# Options which have become so important that they are not optional anymore
#-----------------------------------------------------------------------------
SET(USE_HDF5 ${USE_HDF5_DEFAULT}
  CACHE INTERNAL "Use HDF5 for simulation in/output.")
SET(USE_XERCES ${USE_XERCES_DEFAULT}
  CACHE INTERNAL "Use the Xerces-C XML parser library.")
SET(USE_METIS ${USE_METIS_DEFAULT}
  CACHE INTERNAL "Use Metis Graph partitioning.")

#-----------------------------------------------------------------------------
# CGNS - CFD General Notation System
#-----------------------------------------------------------------------------
OPTION(USE_CGNS "Enable support for the CFD General Notation System."
  ${USE_CGNS_DEFAULT})

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
OPTION(USE_SCPIP "Use the external SCPIP optimizer." ${USE_SCPIP_DEFAULT})

#-----------------------------------------------------------------------------
# Yet another optimization library.
#-----------------------------------------------------------------------------
OPTION(USE_IPOPT "Use the external IPOPT optimizer." ${USE_IPOPT_DEFAULT})

#-----------------------------------------------------------------------------
# SnOpt optimizer library.
#-----------------------------------------------------------------------------
OPTION(USE_SNOPT "Use the external SnOpt optimizer." ${USE_SNOPT_DEFAULT})

#-----------------------------------------------------------------------------
# KNITRO is a commercial optimizer, requiring a license file/server
#-----------------------------------------------------------------------------
OPTION(USE_KNITRO "Use the external KNITRO optimizer." FALSE)

#-----------------------------------------------------------------------------
# Use cgal toggles the inclusion of the CGAL library.
#-----------------------------------------------------------------------------
OPTION(USE_CGAL
  "Activate fast element mapping using CGAL libray."
  ${USE_CGAL_DEFAULT})

OPTION(DEBUG "Build debug version of the code." ${DEBUG_DEFAULT})
OPTION(PROFILING "Turn on compiler profiling (DEBUG/RELEASE)." ${PROFILING_DEFAULT})
OPTION(COVERAGE "Turn on code coverage output." ${COVERAGE_DEFAULT})

OPTION(CFSTOOL "Build CFS++ converter tool (cfstool)." ${CFSTOOL_DEFAULT})

#-----------------------------------------------------------------------------
# Shall we build ParaView?
#-----------------------------------------------------------------------------
OPTION(BUILD_PARAVIEW
  "Build the ParaView postprocessor with HDF5 and OpenFOAM extensions?."
  ${BUILD_PARAVIEW_DEFAULT})

#-----------------------------------------------------------------------------
# Shall we build HDFView?
#-----------------------------------------------------------------------------
OPTION(BUILD_HDFVIEW
  "Build the HDF 4/5 file viewer?"
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

SET(TESTSUITE_DIR "${TESTSUITE_DIR_DEFAULT}" CACHE STRING
  "Path where the CFS++ testsuite is located.")

IF(NOT TESTSUITE_DIR STREQUAL "TESTSUITE_DIR_NOT_FOUND")
  IF(EXISTS ${TESTSUITE_DIR}/TESTSUIT)
    # In order to perform tests we need the CFSTOOL and HDF5!
    SET(CFSTOOL ON CACHE BOOL
      "Build CFS++ converter tool (cfstool)." FORCE)
    
    SET(USE_HDF5 ON CACHE INTERNAL
      "Use HDF5 for simulation in/output." FORCE)
    
    ADD_SUBDIRECTORY(${TESTSUITE_DIR} "${CMAKE_CURRENT_BINARY_DIR}/testsuite")
  ELSE(EXISTS ${TESTSUITE_DIR}/TESTSUIT)
    MESSAGE("Testsuite directory '${TESTSUITE_DIR}' does not contain TESTSUIT. "
      "Please specify a proper location!")
  ENDIF(EXISTS ${TESTSUITE_DIR}/TESTSUIT)
ENDIF(NOT TESTSUITE_DIR STREQUAL "TESTSUITE_DIR_NOT_FOUND")

#-----------------------------------------------------------------------------
# Mark some of the options as advanced to clear up the standard selection
#-----------------------------------------------------------------------------
MARK_AS_ADVANCED(USE_UNV)
MARK_AS_ADVANCED(USE_GMV)
MARK_AS_ADVANCED(USE_MESH)
MARK_AS_ADVANCED(USE_XERCES)
MARK_AS_ADVANCED(CMAKE_BUILD_TYPE)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(EXPLICIT_TEMPLATE_INSTANTIATION)
MARK_AS_ADVANCED(BUILD_PARAVIEW)
MARK_AS_ADVANCED(BUILD_HDFVIEW)
MARK_AS_ADVANCED(COVERAGE)
MARK_AS_ADVANCED(PROFILING)
MARK_AS_ADVANCED(XMLSCHEMA)

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
# Determine if some Fortran system libs are installed (blas, lapack, g2c...)
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/CheckFortranRuntime.cmake")

#-----------------------------------------------------------------------------
# Determine on what UNIX/Linux/Windows distribution we are building on
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/distro.cmake")

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
OPTION(BUILD_TESTING "Build the testing tree." OFF)

IF(BUILD_TESTING)
  SET(BUILDNAME
    "${CFS_ARCH_STR} ${CFS_CXX_COMPILER_NAME} ${CFS_CXX_COMPILER_VER} ${CMAKE_BUILD_TYPE}"
    CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)

  #---------------------------------------------------------------------------
  # Include CMake module for testing support
  #---------------------------------------------------------------------------
  INCLUDE(CTest)
ENDIF(BUILD_TESTING)


#-----------------------------------------------------------------------------
# Set path where shared/static libraries should go
#-----------------------------------------------------------------------------
IF(CFS_ARCH STREQUAL "X86_64")
  SET(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib64/${CFS_ARCH_STR} CACHE INTERNAL
    "Single output directory for building all libraries.")
ELSE(CFS_ARCH STREQUAL "X86_64")
  SET(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR} CACHE INTERNAL
    "Single output directory for building all libraries.")
ENDIF(CFS_ARCH STREQUAL "X86_64")

#-----------------------------------------------------------------------------
# Set path where executables should go
#-----------------------------------------------------------------------------
SET(EXECUTABLE_OUTPUT_PATH ${CFS_BINARY_DIR}/bin/${CFS_ARCH_STR} CACHE INTERNAL
  "Single output directory for building all executables.")


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
  INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIR})
ENDIF(USE_CGAL)

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

IF(USE_ILUPACK)
  ADD_SUBDIRECTORY(source/OLAS/external/ilupack)
ENDIF(USE_ILUPACK)

IF(USE_LAPACK)
  ADD_SUBDIRECTORY(source/OLAS/external/lapack)
ENDIF(USE_LAPACK)

IF(USE_PARDISO)
  ADD_SUBDIRECTORY(source/OLAS/external/pardiso)
ENDIF(USE_PARDISO)

IF(USE_SUITESPARSE)
  ADD_SUBDIRECTORY(source/OLAS/external/umfpack)
  ADD_SUBDIRECTORY(source/OLAS/external/cholmod)
ENDIF(USE_SUITESPARSE)

IF(USE_LIS)
  ADD_SUBDIRECTORY(source/OLAS/external/lis)
ENDIF(USE_LIS)

ADD_SUBDIRECTORY(source/CoupledPDE)
ADD_SUBDIRECTORY(source/DataInOut)
ADD_SUBDIRECTORY(source/DataInOut/DefineFiles)
ADD_SUBDIRECTORY(source/DataInOut/Logging)
ADD_SUBDIRECTORY(source/DataInOut/ParamHandling)

IF(USE_MESH)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/AnsysFile)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/internalMesh)
ENDIF(USE_MESH)

IF(USE_COMSOL)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/COMSOL)
ENDIF(USE_COMSOL)

IF(USE_GMV)
 SET(USE_GMV_INPUT ON 
    CACHE INTERNAL "Use GMV files as input for simulation. Attention! Debugging may not work!")

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
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/xdmf)
ENDIF(USE_HDF5)

IF(USE_GIDPOST)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/GiD)
ENDIF(USE_GIDPOST)

IF(USE_UNV)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/Unverg)
ENDIF(USE_UNV)

IF(USE_ANSYSRST)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/AnsysRST)
ENDIF(USE_ANSYSRST)

ADD_SUBDIRECTORY(source/DataInOut/SimInOut/TextOutput)
ADD_SUBDIRECTORY(source/DataInOut/SimInOut/InfoResultOutput)
IF(NOT MINGW)
  ADD_SUBDIRECTORY(source/DataInOut/SimInOut/Streaming)
ENDIF(NOT MINGW)

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
ADD_SUBDIRECTORY(source/Utils)
ADD_SUBDIRECTORY(source/Utils/mathParser)
ADD_SUBDIRECTORY(source/Utils/EvalIntegrals)
#ADD_SUBDIRECTORY(source/Optimization/Optimizer)
#ADD_SUBDIRECTORY(source/Optimization/Design)
#ADD_SUBDIRECTORY(source/Optimization)
ADD_SUBDIRECTORY(source/main)

IF(CFSTOOL)
  ADD_SUBDIRECTORY(source/cfstool)
ENDIF(CFSTOOL)

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
# Add targets for XML schemata and XSL style sheets
#-----------------------------------------------------------------------------
ADD_SUBDIRECTORY(share/xml)
ADD_SUBDIRECTORY(share/xsl)
