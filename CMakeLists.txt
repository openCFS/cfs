#-----------------------------------------------------------------------------
# Don't build anything unless the version of CMake is high enough.
#-----------------------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.7 FATAL_ERROR)

#-----------------------------------------------------------------------------
# Set CMake policies for newer versions of CMake.
# Enter 'cmake --help-policy CMP000X' to get information about the meaning
# of a certain policy.
#-----------------------------------------------------------------------------
IF(COMMAND CMAKE_POLICY)
  # This policy affects how libraries whose full paths are NOT known are
  # found at link time, but was created due to a change in how CMake
  # deals with libraries whose full paths are known.
  IF(POLICY CMP0003)
    CMAKE_POLICY(SET CMP0003 NEW)
  ENDIF(POLICY CMP0003)
  

  # Issue an error if no CMakeLists.txt exists in a subdir referenced
  # by SUBDIRS or ADD_SUBDIRECTORY 
  IF(POLICY CMP0014)
    CMAKE_POLICY(SET CMP0014 NEW)
  ENDIF(POLICY CMP0014)
ENDIF(COMMAND CMAKE_POLICY)

#-----------------------------------------------------------------------------
# The name of our project is "CFS".  CMakeLists files in this
# project can refer to the root source directory of the project as
# ${CFS_SOURCE_DIR} and to the root binary directory of the project
# as ${CFS_BINARY_DIR}. We use C, C++ and Fortran as languages.
#-----------------------------------------------------------------------------
PROJECT(CFS C CXX Fortran)

#-----------------------------------------------------------------------------
# CFS version number.  Simply the last modification date
#-----------------------------------------------------------------------------
SET(CFS_VERSION_YEAR "11")
SET(CFS_VERSION_MONTH "06")
SET(CFS_VERSION
  "${CFS_VERSION_YEAR}.${CFS_VERSION_MONTH}")

#-----------------------------------------------------------------------------
# Version Name. Every major release gets a label/name
#-----------------------------------------------------------------------------
SET(CFS_NAME "Querulous Quokka")

#-----------------------------------------------------------------------------
# Set hostname of build host
#-----------------------------------------------------------------------------
SITE_NAME(CFS_BUILD_HOST)
SET(CFS_BUILD_HOST "${CFS_BUILD_HOST}"
  CACHE INTERNAL "CFS_BUILD_HOST")

#-----------------------------------------------------------------------------
# Set user who builds CFS++
#-----------------------------------------------------------------------------
SET(CFS_BUILD_USER "$ENV{USER}"
    CACHE INTERNAL "CFS_BUILD_USER")

#-----------------------------------------------------------------------------
# Set default values according to environment or platform
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/platform_defaults.cmake")

#-----------------------------------------------------------------------------
# Switch on/off external libraries and other options for compile
#-----------------------------------------------------------------------------
OPTION(USE_TCL
  "Use Tool Command Language (TCL) for scripting."
  ${USE_TCL_DEFAULT})
OPTION(USE_PYTHON
  "Use Python for scripting."
  ${USE_PYTHON_DEFAULT})
OPTION(USE_MESH
  "Use the mesh file reader for input from ANSYS/GiD."
  ${USE_MESH_DEFAULT})
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
OPTION(USE_ANSYSRST
  "Use ANSYS .RST format for simulation output."
  ${USE_ANSYSRST_DEFAULT})
OPTION(USE_MPCCI
  "Enable support for the Mesh parallel Code Coupling Interface."
  ${USE_MPCCI_DEFAULT})

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
  "BLAS/LAPACK implementation for CFS. May be one of ACML, GOTO, MKL, NETLIB")

SET(CFS_PARDISO  ${CFS_PARDISO_DEFAULT}
  CACHE STRING
  "Pardiso implementation for CFS. May be one of MKL or SCHENK.")

OPTION(USE_PARDISO "Use Pardiso" ${USE_PARDISO_DEFAULT})
OPTION(USE_ARPACK "Use the ARPACK eigensolver library." ${USE_ARPACK_DEFAULT})
OPTION(USE_ILUPACK "Use Preconditioners from ILUPACK." ${USE_ILUPACK_DEFAULT})
OPTION(USE_CHOLMOD "Use the CholMod solver." ${USE_CHOLMOD_DEFAULT})

IF(USE_PARDISO OR
   USE_ARPACK OR
   USE_ILUPACK OR
   USE_CHOLMOD)
  SET(USE_BLAS ON CACHE BOOL
    "Use the Basic Linear Algebra Subroutines library." FORCE)
  SET(USE_LAPACK ON CACHE BOOL
    "Use the LAPACK linear algebra library." FORCE)
ELSE(USE_PARDISO OR
     USE_ARPACK OR
     USE_ILUPACK OR
     USE_CHOLMOD)
  OPTION(USE_LAPACK
    "Use the LAPACK linear algebra library."
    ${USE_LAPACK_DEFAULT})
  OPTION(USE_BLAS
    "Use the Basic Linear Algebra Subroutines library."
    ${USE_BLAS_DEFAULT})
ENDIF(USE_PARDISO OR
      USE_ARPACK OR
      USE_ILUPACK OR
      USE_CHOLMOD)

IF(CFS_PARDISO STREQUAL "MKL")
  SET(CFS_BLAS_LAPACK "MKL" CACHE STRING
    "BLAS/LAPACK implementation for CFS. May be one of ACML, GOTO, MKL, NETLIB"
    FORCE)
ENDIF(CFS_PARDISO STREQUAL "MKL")

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
# Use interpolation toggles the inclusion of the CGAL library.
#-----------------------------------------------------------------------------
OPTION(USE_INTERPOLATION
  "Activate Interpolation between grids."
  ${USE_INTERPOLATION_DEFAULT})

#-----------------------------------------------------------------------------
# The USE_OPENMP option triggers the usage of OpenMP versions of external
# libraries and the compilation of CFS++ with OpenMP compiler switches.
#-----------------------------------------------------------------------------
OPTION(USE_OPENMP
  "Turn on support for OpenMP. Needs GCC >= 4.2 or Intel C++."
  ${USE_OPENMP_DEFAULT})
OPTION(DEBUG "Build debug version of the code." ${DEBUG_DEFAULT})
OPTION(PROFILING "Turn on compiler profiling (DEBUG/RELEASE)." ${PROFILING_DEFAULT})
OPTION(COVERAGE "Turn on code coverage output." ${COVERAGE_DEFAULT})

OPTION(CFSTOOL "Build CFS++ converter tool (cfstool)." ${CFSTOOL_DEFAULT})
OPTION(CPLREADER "Build cplreader for fluid/acoustic coupling."
  ${CPLREADER_DEFAULT})

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
    # In order to perform tests we need the CONVERTER and HDF5!
    SET(CFSTOOL ON CACHE BOOL
      "Build CFS++ converter tool (cfstool)." FORCE)
    
    SET(USE_HDF5 ON CACHE INTERNAL
      "Use HDF5 for simulation in/output." FORCE)
    
    SUBDIRS(${TESTSUITE_DIR})
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
MARK_AS_ADVANCED(USE_MPCCI)
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
# Set compiler flags according to output target and debug setting
# ( if necessesary )
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/compiler.cmake")

#-----------------------------------------------------------------------------
# Determine if some Fortran system libs are installed (blas, lapack, g2c...)
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/FindFortranLibs.cmake")

#-----------------------------------------------------------------------------
# Determine on what UNIX/Linux/Windows distribution we are building on
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/distro.cmake")

#-----------------------------------------------------------------------------
# Include some utility macros (for building ext. libs, etc...) 
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/CFS_macros.cmake")

#-----------------------------------------------------------------------------
# Find executables of a few required programs
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/FindPrograms.cmake")

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
IF(NOT CFSDEPS_FOUND)
  MESSAGE(FATAL_ERROR "CFSDEPS could not be found!")
ENDIF(NOT CFSDEPS_FOUND)

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
EXEC_PROGRAM("date"
  ARGS +%D
  OUTPUT_VARIABLE CFS_CONF_DATE
  RETURN_VALUE RETVAL)

#-----------------------------------------------------------------------------
# Determine URL of subversion repository
#-----------------------------------------------------------------------------
IF(EXISTS "${CFS_SOURCE_DIR}/.git")
  EXEC_PROGRAM("cd ${CFS_SOURCE_DIR} && git svn info | grep URL | cut -d' ' -f2"
    ARGS 
    OUTPUT_VARIABLE CFS_SUBVERSION_REPOS
    RETURN_VALUE RETVAL)
ELSE(EXISTS "${CFS_SOURCE_DIR}/.git")
  FILE(READ "${CFS_SOURCE_DIR}/.svn/entries" SVNSTR)
  STRING(REGEX MATCH "[a-z+]+://[a-zA-Z0-9_+-/@]+:[a-zA-Z0-9_+-/@]+" CFS_SUBVERSION_REPOS "${SVNSTR}")
ENDIF(EXISTS "${CFS_SOURCE_DIR}/.git")

#-----------------------------------------------------------------------------
# Determine global revision of working copy
#-----------------------------------------------------------------------------
EXECUTE_PROCESS(
  COMMAND ${PERL} ${CFS_DEPS_ROOT}/utils/perl/cfsdepsrev.pl
  ${CFS_SOURCE_DIR} ${CFS_BINARY_DIR}/tmp
  OUTPUT_VARIABLE CFS_SUBVERSION_REV
  ERROR_VARIABLE CFS_SUBVERSION_REV
  RESULT_VARIABLE RETVAL)

IF(RETVAL STREQUAL 0)
  STRING(REGEX MATCH "WORKING_COPY_REVISION .*"
    CFS_SUBVERSION_REV "${CFS_SUBVERSION_REV}")
  STRING(REGEX REPLACE "WORKING_COPY_REVISION " ""
    CFS_SUBVERSION_REV "${CFS_SUBVERSION_REV}")
ELSE(RETVAL STREQUAL 0)
  SET(CFS_SUBVERSION_REV "N/A")
  MESSAGE("A problem occurred while trying to determine revision of "
    "CFS++ Subversion repository: ${CFS_SUBVERSION_REV}")
ENDIF(RETVAL STREQUAL 0)


EXECUTE_PROCESS(
  COMMAND ${PERL} ${CFS_DEPS_ROOT}/utils/perl/cfsdepsmodif.pl
  ${CFS_SOURCE_DIR} ${CFS_BINARY_DIR}/tmp
  OUTPUT_VARIABLE CFS_SUBVERSION_MODIF
  ERROR_VARIABLE CFS_SUBVERSION_MODIF
  RESULT_VARIABLE RETVAL)

IF(RETVAL STREQUAL 0)
  STRING(REGEX MATCH "WORKING_COPY_MODIFIED .*"
    CFS_SUBVERSION_MODIF "${CFS_SUBVERSION_MODIF}")
  STRING(REGEX REPLACE "WORKING_COPY_MODIFIED " ""
    CFS_SUBVERSION_MODIF "${CFS_SUBVERSION_MODIF}")
  IF(CFS_SUBVERSION_MODIF MATCHES "modified" OR
      CFS_SUBVERSION_MODIF MATCHES "git")
    SET(CFS_SUBVERSION_REV "${CFS_SUBVERSION_REV} ${CFS_SUBVERSION_MODIF}")
  ENDIF(CFS_SUBVERSION_MODIF MATCHES "modified" OR
    CFS_SUBVERSION_MODIF MATCHES "git")
ELSE(RETVAL STREQUAL 0)
  SET(CFS_SUBVERSION_REV "N/A")
  MESSAGE("A problem occurred while trying to check modified state of "
    "CFS++ Subversion repository: ${CFS_SUBVERSION_MODIFIED}")
ENDIF(RETVAL STREQUAL 0)

#-----------------------------------------------------------------------------
# Make sure the compiler can find include files from CFS.
#-----------------------------------------------------------------------------
INCLUDE_DIRECTORIES(
  ${Boost_INCLUDE_DIR} 
  ${CFS_BINARY_DIR}/include
  ${CFS_SOURCE_DIR}/source)

IF(USE_INTERPOLATION)
  INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIR})
ENDIF(USE_INTERPOLATION)

#-----------------------------------------------------------------------------
# Make sure the linker can find the CFS base libraries once it is built.
#-----------------------------------------------------------------------------
LINK_DIRECTORIES(${CFSDEPS_LIBRARY_DIR} ${LIBRARY_OUTPUT_PATH})

#-----------------------------------------------------------------------------
# Recurse into the "source" subdirectory.  This does not actually
# cause another cmake executable to run.  The same process will
# walk through the project's entire directory structure.
#-----------------------------------------------------------------------------
SUBDIRS(source/OLAS/algsys)
SUBDIRS(source/OLAS/graph)
SUBDIRS(source/OLAS/precond)
SUBDIRS(source/OLAS/solver)
SUBDIRS(source/OLAS/utils)

IF(USE_ARPACK)
  SUBDIRS(source/OLAS/external/arpack)
ENDIF(USE_ARPACK)

IF(USE_ILUPACK)
  SUBDIRS(source/OLAS/external/ilupack)
ENDIF(USE_ILUPACK)

IF(USE_LAPACK)
  SUBDIRS(source/OLAS/external/lapack)
ENDIF(USE_LAPACK)

IF(USE_PARDISO)
  SUBDIRS(source/OLAS/external/pardiso)
ENDIF(USE_PARDISO)

IF(USE_CHOLMOD)
  SUBDIRS(source/OLAS/external/cholmod)
ENDIF(USE_CHOLMOD)

SUBDIRS(source/CoupledPDE)
SUBDIRS(source/DataInOut)
SUBDIRS(source/DataInOut/DefineFiles)
SUBDIRS(source/DataInOut/Logging)
SUBDIRS(source/DataInOut/ParamHandling)

SUBDIRS(source/DataInOut/Scripting)

IF(USE_SCRIPTING_PYTHON)
  SUBDIRS(source/DataInOut/Scripting/python)
ENDIF(USE_SCRIPTING_PYTHON)

IF(USE_SCRIPTING_TCL)
  SUBDIRS(source/DataInOut/Scripting/tcl)
ENDIF(USE_SCRIPTING_TCL)


IF(USE_MESH)
  SUBDIRS(source/DataInOut/SimInOut/AnsysFile)
  SUBDIRS(source/DataInOut/SimInOut/internalMesh)
ENDIF(USE_MESH)

IF(USE_GMV)
  SET(USE_GMV_INPUT ON 
    CACHE INTERNAL "Use GMV files as input for simulation. Attention! Debugging may not work!")

  SET(USE_GMV_OUTPUT ON 
    CACHE INTERNAL "Use GMV files as output for simulation.")

  SUBDIRS(source/DataInOut/SimInOut/gmv)
ENDIF(USE_GMV)

IF(USE_GMSH)
  SUBDIRS(source/DataInOut/SimInOut/gmsh)
ENDIF(USE_GMSH)

IF(USE_HDF5)
  SUBDIRS(source/DataInOut/SimInOut/hdf5)
  SUBDIRS(source/DataInOut/SimInOut/xdmf)
ENDIF(USE_HDF5)

IF(USE_GIDPOST)
  SUBDIRS(source/DataInOut/SimInOut/GiD)
ENDIF(USE_GIDPOST)

IF(USE_UNV)
  SUBDIRS(source/DataInOut/SimInOut/Unverg)
ENDIF(USE_UNV)

IF(USE_ANSYSRST)
  SUBDIRS(source/DataInOut/SimInOut/AnsysRST)
ENDIF(USE_ANSYSRST)

SUBDIRS(source/DataInOut/SimInOut/TextOutput)
SUBDIRS(source/DataInOut/SimInOut/InfoResultOutput)
SUBDIRS(source/DataInOut/SimInOut/Streaming)
SUBDIRS(source/Domain)
SUBDIRS(source/Domain/GridCFS)
SUBDIRS(source/Driver)
SUBDIRS(source/Elements)
SUBDIRS(source/Forms)
SUBDIRS(source/General)
SUBDIRS(source/Materials)
SUBDIRS(source/MatVec)
SUBDIRS(source/ODEDescr)
SUBDIRS(source/ODESolve)
#SUBDIRS(source/ParamIdent)
SUBDIRS(source/PDE)
SUBDIRS(source/Utils)
SUBDIRS(source/Utils/mathParser)
#SUBDIRS(source/Optimization/Optimizer)
#SUBDIRS(source/Optimization/Design)
#SUBDIRS(source/Optimization)
SUBDIRS(source/main)

IF(USE_MPCCI)
  SUBDIRS(source/MpCCIcpl)
ENDIF(USE_MPCCI)

IF(CFSTOOL)
  SUBDIRS(source/cfstool)
ENDIF(CFSTOOL)

IF(CPLREADER)
  SUBDIRS(source/cplreader)
ENDIF(CPLREADER)

#-----------------------------------------------------------------------------
# Build the actual headers with defines from templates.
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/configure_def_headers.cmake")

#-----------------------------------------------------------------------------
# Add targets for building documentation (see 'make help')
#-----------------------------------------------------------------------------
SUBDIRS(share/doc/developer)
SUBDIRS(share/doc/user)

#-----------------------------------------------------------------------------
# Add targets for scripts
#-----------------------------------------------------------------------------
SUBDIRS(share/scripts)

#-----------------------------------------------------------------------------
# Add targets for XML schemata and XSL style sheets
#-----------------------------------------------------------------------------
SUBDIRS(share/xml)
SUBDIRS(share/xsl)
