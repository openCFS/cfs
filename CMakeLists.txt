#-------------------------------------------------------------------------------
# Don't build anything unless the version of CMake is high enough.
#-------------------------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)

#-------------------------------------------------------------------------------
# The name of our project is "CFS".  CMakeLists files in this
# project can refer to the root source directory of the project as
# ${CFS_SOURCE_DIR} and to the root binary directory of the project
# as ${CFS_BINARY_DIR}. We use C, C++ and Fortran as languages.
#-------------------------------------------------------------------------------
PROJECT(CFS C CXX Fortran)

#-----------------------------------------------------------------------------
# CFS version number.  An even minor number corresponds to releases.
#-----------------------------------------------------------------------------
SET(CFS_MAJOR_VERSION 0)
SET(CFS_MINOR_VERSION 6)
SET(CFS_BUILD_VERSION 0)
SET(CFS_VERSION
  "${CFS_MAJOR_VERSION}.${CFS_MINOR_VERSION}.${CFS_BUILD_VERSION}")

#-----------------------------------------------------------------------------
# Set hostname of build host
#-----------------------------------------------------------------------------
SITE_NAME(CFS_BUILD_HOST)
SET(CFS_BUILD_HOST "${CFS_BUILD_HOST}"
  CACHE INTERNAL "CFS_BUILD_HOST")

#-----------------------------------------------------------------------------
# Switch on/off external libraries and other options for compile
#-----------------------------------------------------------------------------
OPTION(USE_LAPACK "Use the LAPACK linear algebra library." ON)
OPTION(USE_ARPACK "Use the ARPACK eigensolver library." ON)
OPTION(USE_PARDISO "Use Pardiso" ON)
OPTION(USE_METIS "Use Metis Graph partitioning." ON)
OPTION(USE_BLAS "Use the basic linear algebra systems library." ON)
OPTION(USE_XERCES "Use the Xerces-C XML parser library." ON)
OPTION(USE_TCL "Use Tool Command Language (TCL) for scripting." ON)
OPTION(USE_PYTHON "Use Python for scripting." OFF)
OPTION(USE_MESH "Use the mesh file reader for input from ANSYS/GiD." ON)
OPTION(USE_GIDPOST "Use the GiD library for simulation output." ON)
OPTION(USE_HDF5 "Use HDF5 for simulation in/output." ON)
OPTION(USE_GMV "Use General Mesh Viewer (GMV) for simulation in/output." ON)
OPTION(USE_UNV "Use Universal File for simulation output." ON)
OPTION(USE_ANSYSRST "Use ANSYS .RST format for simulation output." OFF)
OPTION(USE_TECPLOT "Use Tecplot format for simulation output (mesh only)." OFF)
OPTION(USE_MPCCI "Enable support for the mesh parallel coupling interface."
  OFF)
OPTION(USE_IPOPT "Use the external IPOPT optimizer." OFF)
OPTION(USE_SCPIP "Use the external SCPIP optimizer." ON)

#-----------------------------------------------------------------------------
# Only switch on ILUPACK at lse per default.
#-----------------------------------------------------------------------------
IF(CFS_BUILD_HOST MATCHES "lse")
  OPTION(USE_ILUPACK "Use Preconditioners from ILUPACK." ON)
ELSE(CFS_BUILD_HOST MATCHES "lse")
  OPTION(USE_ILUPACK "Use Preconditioners from ILUPACK." OFF)
ENDIF(CFS_BUILD_HOST MATCHES "lse")

#-----------------------------------------------------------------------------
# Use interpolation toggles the inclusion of the CGAL library.
#-----------------------------------------------------------------------------
OPTION(USE_INTERPOLATION "Activate Interpolation between grids." OFF)

#-----------------------------------------------------------------------------
# USE_ADPATGRID is only queried in the testsuite for the moment.
#-----------------------------------------------------------------------------
OPTION(USE_ADAPTGRID "Use h-Adaptivity for meshes." OFF)

OPTION(USE_OPENMP "Turn on support for OpenMP. Needs GCC >= 4.2 or Intel C++." OFF)

OPTION(DEBUG "Build debug version of the code." ON)
OPTION(PROFILING "Turn on profiling (DEBUG/RELEASE)." OFF)

# Set a sane default value for the CFS++ schema path
SET(XMLSCHEMA "${CFS_SOURCE_DIR}/XML" CACHE STRING
  "Path where the XML schema for the input files is located.")

OPTION(CONVERTER "Build CFS++ converter tool (cfstool)." OFF)
OPTION(CPLREADER "Build cplreader for fluid/acoustic coupling." OFF)

#-----------------------------------------------------------------------------
# Specify path to CFS++ testsuite. If this path exists, CTest compatible
# test cases will be generated in the current binary directory.
#-----------------------------------------------------------------------------
SET(TESTSUITE_DIR "TESTSUITE_DIR_NOT_FOUND" CACHE STRING
  "Path where the CFS++ testsuite is located.")

IF(EXISTS ${TESTSUITE_DIR})
  # In order to perform tests we need the CONVERTER and HDF5!
  SET(CONVERTER ON CACHE BOOL
    "Build CFS++ converter tool (cfstool)." FORCE)

  SET(USE_HDF5 ON CACHE BOOL
    "Use HDF5 for simulation in/output." FORCE)

  SUBDIRS(${TESTSUITE_DIR})
ENDIF(EXISTS ${TESTSUITE_DIR})

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

#-----------------------------------------------------------------------------
# Set some variables according to build type DEBUG / RELEASE
#-----------------------------------------------------------------------------
IF(DEBUG)
  SET(CMAKE_BUILD_TYPE "DEBUG" CACHE STRING
    "Only DEBUG and RELEASE types are valid." FORCE)
  SET(_DEBUG 1)
  
  #---------------------------------------------------------------------------
  # Defines related to tracing of functions and memory usage
  #---------------------------------------------------------------------------
  SET(MEMTRACE 1)
  
  #---------------------------------------------------------------------------
  # Check for certain error conditions (in DEBUG mode)
  #---------------------------------------------------------------------------
  SET(CHECK_INDEX 1)
  SET(CHECK_INITIALIZED 1)
  SET(CHECK_MEM_ALLOC 1)
  SET(CHECK_MEMORY 1)
  SET(CHECK_TYPE_CASTS 1)

ELSE(DEBUG)
  SET(CMAKE_BUILD_TYPE "RELEASE" CACHE STRING
    "Only DEBUG and RELEASE types are valid." FORCE)
  SET(ENABLE_PROFILING OFF CACHE BOOL
    "Turn on profiling (DEBUG only)."FORCE)

  #---------------------------------------------------------------------------
  # Activate expression templates (in RELEASE mode)
  #---------------------------------------------------------------------------
  SET(EXPR_TEMPLATES 1)
ENDIF(DEBUG)

    #---------------------------------------------------------------------------
    # Defines related to tracing of functions and memory usage
    #---------------------------------------------------------------------------
    SET(MEMTRACE 1)

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
# Find executables of a few required programs
#-----------------------------------------------------------------------------
INCLUDE("cmake_modules/FindPrograms.cmake")


#-----------------------------------------------------------------------------
# Set path where shared/static libraries should go
#-----------------------------------------------------------------------------
IF(CFS_ARCH STREQUAL "X86_64")
  SET(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib64 CACHE INTERNAL
    "Single output directory for building all libraries.")
ELSE(CFS_ARCH STREQUAL "X86_64")
  SET(LIBRARY_OUTPUT_PATH ${CFS_BINARY_DIR}/lib CACHE INTERNAL
    "Single output directory for building all libraries.")
ENDIF(CFS_ARCH STREQUAL "X86_64")

#-----------------------------------------------------------------------------
# Set path where executables should go
#-----------------------------------------------------------------------------
SET(EXECUTABLE_OUTPUT_PATH ${CFS_BINARY_DIR}/bin CACHE INTERNAL
  "Single output directory for building all executables.")


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
FILE(READ "${CFS_SOURCE_DIR}/.svn/entries" SVNSTR)
STRING(REGEX MATCH "[a-z+]+://[a-zA-Z0-9_+-/]+" CFS_SUBVERSION_REPOS "${SVNSTR}")

#-----------------------------------------------------------------------------
# Determine global revision of working copy
#-----------------------------------------------------------------------------
IF(SVNVERSION)
EXEC_PROGRAM(${SVNVERSION}
  ARGS ${CFS_SOURCE_DIR}
  OUTPUT_VARIABLE CFS_SUBVERSION_REV
  RETURN_VALUE RETVAL)

IF(NOT RETVAL STREQUAL 0)
  SET(CFS_SUBVERSION_REV "N/A")
ENDIF(NOT RETVAL STREQUAL 0)
ELSE(SVNVERSION)
  SET(CFS_SUBVERSION_REV "N/A")
ENDIF(SVNVERSION)


#-----------------------------------------------------------------------------
# Make sure the compiler can find include files from CFS.
#-----------------------------------------------------------------------------
INCLUDE_DIRECTORIES(
  ${Boost_INCLUDE_DIR} 
  ${CFS_BINARY_DIR}/include
  ${CFS_SOURCE_DIR}/source
  ${CFS_SOURCE_DIR}/source/OLAS)

IF(USE_INTERPOLATION)
  INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIR})
ENDIF(USE_INTERPOLATION)

#-----------------------------------------------------------------------------
# Make sure the linker can find the CFS base libraries once it is built.
#-----------------------------------------------------------------------------
LINK_DIRECTORIES(${CFSDEPS_LIBRARY_DIR} ${CFS_BINARY_DIR}/lib)



#-----------------------------------------------------------------------------
# Recurse into the "source" subdirectory.  This does not actually
# cause another cmake executable to run.  The same process will
# walk through the project's entire directory structure.
#-----------------------------------------------------------------------------
SUBDIRS(source/OLAS/algsys)
SUBDIRS(source/OLAS/graph)
SUBDIRS(source/OLAS/matvec)
#  SUBDIRS(source/multigrid)
SUBDIRS(source/OLAS/precond)
SUBDIRS(source/OLAS/solver)
SUBDIRS(source/OLAS/utils)
#  SUBDIRS(source/main)

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

SUBDIRS(source/CoupledPDE)
SUBDIRS(source/DataInOut)
SUBDIRS(source/DataInOut/CommandLine)
SUBDIRS(source/DataInOut/DefineFiles)
SUBDIRS(source/DataInOut/Logging)
SUBDIRS(source/DataInOut/ParamHandling)

IF(USE_SCRIPTING)
  SUBDIRS(source/DataInOut/Scripting)

  IF(USE_SCRIPTING_PYTHON)
    SUBDIRS(source/DataInOut/Scripting/python)
  ENDIF(USE_SCRIPTING_PYTHON)

  IF(USE_SCRIPTING_TCL)
    SUBDIRS(source/DataInOut/Scripting/tcl)
  ENDIF(USE_SCRIPTING_TCL)
ENDIF(USE_SCRIPTING)


IF(USE_MESH)
  SUBDIRS(source/DataInOut/SimInOut/AnsysFile)
ENDIF(USE_MESH)

IF(USE_GMV)
  OPTION(USE_GMV_INPUT
    "Use GMV files as input for simulation. Attention! Debugging may not work!"
    OFF)

  OPTION(USE_GMV_OUTPUT
    "Use GMV files as output for simulation."
    ON)

  MARK_AS_ADVANCED(USE_GMV_INPUT USE_GMV_OUTPUT)
  SUBDIRS(source/DataInOut/SimInOut/gmv)
ENDIF(USE_GMV)

IF(USE_HDF5)
  SUBDIRS(source/DataInOut/SimInOut/hdf5)
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

IF(USE_TECPLOT)
  SUBDIRS(source/DataInOut/SimInOut/tecplot)
ENDIF(USE_TECPLOT)

SUBDIRS(source/DataInOut/SimInOut/TextOutput)
SUBDIRS(source/Domain)
SUBDIRS(source/Domain/GridCFS)
SUBDIRS(source/Driver)
SUBDIRS(source/Elements)
#  SUBDIRS(source/Estimator)
SUBDIRS(source/Forms)
SUBDIRS(source/General)
SUBDIRS(source/Materials)
SUBDIRS(source/Matrix)
SUBDIRS(source/ODEDescr)
SUBDIRS(source/ODESolve)
SUBDIRS(source/ParamIdent)
SUBDIRS(source/PDE)
SUBDIRS(source/Utils)
SUBDIRS(source/Utils/mathParser)
SUBDIRS(source/Optimization)
SUBDIRS(source/main)

IF(USE_MPCCI)
  SUBDIRS(source/MpCCIcpl)
ENDIF(USE_MPCCI)

IF(CONVERTER)
  SUBDIRS(source/cfstool)
ENDIF(CONVERTER)

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
SUBDIRS(Doc/developer)
SUBDIRS(Doc/user)

#-----------------------------------------------------------------------------
# Add targets for scripts
#-----------------------------------------------------------------------------
SUBDIRS(auxiliary/scripts)

