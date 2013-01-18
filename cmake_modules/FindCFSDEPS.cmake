# - Find CFSDEPS installation 
# This module finds if CFSDEPS is installed and determines where 
# the include files and libraries are. It also determines what the name of
# the library is. This code sets the following variables:
#  
#  CFSDEPS_FOUND     = system has CFSDEPS 
#  CFSDEPS_LIBRARIES = path to the CFSDEPS libraries
#  CFSDEPS_CXX_FLAGS  = Compiler flags for CFSDEPS 
#  CFSDEPS_INCLUDE_DIR      = where to find "CFSDEPS.h"
#  CFSDEPS_DEFINITIONS      = extra defines
#
# AUTHOR
# Simon Triebenbacher simon.triebenbacher@uni-klu.ac.at (02/2009)

#-----------------------------------------------------------------------------
# If user has set environment variables use them. If not use defaults
#-----------------------------------------------------------------------------
SET(CFS_DEPS_ROOT_DUMMY "$ENV{CFS_DEPS_ROOT}")
IF(NOT ${CFS_DEPS_ROOT_DUMMY} STREQUAL "")
  SET(CFS_DEPS_ROOT "${CFS_DEPS_ROOT_DUMMY}" CACHE PATH
    "Root directory of CFSDEPS (environment).")
ELSE(NOT ${CFS_DEPS_ROOT_DUMMY} STREQUAL "")
  SET(CFS_DEPS_ROOT "${CFS_DEPS_ROOT_DEFAULT}" CACHE PATH
    "Root directory of CFSDEPS (default).")
ENDIF(NOT ${CFS_DEPS_ROOT_DUMMY} STREQUAL "")

SET(CFS_DEPS_CD_DUMMY "$ENV{CFS_DEPS_CACHE_DIR}")
IF(NOT ${CFS_DEPS_CD_DUMMY} STREQUAL "")
  SET(CFS_DEPS_CACHE_DIR "${CFS_DEPS_CD_DUMMY}" CACHE PATH
    "Directory for CFSDEPS sources and prebuilt binaries.")
ELSE(NOT ${CFS_DEPS_CD_DUMMY} STREQUAL "")
  SET(CFS_DEPS_CACHE_DIR "${CFS_DEPS_CACHE_DIR_DEFAULT}" CACHE PATH
    "Directory for CFSDEPS sources and prebuilt binaries.")
ENDIF(NOT ${CFS_DEPS_CD_DUMMY} STREQUAL "")

SET(CFS_FORCE_DEPS_DUMMY "$ENV{CFS_FORCE_DEPS_CACHE_DIR}")
IF(NOT ${CFS_FORCE_DEPS_DUMMY} STREQUAL "")
  SET(CFS_FORCE_DEPS_CACHE_DIR ON CACHE PATH
    "Force 'CFS_DEPS_CACHE_DIR/precompiled/forced'.")
ENDIF(NOT ${CFS_FORCE_DEPS_DUMMY} STREQUAL "")
#-----------------------------------------------------------------------------
# Check if the proper files are present in the CFSDEPS directory
#-----------------------------------------------------------------------------
IF(NOT EXISTS "${CFS_DEPS_ROOT}/build_common.pl")
  MESSAGE(FATAL_ERROR "You obviously do not have 'build_common.pl' in "
    "'${CFS_DEPS_ROOT}'. Get CFSDEPS from svn+ssh://lse10/software/cfsdeps/trunk "
    "and place it in /opt/CFSDEPS. You can also just copy the contents of the "
    "directory /home/data/libraries/CFSDEPS. If you want to use different "
    "directories for CFS_DEPS_ROOT and CFS_DEPS_CACHE_DIR, please set "
    "environment variables with the same name before starting CMake or use "
    "its -D switch.")
ENDIF(NOT EXISTS "${CFS_DEPS_ROOT}/build_common.pl")

#-----------------------------------------------------------------------------
# Check if cache directory is present
#-----------------------------------------------------------------------------
IF(NOT CFS_DEPS_CACHE_DIR)
  MESSAGE(FATAL_ERROR "Please set CFS_DEPS_CACHE_DIR. "
    "This dirctory is used to store downloaded sources and prebuilt "
    "binaries, which can be reused for other CFS++ builds."
    "This directory may even be located on a network share.")
ENDIF(NOT CFS_DEPS_CACHE_DIR)

FILE(TO_CMAKE_PATH
  "${CFS_DEPS_CACHE_DIR}"
  CFS_DEPS_CACHE_DIR)

CONFIGURE_FILE("${CFS_DEPS_ROOT}/build_vars.pl.in"
  "${CFS_TEMP_DIR}/build_vars.pl"
  @ONLY )


#-------------------------------------------------------------------------------
# Build MuParser library
#-------------------------------------------------------------------------------
INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindMuParser.cmake")

#-------------------------------------------------------------------------------
# Build zlib library
#-------------------------------------------------------------------------------
INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindZlib.cmake")

#-------------------------------------------------------------------------------
# Build bzip2 library
#-------------------------------------------------------------------------------
INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindBzip2.cmake")

#-------------------------------------------------------------------------------
# Search for CMake 2.8 if older CMake is used
#-------------------------------------------------------------------------------
INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindCMake.cmake")

#-------------------------------------------------------------------------------
# Search for HDF5 library
#-------------------------------------------------------------------------------
IF(USE_HDF5)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindHDF5.cmake")
ENDIF(USE_HDF5)

#-------------------------------------------------------------------------------
# Search for METIS library
#-------------------------------------------------------------------------------
IF(USE_METIS)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindMetis.cmake")
ENDIF(USE_METIS)

#-------------------------------------------------------------------------------
# Search for GiDpost library
#-------------------------------------------------------------------------------
IF(USE_GIDPOST)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindGiDpost.cmake")
ENDIF(USE_GIDPOST)

#-----------------------------------------------------------------------------
# Find FFTW for cfstool fftw feature
#-----------------------------------------------------------------------------
IF(CFSTOOL_FFTW)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindFFTW.cmake")
ENDIF(CFSTOOL_FFTW)

IF(USE_BLAS OR USE_LAPACK)
  #-----------------------------------------------------------------------------
  # Find ACML library
  #-----------------------------------------------------------------------------
  IF(CFS_BLAS_LAPACK STREQUAL "ACML")
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindACML.cmake")
  ENDIF(CFS_BLAS_LAPACK STREQUAL "ACML")

  #-----------------------------------------------------------------------------
  # Find GotoBLAS library
  #-----------------------------------------------------------------------------
  IF(CFS_BLAS_LAPACK STREQUAL "GOTO")
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindLAPACK.cmake")
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindGotoBLAS.cmake")
  ENDIF(CFS_BLAS_LAPACK STREQUAL "GOTO")

  #-----------------------------------------------------------------------------
  # Find Intel Math Kernel library
  #-----------------------------------------------------------------------------
  IF(CFS_BLAS_LAPACK STREQUAL "MKL")
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindIntelMKL.cmake")
  ENDIF(CFS_BLAS_LAPACK STREQUAL "MKL")

  #-----------------------------------------------------------------------------
  # Find Netlib BLAS/LAPACK library
  #-----------------------------------------------------------------------------
  IF(CFS_BLAS_LAPACK STREQUAL "NETLIB")
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindLAPACK.cmake")
  ENDIF(CFS_BLAS_LAPACK STREQUAL "NETLIB")

  #-----------------------------------------------------------------------------
  # Search for Pardiso Library
  #-----------------------------------------------------------------------------
  IF(USE_PARDISO AND CFS_PARDISO STREQUAL "SCHENK")
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindPardiso.cmake")
  ENDIF(USE_PARDISO AND CFS_PARDISO STREQUAL "SCHENK")

  #-----------------------------------------------------------------------------
  # Check which version of the Pardiso API is being used. Pardiso 4.0 intro-
  # duces a new incompatible API due to the new iterative solver feature.
  # However MKL still uses the old API of Pardiso 3.x. This is obviously a
  # software engineering nightmare and the ones responsible for it deserve to
  # be hanged from the highest mast!
  #
  # ATTENTION: If you switch CFS_PARDISO you should also clear
  # PARDISO_API_VER_3 and PARDISO_API_VER_4 from the CMake cache.
  #-----------------------------------------------------------------------------
  IF(USE_PARDISO)
    INCLUDE("cmake_modules/CheckPardisoAPIVersion.cmake")
  ENDIF(USE_PARDISO)
  
  #-----------------------------------------------------------------------------
  # If USE_ARPACK option is defined find ARPACK library
  #-----------------------------------------------------------------------------
  IF(USE_ARPACK)
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindARPACK.cmake")
  ENDIF(USE_ARPACK)
  
  #-----------------------------------------------------------------------------
  # Find ILUPACK library
  #-----------------------------------------------------------------------------
  IF(USE_ILUPACK)
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindLAPACK.cmake")
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindSuiteSparse.cmake")    
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindILUPACK.cmake")
  ENDIF(USE_ILUPACK)

  #-----------------------------------------------------------------------------
  # Find SuiteSparse/CholMod/AMD library
  #-----------------------------------------------------------------------------
  IF(USE_CHOLMOD)
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindSuiteSparse.cmake")
  ENDIF(USE_CHOLMOD)

  #-----------------------------------------------------------------------------
  # Find SnOpt library
  #-----------------------------------------------------------------------------
  IF(USE_SNOPT)
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindSnOpt.cmake")
  ENDIF(USE_SNOPT)

#  MESSAGE("BLAS_LIBRARY ${BLAS_LIBRARY}")
#  MESSAGE("LAPACK_LIBRARY ${LAPACK_LIBRARY}")
#  MESSAGE("PARDISO_LIBRARY ${PARDISO_LIBRARY}")
  
ENDIF(USE_BLAS OR USE_LAPACK)


#-------------------------------------------------------------------------------
# Find Boost
#-------------------------------------------------------------------------------
INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindBoostForCFS.cmake")

#-------------------------------------------------------------------------------
# Our cfs-hdf5 I/O library depends on Boost, that is why we check for it here.
#-------------------------------------------------------------------------------
IF(USE_HDF5)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindCfs-hdf5.cmake")
ENDIF(USE_HDF5)

#-------------------------------------------------------------------------------
# If USE_PYTHON option is defined find Python library
#-------------------------------------------------------------------------------
IF(USE_PYTHON)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindPythonForCFS.cmake")
ENDIF(USE_PYTHON)

#-------------------------------------------------------------------------------
# If USE_TCL option is defined find TCL library
#-------------------------------------------------------------------------------
IF(USE_TCL)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindTCLForCFS.cmake")
ENDIF(USE_TCL)

#-------------------------------------------------------------------------------
# If USE_XERCES option is defined find Xerces library
#-------------------------------------------------------------------------------
IF(USE_XERCES)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindXerces.cmake")
ENDIF(USE_XERCES)

#-----------------------------------------------------------------------------
# Find MpCCI / MPICH
#-----------------------------------------------------------------------------
IF(USE_MPCCI)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindMpCCI.cmake")    
ENDIF(USE_MPCCI)

#-----------------------------------------------------------------------------
# Find CGAL
#-----------------------------------------------------------------------------
IF(USE_INTERPOLATION)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindCGAL.cmake")    
ENDIF(USE_INTERPOLATION)

#-----------------------------------------------------------------------------
# Find IPOPT
#-----------------------------------------------------------------------------
IF(USE_IPOPT)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindIPOPT.cmake")    
ENDIF(USE_IPOPT)

#-----------------------------------------------------------------------------
# Find SCPIP
#-----------------------------------------------------------------------------
IF(USE_SCPIP)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindSCPIP.cmake")    
ENDIF(USE_SCPIP)

#-----------------------------------------------------------------------------
# Find SnOpt
#-----------------------------------------------------------------------------
IF(USE_SNOPT)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindSnOpt.cmake")    
ENDIF(USE_SNOPT)

#-----------------------------------------------------------------------------
# Find ANSYS Customizations
#-----------------------------------------------------------------------------
IF(USE_ANSYSRST)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindANSYSCust.cmake")
ENDIF(USE_ANSYSRST)

#-----------------------------------------------------------------------------
# Find CFX Customizations
#-----------------------------------------------------------------------------
IF(CPLREADER_CFX)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindCFXCust.cmake")
ENDIF(CPLREADER_CFX)

#-----------------------------------------------------------------------------
# Find VTK for OpenFOAM, EnSight and FLUENT readers
#-----------------------------------------------------------------------------
IF(CPLREADER_OPENFOAM OR CPLREADER_ENSIGHT OR CPLREADER_FLUENT)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindVTK.cmake")
ENDIF(CPLREADER_OPENFOAM OR CPLREADER_ENSIGHT OR CPLREADER_FLUENT)

#-----------------------------------------------------------------------------
# Find CGNS for cgns reader
#-----------------------------------------------------------------------------
IF(CPLREADER_CGNS)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindCGNS.cmake")
ENDIF(CPLREADER_CGNS)

#-----------------------------------------------------------------------------
# Find ParaView postprocessor
#-----------------------------------------------------------------------------
IF(BUILD_PARAVIEW)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindParaView.cmake")
ENDIF(BUILD_PARAVIEW)

#-----------------------------------------------------------------------------
# Find HDF file viewer
#-----------------------------------------------------------------------------
IF(BUILD_HDFVIEW)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindHDFView.cmake")
ENDIF(BUILD_HDFVIEW)


EXECUTE_PROCESS(COMMAND "${PERL}" "${CFS_DEPS_ROOT}/utils/perl/cfsdepsmodif.pl"
  ${CFS_DEPS_ROOT} ${CFS_BUILD_DIR}/tmp
  RESULT_VARIABLE retval
  OUTPUT_VARIABLE CFS_DEPS_MODIFIED)

# the boolean $CFS_FORCE_DEPS_CACHE_DIR is a expert feature for people knowing what they do
IF(CFS_DEPS_MODIFIED MATCHES "WORKING_COPY_MODIFIED modified")
  IF(NOT CFS_FORCE_DEPS_CACHE_DIR)
    MESSAGE("The Subversion working copy of CFSDEPS in '${CFS_DEPS_ROOT}' has been modified. Precompiled binaries have therefore been put into '${CFS_BINARY_DIR}/tmp' instead of '${CFS_DEPS_CACHE_DIR}/precompiled'!")
  ENDIF(NOT CFS_FORCE_DEPS_CACHE_DIR)  
ENDIF(CFS_DEPS_MODIFIED MATCHES "WORKING_COPY_MODIFIED modified")

SET(CFSDEPS_FOUND 1)
