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
# Include external project build capability of CMake 2.8
#-----------------------------------------------------------------------------
INCLUDE(ExternalProject)

SET(LSE17_SOURCES_DIR "ftp://lse17.e-technik.uni-erlangen.de:40065/cfsdeps/sources")

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CFLAGS "-arch ${CMAKE_OSX_ARCHITECTURES}")
  SET(CFLAGS "${CFLAGS} -sysroot=${CMAKE_OSX_SYSROOT}")
  SET(CFLAGS "${CFLAGS} -isysroot ${CMAKE_OSX_SYSROOT}")
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

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

CONFIGURE_FILE("${CFS_DEPS_ROOT}/build_vars.pl.in"
  "${CFS_TEMP_DIR}/build_vars.pl"
  @ONLY )

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
SET(MUPARSER_URL "${LSE17_SOURCES_DIR}/muparser")
SET(MUPARSER_ZIP "muparser_v2_2_2.zip")
SET(MUPARSER_MD5 "6d77b5cb8096fe2c50afe36ad41bc14a")

INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/muparser/External_muParser.cmake")

#-------------------------------------------------------------------------------
# Build zlib library
#-------------------------------------------------------------------------------
SET(ZLIB_URL "${LSE17_SOURCES_DIR}/zlib")
SET(ZLIB_GZ "zlib-1.2.6.tar.gz")
SET(ZLIB_MD5 "618e944d7c7cd6521551e30b32322f4a")

INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/zlib/External_zlib.cmake")

#-------------------------------------------------------------------------------
# Build bzip2 library
#-------------------------------------------------------------------------------
SET(BZIP2_URL "${LSE17_SOURCES_DIR}/bzip2")
SET(BZIP2_GZ "bzip2-1.0.6.tar.gz")
SET(BZIP2_MD5 "00b516f4704d4a7cb50a1d97e6e8e15b")

INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/bzip2/External_bzip2.cmake")

#-------------------------------------------------------------------------------
# Search for CMake 2.8 if older CMake is used
#-------------------------------------------------------------------------------
INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindCMake.cmake")

#-------------------------------------------------------------------------------
# Search for HDF5 library
#-------------------------------------------------------------------------------
IF(USE_HDF5)
  SET(HDF5_URL "${LSE17_SOURCES_DIR}/hdf5")
  SET(HDF5_GZ "hdf5-1.8.8.tar.gz")
  SET(HDF5_MD5 "1196e668f5592bfb50d1de162eb16cff")

  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/hdf5/External_HDF5.cmake")
ENDIF(USE_HDF5)

#-------------------------------------------------------------------------------
# Search for METIS library
#-------------------------------------------------------------------------------
IF(USE_METIS)
  SET(METIS_URL "${LSE17_SOURCES_DIR}/metis")
  SET(METIS_GZ "metis-4.0.3.tar.gz")
  SET(METIS_MD5 "d3848b454532ef18dc83e4fb160d1e10")

  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/metis/External_METIS.cmake")
ENDIF(USE_METIS)

#-------------------------------------------------------------------------------
# Search for GiDpost library
#-------------------------------------------------------------------------------
IF(USE_GIDPOST)
  SET(GIDPOST_URL "${LSE17_SOURCES_DIR}/gidpost")
  SET(GIDPOST_ZIP "gidpost1.71.zip")
  SET(GIDPOST_MD5 "df8c3ed913cb8abafa36a47591438538")

  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/gidpost/External_GiDpost.cmake")
ENDIF(USE_GIDPOST)

IF(USE_BLAS OR USE_LAPACK)

  #-----------------------------------------------------------------------------
  # Find Netlib BLAS/LAPACK library
  #-----------------------------------------------------------------------------
  IF(CFS_BLAS_LAPACK STREQUAL "GOTO" OR
      CFS_BLAS_LAPACK STREQUAL "NETLIB" OR
      CFS_BLAS_LAPACK STREQUAL "GOTO" OR
      USE_ILUPACK)
    
    SET(LAPACK_URL "${LSE17_SOURCES_DIR}/lapack")
    SET(LAPACK_GZ "lapack-3.2.1.tgz")
    SET(LAPACK_MD5 "a3202a4f9e2f15ffd05d15dab4ac7857")
    
    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/lapack/External_LAPACK.cmake")
    
  ENDIF(CFS_BLAS_LAPACK STREQUAL "GOTO" OR
    CFS_BLAS_LAPACK STREQUAL "NETLIB" OR
    CFS_BLAS_LAPACK STREQUAL "GOTO" OR
    USE_ILUPACK)


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
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindGotoBLAS.cmake")
  ENDIF(CFS_BLAS_LAPACK STREQUAL "GOTO")

  #-----------------------------------------------------------------------------
  # Find Intel Math Kernel library
  #-----------------------------------------------------------------------------
  IF(CFS_BLAS_LAPACK STREQUAL "MKL")
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindIntelMKL.cmake")
  ENDIF(CFS_BLAS_LAPACK STREQUAL "MKL")

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
    SET(ARPACK_URL "${LSE17_SOURCES_DIR}/arpack")
    SET(ARPACK_GZ "arpack96.tar.gz")
    SET(ARPACK_MD5 "fffaa970198b285676f4156cebc8626e")
  
    SET(ARPACK_PATCH_GZ "patch.tar.gz")
    SET(ARPACK_PATCH_MD5 "14830d758f195f272b8594a493501fa2")
  
    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/arpack/External_ARPACK.cmake")
  ENDIF(USE_ARPACK)
  
  #-----------------------------------------------------------------------------
  # Find SuiteSparse/CholMod/AMD library
  #-----------------------------------------------------------------------------
  IF(USE_CHOLMOD OR USE_ILUPACK)
    SET(SUITESPARSE_URL "${LSE17_SOURCES_DIR}/suitesparse")
    SET(SUITESPARSE_GZ "SuiteSparse-3.7.0.tar.gz")
    SET(SUITESPARSE_MD5 "ecb1d1cc1101cf31f077bab46678e791")

    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/suitesparse/External_SuiteSparse.cmake")
  ENDIF(USE_CHOLMOD OR USE_ILUPACK)

  #-----------------------------------------------------------------------------
  # Find ILUPACK library
  #-----------------------------------------------------------------------------
  IF(USE_ILUPACK)
    SET(ILUPACK_PATH "${CFS_SOURCE_DIR}/cfsdeps/ilupack")
    SET(ILUPACK_GZ "ilupack2.2.1_src.tgz")
    SET(ILUPACK_MD5 "83454bbbbb12bd4efca73df50d2e6d7d")
    
    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/ilupack/External_ILUPACK.cmake")
  ENDIF(USE_ILUPACK)

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
SET(BOOST_URL "${LSE17_SOURCES_DIR}/boost")
SET(BOOST_GZ "boost-1.48.0.tar.gz")
SET(BOOST_MD5 "e8614c0ceecce2a388bce03fcc4d73b4")
INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/boost/External_Boost.cmake")

#-------------------------------------------------------------------------------
# Our cfs-hdf5 I/O library depends on Boost, that is why we check for it here.
#-------------------------------------------------------------------------------
# IF(USE_HDF5)
#   INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindCfs-hdf5.cmake")
# ENDIF(USE_HDF5)

#-------------------------------------------------------------------------------
# If USE_XERCES option is defined find Xerces library
#-------------------------------------------------------------------------------
IF(USE_XERCES)
  SET(XERCES_URL "${LSE17_SOURCES_DIR}/xerces")
  SET(XERCES_GZ "xerces-c-3.1.1.tar.gz")
  SET(XERCES_MD5 "6a8ec45d83c8cfb1584c5a5345cb51ae")
  
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/xerces/External_Xerces-C.cmake")
ENDIF(USE_XERCES)

#-----------------------------------------------------------------------------
# Find CGAL
#-----------------------------------------------------------------------------
IF(USE_INTERPOLATION)
  SET(GMP_URL "${LSE17_SOURCES_DIR}/gmp")
  SET(GMP_BZ2 "gmp-4.2.4.tar.bz2")
  SET(GMP_MD5 "fc1e3b3a2a5038d4d74138d0b9cf8dbe")
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/gmp/External_gmp.cmake")
  
  SET(MPFR_URL "${LSE17_SOURCES_DIR}/mpfr")
  SET(MPFR_BZ2 "mpfr-2.4.0.tar.bz2")
  SET(MPFR_MD5 "f5916d785d4f7e7282057f6a3ebff9ce")
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/mpfr/External_mpfr.cmake")

  SET(CGAL_URL "${LSE17_SOURCES_DIR}/cgal")
  SET(CGAL_GZ "CGAL-3.9.tar.gz")
  SET(CGAL_MD5 "797697130ff9231627521c0a38f16d2f")
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/cgal/External_CGAL.cmake")
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
  SET(SCPIP_PATH "${CFS_SOURCE_DIR}/cfsdeps/scpip")
  SET(SCPIP_BZ2 "scpip.tar.bz2")
  SET(SCPIP_MD5 "8afaf8d8d79981d68b8c726ea508471d")

  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/scpip/External_SCPIP.cmake")
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
  IF(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR WIN32)
    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/ansys_custom/External_ANSYS_custom.cmake")
  ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR WIN32)
ENDIF(USE_ANSYSRST)

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

IF(USE_PYTHON)
  set(PYTHON_MAJOR 2)
  set(PYTHON_MINOR 7)
  set(PYTHON_PATCH 2)
  set(PYTHON_VERSION ${PYTHON_MAJOR}.${PYTHON_MINOR}.${PYTHON_PATCH})
  set(PYTHON_URL "http://paraview.org/files/v3.14/dependencies")
  set(PYTHON_GZ Python-${PYTHON_VERSION}.tgz)
  set(PYTHON_MD5 0ddfe265f1b3d0a8c2459f5bf66894c7)
  
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/python/External_Python.cmake")
ENDIF(USE_PYTHON)

SET(CFSDEPS_FOUND 1)
