#=============================================================================
#
# Find locations of external binary libs (e.g. MKL) and build additional
# external libs from source.
# 
# This module finds and builds libs that CFS++ depends upon and determines
# where the include files and libraries are. 
#  
# AUTHORS
# - Simon Triebenbacher simon.triebenbacher@uni-klu.ac.at (02/2009)
# - Simon Triebenbacher simon.triebenbacher@tuwien.ac.at (01/2013)
#
#=============================================================================

#-----------------------------------------------------------------------------
# Include external project build capability of CMake 2.8
#-----------------------------------------------------------------------------
INCLUDE(ExternalProject)

SET(LSE17_SOURCES_DIR "ftp://lse17.e-technik.uni-erlangen.de:40065")
SET(LSE17_SOURCES_DIR "${LSE17_SOURCES_DIR}/cfsdeps/sources")

#-----------------------------------------------------------------------------
# Set common CFLAGS (and CXXFLAGS) common for all external projects.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# We do not want to see warnings from external projects, since they would
# show up on CDash.
#-----------------------------------------------------------------------------
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CFSDEPS_C_FLAGS "-w")
  SET(CFSDEPS_CXX_FLAGS "-w")
  SET(CFSDEPS_Fortran_FLAGS "-w")
ELSEIF(MSVC)
  STRING(REPLACE "/W3" "/w" CFSDEPS_C_FLAGS "${CMAKE_C_FLAGS_INIT}")
  STRING(REPLACE "/W3" "/w" CFSDEPS_CXX_FLAGS "${CMAKE_CXX_FLAGS_INIT}")
  STRING(REPLACE "/W3" "/w" CFSDEPS_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_INIT}")
ENDIF()

#-----------------------------------------------------------------------------
# On Mac OS X we want to build the external libs for the same SDK and 
# architecture as CFS++.
#-----------------------------------------------------------------------------
IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -arch ${CMAKE_OSX_ARCHITECTURES}")
  SET(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -sysroot=${CMAKE_OSX_SYSROOT}")
  SET(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -isysroot ${CMAKE_OSX_SYSROOT}")

  SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -arch ${CMAKE_OSX_ARCHITECTURES}")
  SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -sysroot=${CMAKE_OSX_SYSROOT}")
  SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -isysroot ${CMAKE_OSX_SYSROOT}")
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

#-----------------------------------------------------------------------------
# If user has set environment variables use them. If not use defaults
#-----------------------------------------------------------------------------
SET(CFS_DEPS_CD_DUMMY "$ENV{CFS_DEPS_CACHE_DIR}")
IF(NOT ${CFS_DEPS_CD_DUMMY} STREQUAL "")
  SET(CFS_DEPS_CACHE_DIR "${CFS_DEPS_CD_DUMMY}" CACHE PATH
    "Directory for CFSDEPS sources and prebuilt binaries.")
ELSE(NOT ${CFS_DEPS_CD_DUMMY} STREQUAL "")
  SET(CFS_DEPS_CACHE_DIR "${CFS_DEPS_CACHE_DIR_DEFAULT}" CACHE PATH
    "Directory for CFSDEPS sources and prebuilt binaries.")
ENDIF(NOT ${CFS_DEPS_CD_DUMMY} STREQUAL "")

#-----------------------------------------------------------------------------
# Check if cache directory is present
#-----------------------------------------------------------------------------
IF(NOT CFS_DEPS_CACHE_DIR)
  MESSAGE(FATAL_ERROR "Please set CFS_DEPS_CACHE_DIR. 
This dirctory is used to store downloaded sources, 
which can be reused for other CFS++ builds.
This directory may even be located on a network share.")
ENDIF(NOT CFS_DEPS_CACHE_DIR)

FILE(TO_CMAKE_PATH
  "${CFS_DEPS_CACHE_DIR}"
  CFS_DEPS_CACHE_DIR)

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
SET(ZLIB_GZ "zlib-1.2.8.tar.gz")
SET(ZLIB_MD5 "44d667c142d7cda120332623eab69f40")

INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/zlib/External_zlib.cmake")

#-------------------------------------------------------------------------------
# Build bzip2 library
#-------------------------------------------------------------------------------
SET(BZIP2_URL "${LSE17_SOURCES_DIR}/bzip2")
SET(BZIP2_GZ "bzip2-1.0.6.tar.gz")
SET(BZIP2_MD5 "00b516f4704d4a7cb50a1d97e6e8e15b")

INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/bzip2/External_bzip2.cmake")

#-------------------------------------------------------------------------------
# Search for HDF5 library
#-------------------------------------------------------------------------------
IF(USE_HDF5)
  SET(HDF5_URL "${LSE17_SOURCES_DIR}/hdf5")
  SET(HDF5_GZ "hdf5-1.8.12.tar.gz")
  SET(HDF5_MD5 "03ad766d225f5e872eb3e5ce95524a08")

  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/hdf5/External_HDF5.cmake")
ENDIF(USE_HDF5)

#-------------------------------------------------------------------------------
# Search for CGNS library
#-------------------------------------------------------------------------------
IF(USE_CGNS)
  SET(CGNS_URL "${LSE17_SOURCES_DIR}/cgns")
  SET(CGNS_GZ "cgnslib_3.1.4-2.tar.gz")
  SET(CGNS_MD5 "e2d57dc5e116ff723ee003eba667b9f9")

  SET(CGNS25_GZ "cgnslib_2.5-5.tar.gz")
  SET(CGNS25_MD5 "ae2a2e79b99d41c63e5ed5f661f70fd9")

  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/cgns/External_CGNS.cmake")
ENDIF(USE_CGNS)

#-------------------------------------------------------------------------------
# Search for STAR-CCM+ I/O library
#-------------------------------------------------------------------------------
IF(USE_CCMIO)
  SET(CCMIO_URL "${LSE17_SOURCES_DIR}/ccmio")
  SET(CCMIO_GZ "libccmio-2.6.1.tar.gz")
  SET(CCMIO_MD5 "f81fbdfb960b1a4f3bcc7feee491efe4")

  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/ccmio/External_CCMIO.cmake")
ENDIF(USE_CCMIO)

#-------------------------------------------------------------------------------
# Search for CFX I/O library
#-------------------------------------------------------------------------------
IF(USE_CFXIO)
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/cfx_custom/External_CFX_Custom.cmake")
ENDIF(USE_CFXIO)

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
  SET(GIDPOST_ZIP "gidpost-2.1.zip")
  SET(GIDPOST_MD5 "a7fe745e40593dc4598920b312663d29")

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
    SET(LAPACK_GZ "lapack-3.4.2.tgz")
    SET(LAPACK_MD5 "61bf1a8a4469d4bdb7604f5897179478")
    
    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/lapack/External_LAPACK.cmake")
    
  ENDIF(CFS_BLAS_LAPACK STREQUAL "GOTO" OR
    CFS_BLAS_LAPACK STREQUAL "NETLIB" OR
    CFS_BLAS_LAPACK STREQUAL "GOTO" OR
    USE_ILUPACK)


  #-----------------------------------------------------------------------------
  # Find ACML library
  #-----------------------------------------------------------------------------
  IF(CFS_BLAS_LAPACK STREQUAL "ACML")
    MESSAGE(FATAL_ERROR "ACML has not been ported to CMake externals yet.")
#    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindACML.cmake")
  ENDIF(CFS_BLAS_LAPACK STREQUAL "ACML")

  #-----------------------------------------------------------------------------
  # Find GotoBLAS library
  #-----------------------------------------------------------------------------
  IF(CFS_BLAS_LAPACK STREQUAL "GOTO")
    MESSAGE(FATAL_ERROR "GotoBLAS has not been ported to CMake externals yet.")
#    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindGotoBLAS.cmake")
  ENDIF(CFS_BLAS_LAPACK STREQUAL "GOTO")

  #-----------------------------------------------------------------------------
  # Find Intel Math Kernel library
  #-----------------------------------------------------------------------------
  IF(CFS_BLAS_LAPACK STREQUAL "MKL")
    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindIntelMKL.cmake")
  ENDIF(CFS_BLAS_LAPACK STREQUAL "MKL")

  IF(CFS_BLAS_LAPACK STREQUAL "APPLE")
    SET(BLAS_LIBRARY "-framework Accelerate")
    SET(LAPACK_LIBRARY "-framework Accelerate")
  ENDIF(CFS_BLAS_LAPACK STREQUAL "APPLE")

  #-----------------------------------------------------------------------------
  # Search for Pardiso Library
  #-----------------------------------------------------------------------------
  IF(USE_PARDISO AND CFS_PARDISO STREQUAL "SCHENK")
    MESSAGE(FATAL_ERROR "Pardiso has not been ported to CMake externals yet.")
#    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindPardiso.cmake")
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
    SET(ARPACK_GZ "arpack-ng_3.1.1.tar.gz")
    SET(ARPACK_MD5 "d65b915736650d8878719d4168e50c36")
  
    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/arpack/External_ARPACK.cmake")
  ENDIF(USE_ARPACK)
  
  #-----------------------------------------------------------------------------
  # Find SuiteSparse/CholMod/UMFPACK/AMD library
  #-----------------------------------------------------------------------------
  IF(USE_SUITESPARSE OR USE_ILUPACK)
    SET(SUITESPARSE_URL "${LSE17_SOURCES_DIR}/suitesparse")
    SET(SUITESPARSE_GZ "SuiteSparse-4.2.1.tar.gz")
    SET(SUITESPARSE_MD5 "4628df9eeae10ae5f0c486f1ac982fce")

    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/suitesparse/External_SuiteSparse.cmake")
  ENDIF(USE_SUITESPARSE OR USE_ILUPACK)

  #-----------------------------------------------------------------------------
  # Find ILUPACK library
  #-----------------------------------------------------------------------------
  IF(USE_ILUPACK)
    SET(ILUPACK_PATH "${CFS_SOURCE_DIR}/cfsdeps/ilupack")
    SET(ILUPACK_GZ "ilupack2.2.1_src.tgz")
    SET(ILUPACK_MD5 "7cb6ba2e854e13d243218d9e9478d13c")
    
    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/ilupack/External_ILUPACK.cmake")
  ENDIF(USE_ILUPACK)

  #-----------------------------------------------------------------------------
  # Find SnOpt library
  #-----------------------------------------------------------------------------
  IF(USE_SNOPT)
    MESSAGE(FATAL_ERROR "SnOpt has not been ported to CMake externals yet.")
#    INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindSnOpt.cmake")
  ENDIF(USE_SNOPT)

#  MESSAGE("BLAS_LIBRARY ${BLAS_LIBRARY}")
#  MESSAGE("LAPACK_LIBRARY ${LAPACK_LIBRARY}")
#  MESSAGE("PARDISO_LIBRARY ${PARDISO_LIBRARY}")
  
ENDIF(USE_BLAS OR USE_LAPACK)

#-------------------------------------------------------------------------------
# Find Library of Iterative Solvers
#-------------------------------------------------------------------------------
IF(USE_LIS)
  SET(LIS_URL "${LSE17_SOURCES_DIR}/lis")
  SET(LIS_GZ "lis-1.3.16.tar.gz")
  SET(LIS_MD5 "7bbbcd2070cca367a98d17767c0ea408")
  
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/lis/External_LIS.cmake")
ENDIF(USE_LIS)

#-------------------------------------------------------------------------------
# Find SuperLU
#-------------------------------------------------------------------------------
IF(USE_SUPERLU)
  SET(SUPERLU_URL "${LSE17_SOURCES_DIR}/superlu")
  SET(SUPERLU_GZ "superlu_4.3.tar.gz")
  SET(SUPERLU_MD5 "b72c6309f25e9660133007b82621ba7c")
  
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/superlu/External_SuperLU.cmake")
ENDIF(USE_SUPERLU)

#-------------------------------------------------------------------------------
# Find Boost
#-------------------------------------------------------------------------------
SET(BOOST_URL "${LSE17_SOURCES_DIR}/boost")
SET(BOOST_GZ "boost_1_52_0.tar.bz2")
SET(BOOST_MD5 "3a855e0f919107e0ca4de4d84ad3f750")
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
IF(USE_CGAL)
  IF(WIN32)
    IF(MINGW AND NOT CMAKE_CROSSCOMPILING OR MSVC)
      SET(MSG "The build of gmp and mpfr is not supported on Windows!")
      SET(MSG "${MSG} It is configure-based and therefore requires a shell")
      SET(MSG "${MSG} interpreter like bash from MSYS.")
      SET(MSG "${MSG} If you need CGAL, you need to cross compile from Linux.")
      MESSAGE(FATAL_ERROR "${MSG}")
    ENDIF()
  ENDIF()

  SET(GMP_URL "${LSE17_SOURCES_DIR}/gmp")
  SET(GMP_BZ2 "gmp-4.2.4.tar.bz2")
  SET(GMP_MD5 "fc1e3b3a2a5038d4d74138d0b9cf8dbe")
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/gmp/External_gmp.cmake")
  
  SET(MPFR_URL "${LSE17_SOURCES_DIR}/mpfr")
  SET(MPFR_BZ2 "mpfr-2.4.0.tar.bz2")
  SET(MPFR_MD5 "f5916d785d4f7e7282057f6a3ebff9ce")
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/mpfr/External_mpfr.cmake")

  SET(CGAL_URL "${LSE17_SOURCES_DIR}/cgal")
  SET(CGAL_BZ2 "CGAL-4.2.tar.bz2")
  SET(CGAL_MD5 "df8e33389a8d9f15eb7eb17200c17002")
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/cgal/External_CGAL.cmake")
ENDIF(USE_CGAL)

#-----------------------------------------------------------------------------
# Find fast box intersection library.
#-----------------------------------------------------------------------------
IF(USE_LIBFBI)
  IF(WIN32)
    IF(NOT MINGW)
      SET(MSG "Only the latest versions of MSVC support the features needed")
      SET(MSG "${MSG} needed to build libfbi. Maybe your version works or")
      SET(MSG "${MSG} not. We just bail out here, to make sure nothing stupid")
      SET(MSG "${MSG} happens.")
      MESSAGE(FATAL_ERROR "${MSG}")
    ENDIF()
  ENDIF()

  SET(LIBFBI_URL "${LSE17_SOURCES_DIR}/spacepart")
  SET(LIBFBI_GZ "libfbi_for_cfs_gitrev_ee570e5e.tgz")
  SET(LIBFBI_MD5 "9484f3573450b20cadc262365eb74b7a")
ENDIF(USE_LIBFBI)
INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/spacepart/External_spacepart.cmake")

#-----------------------------------------------------------------------------
# FLANN - Fast Library for Approximate Nearest Neighbors
#-----------------------------------------------------------------------------
IF(USE_FLANN)
  SET(FLANN_URL "${LSE17_SOURCES_DIR}/flann")
  SET(FLANN_GZ "flann-1.8.4-src.zip")
  SET(FLANN_MD5 "a0ecd46be2ee11a68d2a7d9c6b4ce701")
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/flann/External_FLANN.cmake")
ENDIF(USE_FLANN)

#-----------------------------------------------------------------------------
# Find IPOPT
#-----------------------------------------------------------------------------
IF(USE_IPOPT)
  MESSAGE(FATAL_ERROR "IPOPT has not been ported to CMake externals yet.")
#  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindIPOPT.cmake")    
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
  MESSAGE(FATAL_ERROR "SnOpt has not been ported to CMake externals yet.")
#  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindSnOpt.cmake")    
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
  #---------------------------------------------------------------------------
  # Setup a list of dependencies for ParaView.
  #---------------------------------------------------------------------------
  SET(CFS_PV_DEPENDENCIES 
    hdf5-shared
    cgns
    boost
    zlib
    )

  #---------------------------------------------------------------------------
  # Qt - Let's check if a valid version of Qt is available
  #---------------------------------------------------------------------------
  FIND_PACKAGE(Qt4 4.8.2)
  IF(NOT QT4_FOUND)
    set(QT4_URL "${LSE17_SOURCES_DIR}/qt4")
    set(QT4_GZ qt-everywhere-opensource-src-4.8.2.tar.gz)
    set(QT4_MD5 3c1146ddf56247e16782f96910a8423b)

    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/qt4/External_Qt4.cmake")
  ENDIF()

  #---------------------------------------------------------------------------
  # Build QtCurve style.
  #---------------------------------------------------------------------------
  IF(UNIX)
    INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/qt4/External_QtCurve.cmake")
  ENDIF()

  #---------------------------------------------------------------------------
  # Finally add an external project for ParaViewSuperbuild
  #---------------------------------------------------------------------------
  set(PARAVIEW_SB_URL "${LSE17_SOURCES_DIR}/paraview")
  set(PARAVIEW_SB_GZ pvsuperbuild-4.1.0.tgz)
  set(PARAVIEW_SB_MD5 94e49d03e38955c0cb98c1159f417feb)
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/paraview/External_ParaViewSuperbuild.cmake")
ENDIF(BUILD_PARAVIEW)


#-----------------------------------------------------------------------------
# Find HDF file viewer
#-----------------------------------------------------------------------------
IF(BUILD_HDFVIEW)
  MESSAGE(FATAL_ERROR "HDFView has not been ported to CMake externals yet.")
#  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindHDFView.cmake")
ENDIF(BUILD_HDFVIEW)

IF(USE_PYTHON)
  set(PYTHON_MAJOR 2)
  set(PYTHON_MINOR 7)
  set(PYTHON_PATCH 2)
  set(PYTHON_VERSION ${PYTHON_MAJOR}.${PYTHON_MINOR}.${PYTHON_PATCH})
  set(PYTHON_URL "${LSE17_SOURCES_DIR}/python")
  set(PYTHON_GZ Python-${PYTHON_VERSION}.tgz)
  set(PYTHON_MD5 0ddfe265f1b3d0a8c2459f5bf66894c7)
  
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/python/External_Python.cmake")
ENDIF(USE_PYTHON)

#-------------------------------------------------------------------------------
# The cfsdeps meta target. Issue 'make -jX cfsdeps' to build all required.
# external projects. Especially important for parallel builds!
#-------------------------------------------------------------------------------
ADD_CUSTOM_TARGET(cfsdeps)
ADD_DEPENDENCIES(cfsdeps ${CFSDEPS})
