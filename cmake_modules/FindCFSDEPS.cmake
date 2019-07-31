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
# Include external project build capability of CMake 2.8, cf.
# http://www.kitware.com/media/html/BuildingExternalProjectsWithCMake2.8.html
#-----------------------------------------------------------------------------
INCLUDE(ExternalProject)

#-----------------------------------------------------------------------------
# Include external data capability of CMake 2.8.11, cf.
# http://www.kitware.com/source/home/post/107
#-----------------------------------------------------------------------------
INCLUDE(ExternalData)

SET(CFS_DS_SOURCES_DIR "${CFS_DS_CFSDEPS}/cfsdeps/sources")
SET(CFSDEPS_DIR "${CFS_SOURCE_DIR}/cfsdeps")

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
  SET(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -isysroot ${CMAKE_OSX_SYSROOT}")

  SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -arch ${CMAKE_OSX_ARCHITECTURES}")
  SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -isysroot ${CMAKE_OSX_SYSROOT}")

  IF(CXX_HAS_SYSROOT_FLAG)
    SET(CFSDEPS_C_FLAGS "${CFSDEPS_C_FLAGS} -sysroot=${CMAKE_OSX_SYSROOT}")
    SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -sysroot=${CMAKE_OSX_SYSROOT}")
  ENDIF()
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

FILE(TO_CMAKE_PATH "${CFS_DEPS_CACHE_DIR}" CFS_DEPS_CACHE_DIR)

# for configure projectes we may not use ninja but need make
if("${CMAKE_GENERATOR}" STREQUAL "Ninja")
  message("CMAKE_GENERATOR = ${CMAKE_GENERATOR}")
  find_program(CONFIGURE_MAKE_PROGRAM make)
else()
  set(CONFIGURE_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM} CACHE FILEPATH "program to build configure projects")
endif()
mark_as_advanced(CONFIGURE_MAKE_PROGRAM)


#-------------------------------------------------------------------------------
# Build zlib library
#-------------------------------------------------------------------------------
SET(ZLIB_URL "${CFS_DS_SOURCES_DIR}/zlib")
SET(ZLIB_BASE "zlib")
SET(ZLIB_VER "1.2.8")
SET(ZLIB_GZ "${ZLIB_BASE}-${ZLIB_VER}.tar.gz")
SET(ZLIB_MD5 "44d667c142d7cda120332623eab69f40")

INCLUDE("${CFSDEPS_DIR}/zlib/External_zlib.cmake")

#-------------------------------------------------------------------------------
# Build bzip2 library
#-------------------------------------------------------------------------------
SET(BZIP2_URL "${CFS_DS_SOURCES_DIR}/bzip2")
SET(BZIP2_BASE "bzip2")
SET(BZIP2_VER "1.0.6")
SET(BZIP2_GZ "${BZIP2_BASE}-${BZIP2_VER}.tar.gz")
SET(BZIP2_MD5 "00b516f4704d4a7cb50a1d97e6e8e15b")

INCLUDE("${CFSDEPS_DIR}/bzip2/External_bzip2.cmake")

#-------------------------------------------------------------------------------
# Search for HDF5 library
#-------------------------------------------------------------------------------
IF(USE_HDF5)
  SET(HDF5_URL "${CFS_DS_SOURCES_DIR}/hdf5")
  SET(HDF5_BASE "hdf5")
  IF(APPLE)
    # macOS 10.12 requires gcc as clang does not catch exceptions. gcc comes with brew as 6.2.0 which
    # is not able to compile hdf5-1.8.12 but works with 1.8.17. However 1.8.17 requires cmake >= 3.0
    # which is unconvenient for many CFS developers and also requires the following changes for the mingw
    # cross-compiler (Windows on Linux)  
    # - hdf5-cross-compile.patch shall be replaced by hdf5-cross-compile.hdf5-1.8.17.patch (mosty case sensitive stuff)
    # - hdf5-mingw.patch becomes obsolete as the patch now comes from upstream
    # - TryRun* needs to be checked for the prefix (H5 -> HDF5) ...
    SET(HDF5_VER "1.8.17")
    SET(HDF5_MD5 "34bd1afa5209259201a41964100d6203") # 1.8.17
  ELSE()
    SET(HDF5_VER "1.8.12")
    SET(HDF5_MD5 "03ad766d225f5e872eb3e5ce95524a08")
  ENDIF()

  SET(HDF5_BZ2 "${HDF5_BASE}-${HDF5_VER}.tar.bz2")
  INCLUDE("${CFSDEPS_DIR}/hdf5/External_HDF5.cmake")
ENDIF(USE_HDF5)

#-------------------------------------------------------------------------------
# Search for CGNS library
#-------------------------------------------------------------------------------
IF(USE_CGNS)
  SET(CGNS_URL "${CFS_DS_SOURCES_DIR}/cgns")
  SET(CGNS_BASE "v")
  SET(CGNS_VER "3.3.0")
  SET(CGNS_GZ "${CGNS_BASE}${CGNS_VER}.tar.gz")
  SET(CGNS_MD5 "64e5e8d97144c1462bee9ea6b2a81d7f")

  INCLUDE("${CFSDEPS_DIR}/cgns/External_CGNS.cmake")
ENDIF(USE_CGNS)

#-------------------------------------------------------------------------------
# Search for STAR-CCM+ I/O library
#-------------------------------------------------------------------------------
IF(USE_CCMIO)
  SET(CCMIO_URL "${CFS_DS_SOURCES_DIR}/ccmio")
  SET(CCMIO_BASE "libccmio")
  SET(CCMIO_VER "2.6.1")
  SET(CCMIO_GZ "${CCMIO_BASE}-${CCMIO_VER}.tar.gz")
  SET(CCMIO_MD5 "f81fbdfb960b1a4f3bcc7feee491efe4")

  INCLUDE("${CFSDEPS_DIR}/ccmio/External_CCMIO.cmake")
ENDIF(USE_CCMIO)

#-------------------------------------------------------------------------------
# Search for CFX I/O library
#-------------------------------------------------------------------------------
IF(USE_CFXIO)
  INCLUDE("${CFSDEPS_DIR}/cfx_custom/External_CFX_Custom.cmake")
ENDIF(USE_CFXIO)

#-------------------------------------------------------------------------------
# Search for METIS library
#-------------------------------------------------------------------------------
IF(USE_METIS)
  SET(METIS_URL "${CFS_DS_SOURCES_DIR}/metis")
  SET(METIS_BASE "metis")
  SET(METIS_VER "4.0.3")
  SET(METIS_GZ "${METIS_BASE}-${METIS_VER}.tar.gz")
  SET(METIS_MD5 "d3848b454532ef18dc83e4fb160d1e10")

  INCLUDE("${CFSDEPS_DIR}/metis/External_METIS.cmake")
ENDIF(USE_METIS)

#-------------------------------------------------------------------------------
# Search for GiDpost library
#-------------------------------------------------------------------------------
IF(USE_GIDPOST)
  SET(GIDPOST_URL "${CFS_DS_SOURCES_DIR}/gidpost")
  SET(GIDPOST_BASE "gidpost")
  SET(GIDPOST_VER "2.1")
  SET(GIDPOST_ZIP "${GIDPOST_BASE}-${GIDPOST_VER}.zip")
  SET(GIDPOST_MD5 "a7fe745e40593dc4598920b312663d29")

  INCLUDE("${CFSDEPS_DIR}/gidpost/External_GiDpost.cmake")
ENDIF(USE_GIDPOST)

#-----------------------------------------------------------------------------
# Find Netlib BLAS/LAPACK library
# MKL contains blas and lapack, OpenBLAS contains blas and somehow also lapack?!
#-----------------------------------------------------------------------------
IF(CFS_BLAS_LAPACK STREQUAL "NETLIB" OR USE_ILUPACK  )
    
  SET(LAPACK_URL "${CFS_DS_SOURCES_DIR}/lapack")
  SET(LAPACK_BASE "lapack")
  SET(LAPACK_VER "3.4.2")
  SET(LAPACK_GZ "${LAPACK_BASE}-${LAPACK_VER}.tgz")
  SET(LAPACK_MD5 "61bf1a8a4469d4bdb7604f5897179478")
    
  INCLUDE("${CFSDEPS_DIR}/lapack/External_LAPACK.cmake")
    
ENDIF(CFS_BLAS_LAPACK STREQUAL "NETLIB"  OR USE_ILUPACK )

#-----------------------------------------------------------------------------
# Find OpenBLAS/LAPACK library
# see NETLIB comment
#-----------------------------------------------------------------------------
if(CFS_BLAS_LAPACK STREQUAL "OPENBLAS")
    
  set(OPENBLAS_URL "${CFS_DS_SOURCES_DIR}/openblas")
  set(OPENBLAS_BASE "OpenBLAS")
  # the latest revision 0.2.20 does not compile with -DDYNAMIC_ARCH=1 but the current development version does
  # set(OPENBLAS_VER "0.2.20")
  # set(OPENBLAS_GZ "v${OPENBLAS_VER}.tar.gz")
  set(OPENBLAS_REV "874df65") # 1.2.2018
  set(OPENBLAS_ZIP "${OPENBLAS_REV}.zip")
  # this is the filename on https://github.com/xianyi/OpenBLAS/archive, the sourceforge link is with spaces
  set(OPENBLAS_MD5 "0594f38346c725b5c88fcb648c4af1e8")
  INCLUDE("${CFSDEPS_DIR}/openblas/External_OpenBLAS.cmake")
endif(CFS_BLAS_LAPACK STREQUAL "OPENBLAS")

#-----------------------------------------------------------------------------
# Find Intel Math Kernel library
# see NETLIB comment
#-----------------------------------------------------------------------------
IF(CFS_BLAS_LAPACK STREQUAL "MKL")
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindIntelMKL.cmake")
ENDIF(CFS_BLAS_LAPACK STREQUAL "MKL")

#-----------------------------------------------------------------------------
# Check which version of the Pardiso API is being used. Pardiso 4.0 intro-
# duces a new incompatible API due to the new iterative solver feature.
# However MKL still uses the old API of Pardiso 3.x. This is obviously a
# software engineering nightmare and the ones responsible for it deserve to
# be hanged from the highest mast!
#
# ATTENTION: If you switch CFS_PARDISO you should also clear
# PARDISO_API_VER_3 and PARDISO_API_VER_4 from the CMake cache.
# 
# The non-MKL Pardiso version hasn't been used for quite a while. 
# Generally one needs to switch off USE_PARDISO with CFS_BLAS_LAPACK not MKL!
#-----------------------------------------------------------------------------
IF(USE_PARDISO)
  INCLUDE("cmake_modules/CheckPardisoAPIVersion.cmake")
ENDIF(USE_PARDISO)
  
#-----------------------------------------------------------------------------
# If USE_ARPACK option is defined find ARPACK library
#-----------------------------------------------------------------------------
IF(USE_ARPACK)
  SET(ARPACK_URL "${CFS_DS_SOURCES_DIR}/arpack")
  SET(ARPACK_BASE "arpack")
  SET(ARPACK_VER "ng-3.2.0")
  SET(ARPACK_GZ "${ARPACK_BASE}-${ARPACK_VER}.tar.gz")
  SET(ARPACK_MD5 "0ae8a0bb796370b06647d9e005c0f3ea")
  
  INCLUDE("${CFSDEPS_DIR}/arpack/External_ARPACK.cmake")
ENDIF(USE_ARPACK)
  
#-----------------------------------------------------------------------------
# Find SuiteSparse/CholMod/UMFPACK/AMD library
#-----------------------------------------------------------------------------
IF(USE_SUITESPARSE OR USE_ILUPACK)
  SET(SUITESPARSE_URL "${CFS_DS_SOURCES_DIR}/suitesparse")
  SET(SUITESPARSE_BASE "SuiteSparse")
  SET(SUITESPARSE_VER "4.2.1")
  SET(SUITESPARSE_GZ "${SUITESPARSE_BASE}-${SUITESPARSE_VER}.tar.gz")
  SET(SUITESPARSE_MD5 "4628df9eeae10ae5f0c486f1ac982fce")

  INCLUDE("${CFSDEPS_DIR}/suitesparse/External_SuiteSparse.cmake")
ENDIF(USE_SUITESPARSE OR USE_ILUPACK)

#-----------------------------------------------------------------------------
# Find ILUPACK library
#-----------------------------------------------------------------------------
IF(USE_ILUPACK_PARALLEL)
  #Since the latest version of ilupack requires GCC > 5.0 or the latest ICC compilers
  IF((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND 
   ((CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.0") OR (CMAKE_C_COMPILER_VERSION VERSION_LESS "5.0") OR (CFS_FORTRAN_COMPILER_VER VERSION_LESS "5.0") ))
    MESSAGE(FATAL_ERROR "Ilupack can be compiled only when gcc,g++ and gfortran compiler versions are greater than 5")
  ENDIF()
  # TODO: For intel compilers still one needs to figure out the proper compiler versions
# 
  SET(ILUPACK_PATH "${CFS_BINARY_DIR}/cfsdeps/ilupack")
  SET(ILUPACK_VER "2.4_parallel_0831")
  SET(ILUPACK_GZ "ilupack-${ILUPACK_VER}_src.tgz")
  SET(ILUPACK_MD5 "0a5792597f8120d71e221de601440137")
  INCLUDE("${CFSDEPS_DIR}/ilupack/External_ILUPACK.cmake")
  #ADD_DEPENDENCIES(ilupack metis) # ilupack has its own (parallel) metis additional to cfs-metis
ELSEIF(USE_ILUPACK)
  SET(ILUPACK_PATH "${CFS_BINARY_DIR}/cfsdeps/ilupack")
  SET(ILUPACK_BASE "ilupack")
  SET(ILUPACK_VER "2.2.1")
  SET(ILUPACK_GZ "${ILUPACK_BASE}${ILUPACK_VER}_src.tgz")
  SET(ILUPACK_MD5 "7cb6ba2e854e13d243218d9e9478d13c")
  INCLUDE("${CFSDEPS_DIR}/ilupack/External_ILUPACK.cmake")
ENDIF()




#  MESSAGE("BLAS_LIBRARY ${BLAS_LIBRARY}")
#  MESSAGE("LAPACK_LIBRARY ${LAPACK_LIBRARY}")
#  MESSAGE("PARDISO_LIBRARY ${PARDISO_LIBRARY}")
  

#-------------------------------------------------------------------------------
# Find Library of Iterative Solvers
#-------------------------------------------------------------------------------
IF(USE_LIS)
  SET(LIS_URL "${CFS_DS_SOURCES_DIR}/lis")
  SET(LIS_BASE "lis")
  SET(LIS_VER "1.7.28")
  SET(LIS_GZ "${LIS_BASE}-${LIS_VER}.tar.gz")
  SET(LIS_MD5 "59fbaea19fbf443d9a991ebacd7e58cd")
  
  INCLUDE("${CFSDEPS_DIR}/lis/External_LIS.cmake")
ENDIF(USE_LIS)

#-------------------------------------------------------------------------------
# Find Library for FFT and IFFT 
#-------------------------------------------------------------------------------
IF(USE_FFTW)
  SET(FFTW_URL "${CFS_DS_SOURCES_DIR}/fftw")
  SET(FFTW_BASE "fftw")
  SET(FFTW_VER "3.3.4")
  SET(FFTW_GZ "${FFTW_BASE}-${FFTW_VER}.tar.gz")
  SET(FFTW_MD5 "2edab8c06b24feeb3b82bbb3ebf3e7b3")
  
  INCLUDE("${CFSDEPS_DIR}/fftw/External_FFTW.cmake")
ENDIF(USE_FFTW)

#-------------------------------------------------------------------------------
# Find SuperLU
#-------------------------------------------------------------------------------
IF(USE_SUPERLU)
  SET(SUPERLU_URL "${CFS_DS_SOURCES_DIR}/superlu")
  SET(SUPERLU_BASE "superlu")
  SET(SUPERLU_VER "4.3")
  SET(SUPERLU_GZ "${SUPERLU_BASE}_${SUPERLU_VER}.tar.gz")
  SET(SUPERLU_MD5 "b72c6309f25e9660133007b82621ba7c")
  
  INCLUDE("${CFSDEPS_DIR}/superlu/External_SuperLU.cmake")
ENDIF(USE_SUPERLU)

#-------------------------------------------------------------------------------
# Find Boost
#-------------------------------------------------------------------------------
SET(BOOST_BASE "boost")
SET(BOOST_MAJOR_VER 1)
SET(BOOST_MINOR_VER 66)
SET(BOOST_VER "${BOOST_MAJOR_VER}.${BOOST_MINOR_VER}")
SET(BOOST_URL "${CFS_DS_SOURCES_DIR}/boost")
SET(BOOST_GZ "${BOOST_BASE}_${BOOST_MAJOR_VER}_${BOOST_MINOR_VER}_0.tar.bz2")
SET(BOOST_MD5 "b2dfbd6c717be4a7bb2d88018eaccf75") #1.66
#SET(BOOST_MD5 "7fbd1890f571051f2a209681d57d486a") # 1.68
INCLUDE("${CFSDEPS_DIR}/boost/External_Boost.cmake")

#-------------------------------------------------------------------------------
# Build MuParser library
#-------------------------------------------------------------------------------
SET(MUPARSER_URL "${CFS_DS_SOURCES_DIR}/muparser")
SET(MUPARSER_BASE "muparser")
# https://github.com/beltoforion/muparser/archive/v2.2.6.1.tar.gz
# revision 388b3f9 contains a necessary feature request not available in 2.2.6.1
# someone shall switch to a real revision once strfunc4-5 are there
# https://codeload.github.com/beltoforion/muparser/zip/388b3f9
SET(MUPARSER_VER "2.2.6.1")
# SET(MUPARSER_VER "388b3f9")
SET(MUPARSER_TGZ "v${MUPARSER_VER}.tar.gz") # # 2.2.6.1
#SET(MUPARSER_TGZ "${MUPARSER_VER}") # no extension for the revision
SET(MUPARSER_MD5 "410d29b4c58d1cdc2fc9ed1c1c7f67fe") # 2.2.6.1
#SET(MUPARSER_MD5 "07fbf24e44e8d94399dfb3bda4f454ba") # 388b3f9

INCLUDE("${CFSDEPS_DIR}/muparser/External_muParser.cmake")

#-------------------------------------------------------------------------------
# Xerces library or libxml2, triggered by XML_READER
#-------------------------------------------------------------------------------
IF(USE_XERCES)
  SET(XERCES_URL "${CFS_DS_SOURCES_DIR}/xerces")
  SET(XERCES_BASE "xerces")
  SET(XERCES_VER "3.1.3")
  SET(XERCES_GZ "${XERCES_BASE}-c-${XERCES_VER}.tar.gz")
  SET(XERCES_MD5 "70320ab0e3269e47d978a6ca0c0e1e2d") 
  INCLUDE("${CFSDEPS_DIR}/xerces/External_Xerces-C.cmake")
ENDIF(USE_XERCES)

#-------------------------------------------------------------------------------
# libxml2 is an alternative for Xerces
#-------------------------------------------------------------------------------
IF(USE_LIBXML2)
  SET(LIBXML2_VER "2.9.4")
  SET(LIBXML2_GZ "libxml2-${LIBXML2_VER}.tar.gz")
  SET(LIBXML2_MD5 "ae249165c173b1ff386ee8ad676815f5") 
  INCLUDE("${CFSDEPS_DIR}/libxml2/External_LibXml2.cmake")
ENDIF(USE_LIBXML2)

#-----------------------------------------------------------------------------
# Find VTK - used for Ensight only
#-----------------------------------------------------------------------------
IF(USE_VTK)
  SET(VTK_URL "${CFS_DS_SOURCES_DIR}/vtk")
  SET(VTK_TAR "VTK-7.1.1.tar.gz")
  SET(VTK_VERSION "7.1") # must be 2 digits for include-dir to be correct
  SET(VTK_MD5 "daee43460f4e95547f0635240ffbc9cb")
  INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/vtk/External_VTK.cmake")
ENDIF(USE_VTK)

#-----------------------------------------------------------------------------
# Find CGAL
#-----------------------------------------------------------------------------
IF(USE_CGAL)
  SET(MSG "The build of gmp and mpfr is only supported for MSYS on Windows!")
  SET(MSG "${MSG} It is configure-based and therefore requires a shell")
  SET(MSG "${MSG} interpreter like bash from MSYS. If you need CGAL, you need")
  SET(MSG "${MSG} to use an MSYS environment or cross compile from Linux.")     

  IF(WIN32)
    IF(MINGW AND NOT CMAKE_CROSSCOMPILING OR MSVC)
      IF(NOT $ENV{MSYSTEM} STREQUAL "MINGW32")
        MESSAGE(FATAL_ERROR "${MSG}")
      ENDIF()
    ENDIF()
  ENDIF()

  SET(GMP_URL "${CFS_DS_SOURCES_DIR}/gmp")
  SET(GMP_BASE "gmp")
  SET(GMP_VER "6.1.2")
  SET(GMP_BZ2 "${GMP_BASE}-${GMP_VER}.tar.bz2")
  SET(GMP_MD5 "8ddbb26dc3bd4e2302984debba1406a5")
  INCLUDE("${CFSDEPS_DIR}/gmp/External_gmp.cmake")
  
  SET(MPFR_URL "${CFS_DS_SOURCES_DIR}/mpfr")
  SET(MPFR_BASE "mpfr")
  SET(MPFR_VER "3.1.5")
  SET(MPFR_BZ2 "${MPFR_BASE}-${MPFR_VER}.tar.bz2")
  SET(MPFR_MD5 "b1d23a55588e3b2a13e3be66bc69fd8d")
  INCLUDE("${CFSDEPS_DIR}/mpfr/External_mpfr.cmake")

  SET(CGAL_URL "${CFS_DS_SOURCES_DIR}/cgal")
  SET(CGAL_BASE "CGAL")
  SET(CGAL_VER "4.2")
  SET(CGAL_BZ2 "${CGAL_BASE}-${CGAL_VER}.tar.bz2")
  SET(CGAL_MD5 "df8e33389a8d9f15eb7eb17200c17002")
  INCLUDE("${CFSDEPS_DIR}/cgal/External_CGAL.cmake")
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

  SET(LIBFBI_URL "${CFS_DS_SOURCES_DIR}/spacepart")
  SET(LIBFBI_BASE "libfbi")
  SET(LIBFBI_VER "for_cfs_gitrev_ee570e5e")
  SET(LIBFBI_GZ "${LIBFBI_BASE}_${LIBFBI_VER}.tgz")
  SET(LIBFBI_MD5 "9484f3573450b20cadc262365eb74b7a")
ENDIF(USE_LIBFBI)
INCLUDE("${CFSDEPS_DIR}/spacepart/External_spacepart.cmake")

#-----------------------------------------------------------------------------
# FEAST - FEAST Eigenvalue Solver
#-----------------------------------------------------------------------------
IF(USE_FEAST_COMMUNITY)
  SET(FEAST_URL "${CFS_DS_SOURCES_DIR}/feast")
  SET(FEAST_BASE "feast")
  SET(FEAST_VER "3.0")
  SET(FEAST_GZ "${FEAST_BASE}_${FEAST_VER}.tgz")
  SET(FEAST_MD5 "f03819c19a8724d0095dd24eae7ba43a")
  INCLUDE("${CFSDEPS_DIR}/feast/External_FEAST.cmake")
ENDIF()

#-----------------------------------------------------------------------------
# FLANN - Fast Library for Approximate Nearest Neighbors
#-----------------------------------------------------------------------------
IF(USE_FLANN)
  SET(FLANN_URL "${CFS_DS_SOURCES_DIR}/flann")
  SET(FLANN_BASE "flann")
  SET(FLANN_VER "1.9.1")
  SET(FLANN_ZIP "${FLANN_BASE}-${FLANN_VER}.zip")
  SET(FLANN_MD5 "4a6cc62db8ed09dd8a0c6537f6720f12")
  INCLUDE("${CFSDEPS_DIR}/flann/External_FLANN.cmake")
ENDIF(USE_FLANN)

#-----------------------------------------------------------------------------
# Find SCPIP - A special optimizer for topology optimization. 
# This is not open source, so check with Christoph Zillober, Uni-Wuerzburg first 
#-----------------------------------------------------------------------------
IF(USE_SCPIP)
  SET(SCPIP_URL "${CFS_DS_SOURCES_DIR}/scpip")
  SET(SCPIP_PATH "${CFS_BINARY_DIR}/cfsdeps/scpip")
  SET(SCPIP_BASE "scpip")
  SET(SCPIP_VER "")
  SET(SCPIP_BZ2 "${SCPIP_BASE}.tar.bz2")
  SET(SCPIP_MD5 "8afaf8d8d79981d68b8c726ea508471d")  
  INCLUDE("${CFSDEPS_DIR}/scpip/External_SCPIP.cmake")
ENDIF(USE_SCPIP)


#-----------------------------------------------------------------------------
# Find SNOPT - A general purpose commercial optimizer 
#-----------------------------------------------------------------------------
IF(USE_SNOPT)
  SET(SNOPT_URL "${CFS_DS_SOURCES_DIR}/snopt")
  SET(SNOPT_BASE "snopt")
  SET(SNOPT_VER "7.2.8")
  SET(SNOPT_MD5 "9e75be8400eb878b9cb3d489084af196")
  SET(SNOPT_PATH "${CFS_BINARY_DIR}/cfsdeps/snopt")
  # as snopt is commerical propriatary copy we hand an encrypted source file
  # you need to provide CFS_KEY_SNOPT with the password
  SET(SNOPT_ZIP "${SNOPT_BASE}-${SNOPT_VER}-cfsdeps.zip")
  INCLUDE("${CFSDEPS_DIR}/snopt/External_SNOPT.cmake")
ENDIF(USE_SNOPT)

#-----------------------------------------------------------------------------
# Find IPOPT - A general purpos open source optimizer 
#-----------------------------------------------------------------------------
IF(USE_IPOPT)
  SET(IPOPT_PATH "${CFS_BINARY_DIR}/cfsdeps/ipopt")
  SET(IPOPT_BASE "Ipopt")
  SET(IPOPT_VER "3.11.9")
  SET(IPOPT_TGZ "${IPOPT_BASE}-${IPOPT_VER}.tgz")
  SET(IPOPT_MD5 "657fa0f2f301f0d7b2a4e5b43e2370f5") 
 
  INCLUDE("${CFSDEPS_DIR}/ipopt/External_IPOPT.cmake")
ENDIF(USE_IPOPT)

#-----------------------------------------------------------------------------
# Find SGPP - A toolbox for sparse grid interpolation 
#-----------------------------------------------------------------------------
IF(USE_SGPP)
  SET(SGPP_PATH "${CFS_BINARY_DIR}/cfsdeps/sgpp")
  SET(SGPP_BASE "sgopt")
  SET(SGPP_VER "2016-03-04_166a3d9")
  SET(SGPP_ZIP "${SGPP_BASE}_${SGPP_VER}.zip")
  INCLUDE("${CFSDEPS_DIR}/sgpp/External_SGPP.cmake")
ENDIF(USE_SGPP)

#-------------------------------------------------------------------------------
# Build ZeroMQ distributed messaging library
#-------------------------------------------------------------------------------
IF(USE_ZEROMQ)
  SET(ZEROMQ_URL "${CFS_DS_SOURCES_DIR}/zeromq")
  SET(ZEROMQ_GZ "zeromq-4.1.3.tar.gz")
  SET(ZEROMQ_MD5 "d0824317348cfb44b8692e19cc73dc3a")

  INCLUDE("${CFSDEPS_DIR}/zeromq/External_zeromq.cmake")
ENDIF(USE_ZEROMQ)


#-----------------------------------------------------------------------------
# Find ANSYS Customizations
#-----------------------------------------------------------------------------
IF(USE_ANSYSRST)
  IF(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR WIN32)
    INCLUDE("${CFSDEPS_DIR}/ansys_custom/External_ANSYS_custom.cmake")
  ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR WIN32)
ENDIF(USE_ANSYSRST)

#-----------------------------------------------------------------------------
# Find HDF file viewer
#-----------------------------------------------------------------------------
IF(BUILD_HDFVIEW)
  MESSAGE(FATAL_ERROR "HDFView has not been ported to CMake externals yet.")
ENDIF(BUILD_HDFVIEW)

#-------------------------------------
# External anaconda 3 as a service for test machines at TU-Wien
#-------------------------------------
if(BUILD_ANACONDA3)
  SET(ANACONDA3_URL "${CFS_DS_SOURCES_DIR}/anaconda3")
  SET(ANACONDA3_SH "Anaconda3-4.2.0-Linux-x86_64.sh")
  SET(ANACONDA3_MD5 "4692f716c82deb9fa6b59d78f9f6e85c")
  INCLUDE("${CFSDEPS_DIR}/anaconda3/External_anaconda3.cmake")
endif(BUILD_ANACONDA3)

# PETSc requires mpi
if(USE_PETSC)
  SET(PETSC_VER "3.9.3")
  SET(PETSC_TGZ "petsc-${PETSC_VER}.tar.gz")
  SET(PETSC_MD5 "7b71d705f66f9961cb0e2da3f9da79a1")
  
  INCLUDE("${CFSDEPS_DIR}/petsc/External_PETSC.cmake")
endif(USE_PETSC)

# hwloc is a build dependency for ghost/phist but not explicitly used, therefore BUILD_HWLOC 
if(BUILD_HWLOC)
  SET(HWLOC_VER "1.11.8") # note that 1.11 is hardcoded in External_HWLOC!
  SET(HWLOC_TGZ "hwloc-${HWLOC_VER}.tar.gz")
  SET(HWLOC_MD5 "a0fa1c9109a4d8b4b6568e62cc9b6e30") 
  
  INCLUDE("${CFSDEPS_DIR}/hwloc/External_HWLOC.cmake")
endif(BUILD_HWLOC)

# ghost is required for phist or could be used standalone
if(BUILD_GHOST)
  # we use the cfs-fork of ghost and download the stuff via bitbuket
  # we could also use a subversion mirror on github but only for ghost, not for phist
  # svn co https://github.com/RRZE-HPC/GHOST/trunk@r<REVSION>
  set(GHOST_REV "a3b75fc52c7e") # subversion revision numbers are are more easily handable :(
  set(GHOST_MD5 "47a7f3b21cf8e8f16c69a6889a184f60")
  set(GHOST_ZIP "${GHOST_REV}.zip")
  # https://bitbucket.org/fabian_wein/cfs_ghost/get/840f2717f849.zip -> fabian_wein-cfs_ghost-840f2717f849
  # https://bitbucket.org/essex/ghost/get/f3c78b57e836.zip -> essex-ghost-f3c78b57e836
  set(GHOST_BB_USER "essex")
  set(GHOST_BB_PROJECT "ghost")
  include("${CFSDEPS_DIR}/ghost/External_GHOST.cmake")
  
  ADD_DEPENDENCIES(ghost hwloc)
endif(BUILD_GHOST)

# phist provides a ghost (=cuda if available) based EV-solver
if(USE_PHIST_EV OR USE_PHIST_CG)
  set(PHIST_REV "853be38a7078") 
  set(PHIST_MD5 "a07fcc37d45a49d5e2e8fc2b31f334a7")
  set(PHIST_ZIP "${PHIST_REV}.zip")
  set(PHIST_BB_USER "essex")
  set(PHIST_BB_PROJECT "phist")
  include("${CFSDEPS_DIR}/phist/External_PHIST.cmake")
  
  ADD_DEPENDENCIES(phist ghost)
endif()


#-------------------------------------------------------------------------------
# The cfsdeps meta target. Issue 'make -jX cfsdeps' to build all required.
# external projects. Especially important for parallel builds!
#-------------------------------------------------------------------------------
ADD_CUSTOM_TARGET(cfsdeps)
ADD_DEPENDENCIES(cfsdeps ${CFSDEPS})
