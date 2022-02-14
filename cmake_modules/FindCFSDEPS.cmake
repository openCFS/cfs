#============================================================================
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

SET(CFS_DS_SOURCES_DIR "${CFS_FAU_MIRROR}/sources")
SET(CFSDEPS_DIR "${CFS_SOURCE_DIR}/cfsdeps")

#-----------------------------------------------------------------------------
# Set common CFLAGS (and CXXFLAGS) common for all external projects.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# We do not want to see warnings from external projects, since they would
# show up on CDash.
#-----------------------------------------------------------------------------
if(CFS_CXX_COMPILER_NAME STREQUAL "GCC" OR CFS_CXX_COMPILER_NAME STREQUAL "CLANG")
  if(NOT CFS_OPT_FLAGS)
    message(STATUS "CFS_OPT_FLAGS not set, check order with compile.cmake")
  endif()
  set(CFSDEPS_C_FLAGS "${CFS_OPT_FLAGS} -w")
  set(CFSDEPS_CXX_FLAGS "${CFS_OPT_FLAGS} -w ${CFSDEPS_CXX_FLAGS}")
  if(USE_CGAL) # remove when we use header only CGAL
    set(CFSDEPS_C_FLAGS "-frounding-math ${CFSDEPS_C_FLAGS}")
    set(CFSDEPS_CXX_FLAGS "-frounding-math ${CFSDEPS_CXX_FLAGS}")
  endif()
endif()
if(CFS_FORTRAN_COMPILER_NAME STREQUAL "GCC" OR CFS_CXX_COMPILER_NAME STREQUAL "FLANG")
  set(CFSDEPS_Fortran_FLAGS "{CFS_OPT_FLAGS} -w")
endif()  

# TODO: Intel is missing but there is a lot CFSDEPS_ stuff for intel in compiler.cmake
#message(STATUS "CFS_OPT_FLAGS = ${CFS_OPT_FLAGS}")
#message(STATUS "CMAKE_COMPILER_IS_GNUCXX = ${CMAKE_COMPILER_IS_GNUCXX}")
#message(STATUS "CFS_CXX_COMPILER_NAME = ${CFS_CXX_COMPILER_NAME}")
#message(STATUS "CFS_FORTRAN_COMPILER_NAME = ${CFS_FORTRAN_COMPILER_NAME}")
#message(STATUS "CFSDEPS_CXX_FLAGS = ${CFSDEPS_CXX_FLAGS}")

# handle gfortran >= 10.
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU" AND (NOT ${FC_VERSION} VERSION_LESS 10))
  # was once --std=legacy
  # see https://github.com/Reference-LAPACK/lapack/issues/353
  set(CFSDEPS_Fortran_FLAGS "${CFSDEPS_Fortran_FLAGS} -fallow-argument-mismatch")
endif()  

if(MSVC)
  STRING(REPLACE "/W3" "/w" CFSDEPS_C_FLAGS "${CMAKE_C_FLAGS_INIT}")
  STRING(REPLACE "/W3" "/w" CFSDEPS_CXX_FLAGS "${CMAKE_CXX_FLAGS_INIT}")
  STRING(REPLACE "/W3" "/w" CFSDEPS_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_INIT}")
endif()

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
  find_program(CONFIGURE_MAKE_PROGRAM make)
else()
  set(CONFIGURE_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM} CACHE FILEPATH "program to build configure projects")
endif()
mark_as_advanced(CONFIGURE_MAKE_PROGRAM)

# The Python lib stuff is not built but found on the system when we embedd python
if(USE_EMBEDDED_PYTHON)
  # see find_package(PythonInterp) in FindPrograms.cmake
  
  # sets PYTHON_LIBRARY and PYTHON_INCLUDE_DIR, which can be also set via -D
  find_package(PythonLibs)

  # it might happen, when configured first w/o USE_EMBEDDED_PYTHON, that the cachend two variables are empty
  # and not overwritten by find_package - cmake is such a mess :( Setting via ccmake helps  
  # message(STATUS "PYTHON_INCLUDE_DIR = ${PYTHON_INCLUDE_DIR}")
  # message(STATUS "PYTHON_LIBRARY = ${PYTHON_LIBRARY}")
  
  # PYTHON_SITE_PACKAGES_DIR needs to be set to /usr/lib64/python3.8/site-packages
  # PYTHON_LIBRARY is /usr/lib64/libpython3.8.so, so it does not help
  # there is no uniform cmake python support. There are several FindPythonModule.cmake on the web but not in the standard distribution
  # here we extract the core from https://git.trustedfirmware.org/TF-M/trusted-firmware-m.git/tree/cmake/FindPythonModules.cmake
  # you need to include like #include <numpy/core/include/numpy/arrayobject.h>
  execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print(numpy.__file__);"
      RESULT_VARIABLE PY_NUMPY_STATUS # 0 is good, 1 is bad
      OUTPUT_VARIABLE PY_NUMPY_PATH # shall be /usr/lib64/python3.8/site-packages/numpy/__init__.py
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(NOT PY_NUMPY_STATUS) # 0 is good
    get_filename_component(PY_NUMPY_DIR ${PY_NUMPY_PATH} DIRECTORY) #  /usr/lib64/python3.8/site-packages/numpy  
    get_filename_component(PY_SITE_DIR ${PY_NUMPY_DIR} DIRECTORY) # /usr/lib64/python3.8/site-packages
  endif()

  # ignored if already set manually  
  set(PYTHON_SITE_PACKAGES_DIR ${PY_SITE_DIR} CACHE PATH "base for python site-packages")     
  mark_as_advanced(PYTHON_SITE_PACKAGES_DIR)

  # fun stuff: If you encounter the error message Pyhton.h not found, uncomment the message and try again :)
  # message(STATUS "add include_directories(${PYTHON_INCLUDE_DIR})")
  include_directories(${PYTHON_INCLUDE_DIR})

  if(NOT PYTHON_SITE_PACKAGES_DIR)
    message(FATAL_ERROR "cannot determine numpy for ${PYTHON_EXECUTABLE} required for USE_EMBEDDED_PYTHON. Possibly set advanced PYTHON_SITE_PACKAGES_DIR")
  endif()
endif()

#-------------------------------------------------------------------------------
# Build zlib library
#-------------------------------------------------------------------------------
SET(ZLIB_VER "1.2.8")
SET(ZLIB_GZ "zlib-${ZLIB_VER}.tar.gz")
SET(ZLIB_MD5 "44d667c142d7cda120332623eab69f40")

INCLUDE("${CFSDEPS_DIR}/zlib/External_zlib.cmake")

#-------------------------------------------------------------------------------
# Build bzip2 library
#-------------------------------------------------------------------------------
SET(BZIP2_URL "${CFS_DS_SOURCES_DIR}/bzip2")
SET(BZIP2_VER "1.0.6")
SET(BZIP2_GZ "bzip2-${BZIP2_VER}.tar.gz")
SET(BZIP2_MD5 "00b516f4704d4a7cb50a1d97e6e8e15b")

INCLUDE("${CFSDEPS_DIR}/bzip2/External_bzip2.cmake")

#-------------------------------------------------------------------------------
# Search for HDF5 library
#-------------------------------------------------------------------------------
  IF(APPLE)
    SET(HDF5_VER "1.8.17")
    SET(HDF5_MD5 "34bd1afa5209259201a41964100d6203") # 1.8.17
  ELSE()
    SET(HDF5_VER "1.8.12")
    SET(HDF5_MD5 "03ad766d225f5e872eb3e5ce95524a08")
  ENDIF()
  
  SET(HDF5_BZ2 "hdf5-${HDF5_VER}.tar.bz2")
  INCLUDE("${CFSDEPS_DIR}/hdf5/External_HDF5.cmake")

#-------------------------------------------------------------------------------
# Search for CGNS library
#-------------------------------------------------------------------------------
IF(USE_CGNS)
  SET(CGNS_VER "3.3.0")
  SET(CGNS_GZ "v${CGNS_VER}.tar.gz")
  SET(CGNS_MD5 "64e5e8d97144c1462bee9ea6b2a81d7f")

  INCLUDE("${CFSDEPS_DIR}/cgns/External_CGNS.cmake")
ENDIF(USE_CGNS)

#-------------------------------------------------------------------------------
# Search for METIS library
#-------------------------------------------------------------------------------
IF(USE_METIS)
  SET(METIS_VER "4.0.3")
  SET(METIS_GZ "metis-${METIS_VER}.tar.gz")
  SET(METIS_MD5 "d3848b454532ef18dc83e4fb160d1e10")

  INCLUDE("${CFSDEPS_DIR}/metis/External_METIS.cmake")
ENDIF(USE_METIS)

#-------------------------------------------------------------------------------
# Search for GiDpost library
#-------------------------------------------------------------------------------
IF(USE_GIDPOST)
  SET(GIDPOST_VER "2.1")
  SET(GIDPOST_ZIP "gidpost-${GIDPOST_VER}.zip")
  SET(GIDPOST_MD5 "a7fe745e40593dc4598920b312663d29")

  INCLUDE("${CFSDEPS_DIR}/gidpost/External_GiDpost.cmake")
ENDIF(USE_GIDPOST)

#-----------------------------------------------------------------------------
# Find Netlib BLAS/LAPACK library
# MKL contains blas and lapack, OpenBLAS contains blas and somehow also lapack?!
#-----------------------------------------------------------------------------
IF(USE_BLAS_LAPACK STREQUAL "NETLIB" OR USE_ILUPACK  )
  SET(LAPACK_VER "3.9.0")
  SET(LAPACK_GZ "v${LAPACK_VER}.tar.gz")
  SET(LAPACK_MD5 "0b251e2a8d5f949f99b50dd5e2200ee2")
    
  INCLUDE("${CFSDEPS_DIR}/lapack/External_LAPACK.cmake")
    
ENDIF(USE_BLAS_LAPACK STREQUAL "NETLIB"  OR USE_ILUPACK )

#-----------------------------------------------------------------------------
# Find OpenBLAS/LAPACK library
# see NETLIB comment
#-----------------------------------------------------------------------------
if(USE_BLAS_LAPACK STREQUAL "OPENBLAS")
  set(OPENBLAS_VER "0.3.17")
  set(OPENBLAS_GZ "v${OPENBLAS_VER}.tar.gz")
  set(OPENBLAS_MD5 "5429954163bcbaccaa13e11fe30ca5b6")
  INCLUDE("${CFSDEPS_DIR}/openblas/External_OpenBLAS.cmake")
endif(USE_BLAS_LAPACK STREQUAL "OPENBLAS")

#-----------------------------------------------------------------------------
# Find Intel Math Kernel library
# see NETLIB comment
#-----------------------------------------------------------------------------
IF(USE_BLAS_LAPACK STREQUAL "MKL")
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindIntelMKL.cmake")
ENDIF(USE_BLAS_LAPACK STREQUAL "MKL")

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
# Generally one needs to switch off USE_PARDISO with USE_BLAS_LAPACK not MKL!
#-----------------------------------------------------------------------------
IF(USE_PARDISO)
  INCLUDE("cmake_modules/CheckPardisoAPIVersion.cmake")
ENDIF(USE_PARDISO)
  
#-----------------------------------------------------------------------------
# If USE_ARPACK option is defined find ARPACK library
#-----------------------------------------------------------------------------
IF(USE_ARPACK)
  # remove --std=legacy in External_ARPACK when > 3.7.0
  SET(ARPACK_VER "3.7.0")
  SET(ARPACK_GZ "${ARPACK_VER}.tar.gz")
  SET(ARPACK_MD5 "6fc6c6bf78dbd4f144595ef0675c8430")
  
  INCLUDE("${CFSDEPS_DIR}/arpack/External_ARPACK.cmake")
ENDIF(USE_ARPACK)
  
#-----------------------------------------------------------------------------
# Find SuiteSparse/CholMod/UMFPACK/AMD library
#-----------------------------------------------------------------------------
IF(USE_SUITESPARSE OR USE_ILUPACK)
  SET(SUITESPARSE_VER "4.2.1")
  SET(SUITESPARSE_GZ "SuiteSparse-${SUITESPARSE_VER}.tar.gz")
  SET(SUITESPARSE_MD5 "4628df9eeae10ae5f0c486f1ac982fce")

  INCLUDE("${CFSDEPS_DIR}/suitesparse/External_SuiteSparse.cmake")
ENDIF()

#-----------------------------------------------------------------------------
# Find ILUPACK library
#-----------------------------------------------------------------------------
# TODO: skip ILUPACK_PARALLEL, it makes not really sense
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
  SET(ILUPACK_ZIP "ilupack-${ILUPACK_VER}_src.tgz")
  SET(ILUPACK_MD5 "0a5792597f8120d71e221de601440137")
  INCLUDE("${CFSDEPS_DIR}/ilupack/External_ILUPACK.cmake")
  #ADD_DEPENDENCIES(ilupack metis) # ilupack has its own (parallel) metis additional to cfs-metis
ELSEIF(USE_ILUPACK)
  SET(ILUPACK_VER "2.2.1")
  SET(ILUPACK_ZIP "ilupack-${ILUPACK_VER}-cfsdeps.zip")
  SET(ILUPACK_MD5 "79887591fb508013d874c97d4b252ed4")
  INCLUDE("${CFSDEPS_DIR}/ilupack/External_ILUPACK.cmake")
ENDIF()

#-------------------------------------------------------------------------------
# Find Library of Iterative Solvers
#-------------------------------------------------------------------------------
IF(USE_LIS)
  SET(LIS_VER "2.0.24")
  SET(LIS_ZIP "lis-${LIS_VER}.zip")
  SET(LIS_MD5 "1c7a84d3204dcefb3af2627eba58e97f")
  
  INCLUDE("${CFSDEPS_DIR}/lis/External_LIS.cmake")
ENDIF(USE_LIS)

#-------------------------------------------------------------------------------
# Find Library for FFT and IFFT 
#-------------------------------------------------------------------------------
IF(USE_FFTW)
  SET(FFTW_VER "3.3.4")
  SET(FFTW_GZ "fftw-${FFTW_VER}.tar.gz")
  SET(FFTW_MD5 "2edab8c06b24feeb3b82bbb3ebf3e7b3")
  
  INCLUDE("${CFSDEPS_DIR}/fftw/External_FFTW.cmake")
ENDIF(USE_FFTW)

#-------------------------------------------------------------------------------
# Find SuperLU
#-------------------------------------------------------------------------------
IF(USE_SUPERLU)
  SET(SUPERLU_VER "4.3")
  SET(SUPERLU_GZ "superlu_${SUPERLU_VER}.tar.gz")
  SET(SUPERLU_MD5 "b72c6309f25e9660133007b82621ba7c")
  
  INCLUDE("${CFSDEPS_DIR}/superlu/External_SuperLU.cmake")
ENDIF(USE_SUPERLU)

#-------------------------------------------------------------------------------
# Find Boost
#-------------------------------------------------------------------------------
SET(BOOST_MAJOR_VER 1)
SET(BOOST_MINOR_VER 73)
SET(BOOST_VER "${BOOST_MAJOR_VER}.${BOOST_MINOR_VER}")
SET(BOOST_GZ "boost_${BOOST_MAJOR_VER}_${BOOST_MINOR_VER}_0.tar.bz2")
SET(BOOST_MD5 "b2dfbd6c717be4a7bb2d88018eaccf75") #1.66
#SET(BOOST_MD5 "7fbd1890f571051f2a209681d57d486a") # 1.68
SET(BOOST_MD5 "9273c8c4576423562bbe84574b07b2bd") # 1.73
INCLUDE("${CFSDEPS_DIR}/boost/External_Boost.cmake")

#-------------------------------------------------------------------------------
# Build MuParser library
#-------------------------------------------------------------------------------
# https://github.com/beltoforion/muparser/archive/v2.2.6.1.tar.gz
# revision 388b3f9 contains a necessary feature request not available in 2.2.6.1
# someone shall switch to a real revision once strfunc4-5 are there
# https://codeload.github.com/beltoforion/muparser/zip/388b3f9
SET(MUPARSER_VER "2.2.6.1")
SET(MUPARSER_TGZ "v${MUPARSER_VER}.tar.gz") 
SET(MUPARSER_MD5 "410d29b4c58d1cdc2fc9ed1c1c7f67fe") # 2.2.6.1

INCLUDE("${CFSDEPS_DIR}/muparser/External_muParser.cmake")

#-------------------------------------------------------------------------------
# Xerces library or libxml2, triggered by CFS_XML_READER
#-------------------------------------------------------------------------------
IF(USE_XERCES)
  SET(XERCES_VER "3.1.3")
  SET(XERCES_GZ "xerces-c-${XERCES_VER}.tar.gz")
  SET(XERCES_MD5 "70320ab0e3269e47d978a6ca0c0e1e2d") 
  INCLUDE("${CFSDEPS_DIR}/xerces/External_Xerces-C.cmake")
ENDIF()

#-------------------------------------------------------------------------------
# libxml2 is an alternative for Xerces
#-------------------------------------------------------------------------------
IF(USE_LIBXML2)
  SET(LIBXML2_VER "2.9.4")
  SET(LIBXML2_GZ "libxml2-${LIBXML2_VER}.tar.gz")
  SET(LIBXML2_MD5 "ae249165c173b1ff386ee8ad676815f5") 
  INCLUDE("${CFSDEPS_DIR}/libxml2/External_LibXml2.cmake")
ENDIF()

#-----------------------------------------------------------------------------
# Find VTK - used for Ensight only
#-----------------------------------------------------------------------------
IF(USE_VTK)
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
    MESSAGE(FATAL_ERROR "${MSG}")
   ENDIF()

  SET(GMP_VER "6.1.2")
  SET(GMP_BZ2 "gmp-${GMP_VER}.tar.bz2")
  SET(GMP_MD5 "8ddbb26dc3bd4e2302984debba1406a5")
  INCLUDE("${CFSDEPS_DIR}/gmp/External_gmp.cmake")
  
  SET(MPFR_VER "3.1.5")
  SET(MPFR_BZ2 "mpfr-${MPFR_VER}.tar.bz2")
  SET(MPFR_MD5 "b1d23a55588e3b2a13e3be66bc69fd8d")
  INCLUDE("${CFSDEPS_DIR}/mpfr/External_mpfr.cmake")

  SET(CGAL_VER "4.2")
  SET(CGAL_BZ2 "CGAL-${CGAL_VER}.tar.bz2")
  SET(CGAL_MD5 "df8e33389a8d9f15eb7eb17200c17002")
  INCLUDE("${CFSDEPS_DIR}/cgal/External_CGAL.cmake")
ENDIF(USE_CGAL)

#-----------------------------------------------------------------------------
# Find fast box intersection library.
#-----------------------------------------------------------------------------
IF(USE_LIBFBI)
  IF(WIN32)
    SET(MSG "Only the latest versions of MSVC support the features needed")
    SET(MSG "${MSG} needed to build libfbi. Maybe your version works or")
    SET(MSG "${MSG} not. We just bail out here, to make sure nothing stupid")
    SET(MSG "${MSG} happens.")
    MESSAGE(FATAL_ERROR "${MSG}")
  ENDIF()

  SET(LIBFBI_VER "for_cfs_gitrev_ee570e5e")
  SET(LIBFBI_GZ "libfbi_${LIBFBI_VER}.tgz")
  SET(LIBFBI_MD5 "9484f3573450b20cadc262365eb74b7a")
ENDIF(USE_LIBFBI)
INCLUDE("${CFSDEPS_DIR}/spacepart/External_spacepart.cmake")

#-----------------------------------------------------------------------------
# FEAST - FEAST Eigenvalue Solver
#-----------------------------------------------------------------------------
IF(USE_FEAST_COMMUNITY)
  SET(FEAST_VER "3.0") # note that this is ignored in feast/CMakeLists.txt
  SET(FEAST_GZ "feast_${FEAST_VER}.tgz")
  SET(FEAST_MD5 "f03819c19a8724d0095dd24eae7ba43a")
  INCLUDE("${CFSDEPS_DIR}/feast/External_FEAST.cmake")
ENDIF()

#-----------------------------------------------------------------------------
# FLANN - Fast Library for Approximate Nearest Neighbors
#-----------------------------------------------------------------------------
IF(USE_FLANN)
  SET(FLANN_VER "1.9.1")
  SET(FLANN_ZIP "flann-${FLANN_VER}.zip")
  SET(FLANN_MD5 "4a6cc62db8ed09dd8a0c6537f6720f12")
  INCLUDE("${CFSDEPS_DIR}/flann/External_FLANN.cmake")
ENDIF(USE_FLANN)

#-----------------------------------------------------------------------------
# Find SCPIP - A special optimizer for topology optimization. 
# This is not open source, so check with Christoph Zillober, Uni-Wuerzburg first 
#-----------------------------------------------------------------------------
if(USE_SCPIP)
  set(SCPIP_ZIP "scpip-cfsdeps.zip")
  set(SCPIP_MD5 "a2767521596ed925e53b9fea6df4d77e")  
  include("${CFSDEPS_DIR}/scpip/External_SCPIP.cmake")
endif(USE_SCPIP)

#-----------------------------------------------------------------------------
# Find SNOPT - A commercial general purpose commercial optimizer 
#-----------------------------------------------------------------------------
if(USE_SNOPT)
  # you need to have a license for the source. 
  # with keys for CFS_DOWNLOAD_SNOPT and CFS_KEY_SNOPT a special licsence can be used
  set(SNOPT_VER "7.2.8")
  set(SNOPT_MD5 "9e75be8400eb878b9cb3d489084af196")
  set(SNOPT_ZIP "snopt-${SNOPT_VER}-cfsdeps.zip")
  include("${CFSDEPS_DIR}/snopt/External_SNOPT.cmake")
endif(USE_SNOPT)

#-----------------------------------------------------------------------------
# Find IPOPT - A general purpos open source optimizer 
#-----------------------------------------------------------------------------
IF(USE_IPOPT)
  SET(IPOPT_VER "3.14.2")
  SET(IPOPT_ZIP "Ipopt-${IPOPT_VER}.zip")
  SET(IPOPT_MD5 "83b94837024ef3a30688cf094cff86d7") 
 
  INCLUDE("${CFSDEPS_DIR}/ipopt/External_IPOPT.cmake")
ENDIF(USE_IPOPT)

#-----------------------------------------------------------------------------
# Find SGPP - A toolbox for sparse grid interpolation 
#-----------------------------------------------------------------------------
IF(USE_SGPP)
  # this version is outdated, the current one is opensource https://sgpp.sparsegrids.org/
  message(FATAL_ERROR "switch to open source sgpp")

  SET(SGPP_VER "2016-03-04_166a3d9")
  SET(SGPP_ZIP "sgopt_${SGPP_VER}.zip")
  INCLUDE("${CFSDEPS_DIR}/sgpp/External_SGPP.cmake")
ENDIF(USE_SGPP)

#-----------------------------------------------------------------------------
# Find HDF file viewer
#-----------------------------------------------------------------------------
IF(BUILD_HDFVIEW)
  MESSAGE(FATAL_ERROR "HDFView has not been ported to CMake externals yet.")
ENDIF(BUILD_HDFVIEW)

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
  SET(HWLOC_MD5 "a0fa1c9109a4d8b4b6568e62cc9b6e30") # 1.11.8 
  
  INCLUDE("${CFSDEPS_DIR}/hwloc/External_HWLOC.cmake")
endif(BUILD_HWLOC)

# ghost is required for phist or could be used standalone
if(BUILD_GHOST)
  # we use the cfs-fork of ghost and download the stuff via bitbuket
  # we could also use a subversion mirror on github but only for ghost, not for g
  # svn co https://github.com/RRZE-HPC/GHOST/trunk@r<REVSION>
  # set(GHOST_REV "965be2d1aa20") # subversion revision numbers are are more easily handable :(
  # set(GHOST_MD5 "1f441c4b82aaf0e9ff507857ffbd8c0e")
  set(GHOST_REV "0e54b108ada8") # does not work with cfs (06.03.20)
  set(GHOST_MD5 "b0eca2287f12c0dd2c5fda03a99c7de4")
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
  #set(PHIST_REV "8a22be1e42aa") 
  #set(PHIST_MD5 "076a7bc70040a375f285c3e9fee2112d")
  set(PHIST_REV "1935ec0accc1")  # does not work with cfs (06.03.20)
  set(PHIST_MD5 "cfec4ad70a3838894742ac7732f56070")
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

