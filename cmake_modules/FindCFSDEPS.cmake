#============================================================================
#
# Find locations of external binary libs (e.g. MKL) and build additional
# external libs from source.
# 
# This module finds and builds libs that openCFS depends upon and determines
# where the include files and libraries are. 
#
#=============================================================================

include(ExternalProject) # cmake external project
include("cmake_modules/DependencyTools.cmake") # our own helper for cfsdeps handling (pseudo object oriented)

SET(CFS_DS_SOURCES_DIR "${CFS_FAU_MIRROR}/sources")
SET(CFSDEPS_DIR "${CFS_SOURCE_DIR}/cfsdeps")

# CFSDEPS_*_FLAGS are set in compiler.cmake

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
which can be reused for other openCFS builds.
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
include("${CFSDEPS_DIR}/zlib/External_zlib.cmake")


#-------------------------------------------------------------------------------
# Search for HDF5 library
#-------------------------------------------------------------------------------
include("${CFSDEPS_DIR}/hdf5/External_HDF5.cmake")

#-------------------------------------------------------------------------------
# Search for CGNS library
#-------------------------------------------------------------------------------
IF(USE_CGNS)
  SET(CGNS_VER "4.3.0")
  SET(CGNS_GZ "v${CGNS_VER}.zip")
  SET(CGNS_MD5 "37512acaac66c368b454658dc8a806ff")

  INCLUDE("${CFSDEPS_DIR}/cgns/External_CGNS.cmake")
ENDIF(USE_CGNS)

if(USE_METIS)
  include("${CFSDEPS_DIR}/metis/External_METIS.cmake")
endif()

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
if(USE_BLAS_LAPACK STREQUAL "NETLIB")
  include("${CFSDEPS_DIR}/netlib/External_Netlib.cmake")
endif()

#-----------------------------------------------------------------------------
# Find OpenBLAS/LAPACK library
# see NETLIB comment
#-----------------------------------------------------------------------------
if(USE_BLAS_LAPACK STREQUAL "OPENBLAS")
  include("${CFSDEPS_DIR}/openblas/External_OpenBLAS.cmake")
endif()

#-----------------------------------------------------------------------------
# Find Intel Math Kernel library
# see NETLIB comment
#-----------------------------------------------------------------------------
if(USE_BLAS_LAPACK STREQUAL "MKL")
  include("${CFS_SOURCE_DIR}/cmake_modules/FindIntelMKL.cmake")
endif()

# Apple's Accelerate Framework is a system lib - nothing to build
if(USE_BLAS_LAPACK STREQUAL "ACCELERATE")
  include("${CFS_SOURCE_DIR}/cmake_modules/FindAppleAccelerate.cmake")
endif()


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
  
if(USE_ARPACK)
  include("${CFSDEPS_DIR}/arpack/External_ARPACK.cmake")
endif()
  
if(USE_SUITESPARSE)
  include("${CFSDEPS_DIR}/suitesparse/External_SuiteSparse.cmake")
endif()

if(USE_EIGEN)
  include("${CFSDEPS_DIR}/eigen/External_EIGEN.cmake")
endif()

# Find Library of Iterative Solvers
if(USE_LIS)
  include("${CFSDEPS_DIR}/lis/External_LIS.cmake")
endif(USE_LIS)

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
include("${CFSDEPS_DIR}/boost/External_Boost.cmake")

#-------------------------------------------------------------------------------
# Build MuParser library
#-------------------------------------------------------------------------------
include("${CFSDEPS_DIR}/muparser/External_muParser.cmake")

#-------------------------------------------------------------------------------
# Xerces library or libxml2, triggered by CFS_XML_READER
#-------------------------------------------------------------------------------
if(USE_XERCES)
  include("${CFSDEPS_DIR}/xerces/External_Xerces-C.cmake")
endif()

#-------------------------------------------------------------------------------
# libxml2 is an alternative for Xerces
#-------------------------------------------------------------------------------
if(USE_LIBXML2)
  include("${CFSDEPS_DIR}/libxml2/External_LibXml2.cmake")
endif()

#-----------------------------------------------------------------------------
# Find VTK - used for Ensight only
#-----------------------------------------------------------------------------
if(USE_VTK)
  include("${CFS_SOURCE_DIR}/cfsdeps/vtk/External_VTK.cmake")
endif(USE_VTK)

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

  SET(CGAL_URL "${CFS_DS_SOURCES_DIR}/cgal")
  SET(CGAL_VER "4.9.1")
  SET(CGAL_BZ2 "CGAL-${CGAL_VER}.tar.xz")
  SET(CGAL_MD5 "820ef17ffa7ed87af6cc9918a961d966")
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
if(USE_FEAST)
  include("${CFSDEPS_DIR}/feast/External_FEAST.cmake")
endif()

# FLANN - Fast Library for Approximate Nearest Neighbors
if(USE_FLANN)
  include("${CFSDEPS_DIR}/flann/External_FLANN.cmake")
endif()

#-----------------------------------------------------------------------------
# Find SCPIP - A special optimizer for topology optimization. 
# This is not open source, so check with Christoph Zillober, Uni-Wuerzburg first 
#-----------------------------------------------------------------------------
if(USE_SCPIP)
  include("${CFSDEPS_DIR}/scpip/External_SCPIP.cmake")
endif()

#-----------------------------------------------------------------------------
# Find SNOPT - A commercial general purpose commercial optimizer 
#-----------------------------------------------------------------------------
if(USE_SNOPT)
  include("${CFSDEPS_DIR}/snopt/External_SNOPT.cmake")
endif()

#-----------------------------------------------------------------------------
# Find Ipopt - A general purpose open source optimizer 
#-----------------------------------------------------------------------------
if(USE_IPOPT)
  if(WIN32)
    include("${CFSDEPS_DIR}/ipopt/External_IPOPT_Win.cmake")
  else()
    include("${CFSDEPS_DIR}/ipopt/External_IPOPT.cmake")
    # if no MKL is available, use HSL via ThirdParty-HSL. creates dynamic coinlib for Ipopt
    if(NOT(USE_BLAS_LAPACK STREQUAL "MKL"))
      include("${CFSDEPS_DIR}/hsl/External_HSL.cmake")
    endif()
  endif()
endif()

#-----------------------------------------------------------------------------
# Find external SGP (Sequential Global Programming) library - An open source optimizer for structural optimization problems
# Currently: Simulatneous topology + local orientation optimization 
#-----------------------------------------------------------------------------
if(USE_SGP)
  include("${CFSDEPS_DIR}/sgp/External_SGP.cmake")
endif()

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

