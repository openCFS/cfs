#-------------------------------------------------------------------------------
# PETSc (Portable, Extensible Toolkit for Scientific Computation)
# https://www.mcs.anl.gov/petsc/
# PETSc is a parallel libraray based on MPI only (no OpenMP!). 
# It hides the complexity of MPI to some extend but also 
# provide a large selection of algebraic solvers and preconditioners.
#
# PETSc is also efficient to be run on a single shared memory system!
# 
# See comments in FindCFSDEPS.cmake close to USE_MPI

set(PETSC_PREFIX  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/petsc")
set(PETSC_SOURCE  "${PETSC_PREFIX}/src/petsc")
set(PETSC_INSTALL  "${PETSC_PREFIX}/install")

# petsc is build by configure, therefore no cmake args needed
# no patches are required

SET(MIRRORS
  "http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/${PETSC_TGZ}"
  "${CFS_DS_SOURCES_DIR}/petsc/${PETSC_TGZ}")

SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/petsc/${PETSC_TGZ}")
SET(MD5_SUM ${PETSC_MD5})

SET(DLFN "${PETSC_PREFIX}/petsc-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

# After the installation we copy to cfs
SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/petsc/petsc-post_install.cmake.in")
SET(PI "${PETSC_PREFIX}/petsc-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

# we need to add our MPI_BASE (see mpi.cmake) and assume it common for CC, CXX and FC
PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "petsc" "${PETSC_VER}_${MPI_BASE}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${PETSC_INSTALL}")

SET(ZIPFROMCACHE "${PETSC_PREFIX}/petsc-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${PETSC_PREFIX}/petsc-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

# Determine paths of petsc libraries.
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
SET(PETSC_LIBRARY "${CFS_BINARY_DIR}/${LIB_SUFFIX}/libpetsc.so" CACHE FILEPATH "petsc library")
SET(PETSC_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

MARK_AS_ADVANCED(PETSC_LIBRARY)
MARK_AS_ADVANCED(PETSC_INCLUDE_DIR)

#-------------------------------------------------------------------------------
# The petsc external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(petsc
    PREFIX "${PETSC_PREFIX}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # special handling for mkl
  IF(USE_BLAS_LAPACK MATCHES "MKL")
    SET(PETSC_BLAS " --with-blas-lapack-dir=${MKL_ROOT}")
  ELSE()
    # shall work at least for openblas
    set(PETSC_BLAS " --with-blaslapack-lib=${LAPACK_LIB}")
  ENDIF()       
 
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    SET(PETSC_DEBUG 1)
  ELSE()
    SET(PETSC_DEBUG 0)
  ENDIF()
  ExternalProject_Add(petsc
    PREFIX "${PETSC_PREFIX}"
    SOURCE_DIR "${PETSC_SOURCE}"
    URL ${LOCAL_FILE}
    URL_MD5 ${PETSC_MD5}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    # petsc has a python2 configure, maybe later also python3 compatible?
    # python2 config/configure.py  CC=$CC CXX=$CXX FC=$FC --prefix=/home/fwein/tmp/petsc-3.8.3/killme
    # CC, .. shall be mpicc to be verified in mpi.cmake. Unfortunately there is no assert() in cmake :(
    #Coudn't avoid copy paste for configure command as the configure command messes up things quotes
    CONFIGURE_COMMAND python2 ${PETSC_SOURCE}/configure --with-cc=${CMAKE_C_COMPILER} --with-cxx=${CMAKE_CXX_COMPILER} --with-fc=${CMAKE_Fortran_COMPILER} --with-debugging=${PETSC_DEBUG} COPTFLAGS=-O3 CXXOPTFLAGS=-O3 FOPTFLAGS=-O3 --prefix=${PETSC_INSTALL} --with-valgrind=0 --with-blas-lapack-dir=${MKL_ROOT_DIR}
    INSTALL_COMMAND make install
    BUILD_BYPRODUCTS ${PETSC_LIBRARY}
  )

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(petsc cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${PETSC_PREFIX}
  )

  #-------------------------------------------------------------------------------
  # Execute the stuff from petsc-post_install.cmake after installation
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(petsc post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(petsc cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

# Add project to global list of CFSDEPS, this allows "make petsc"
SET(CFSDEPS ${CFSDEPS} petsc)
