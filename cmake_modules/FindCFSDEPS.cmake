# This module finds and builds libs that openCFS depends from source with individual configuration

include(ExternalProject) # cmake external project
include("cmake_modules/DependencyTools.cmake") # our own helper for cfsdeps handling (pseudo object oriented)

if(NOT CMAKE_GENERATOR MATCHES "NMake Makefiles")
  message(STATUS "Most cfsdeps will be build with ${CFS_DEPS_BUILD_THREADS} threads.")
endif()

set(CFS_DS_SOURCES_DIR "${CFS_FAU_MIRROR}/sources")
set(CFSDEPS_DIR "${CFS_SOURCE_DIR}/cfsdeps")

# CFSDEPS_*_FLAGS are set in compiler.cmake

# If user has set environment variables use them. If not use defaults
set(CFS_DEPS_CD_DUMMY "$ENV{CFS_DEPS_CACHE_DIR}")
if(NOT ${CFS_DEPS_CD_DUMMY} STREQUAL "")
  set(CFS_DEPS_CACHE_DIR "${CFS_DEPS_CD_DUMMY}" CACHE PATH "Directory for CFSDEPS sources and prebuilt binaries.")
else()
  set(CFS_DEPS_CACHE_DIR "${CFS_DEPS_CACHE_DIR_DEFAULT}" CACHE PATH "Directory for CFSDEPS sources and prebuilt binaries.")
endif()

assert_set(CFS_DEPS_CACHE_DIR)

file(TO_CMAKE_PATH "${CFS_DEPS_CACHE_DIR}" CFS_DEPS_CACHE_DIR)

# The Python lib stuff is not built but linked with the system python when we embed python
if(USE_EMBEDDED_PYTHON)
  # we don't need the executable for embedded python - the testsuite finds it's own executable
  # use -DPython_ROOT_DIR=... to help in case or check cmake's FindPython for other hints
  # we need at least Python 3.8 but only with cmake 3.19 this requirement can be enforced
  set(Python_FIND_VIRTUALENV "FIRST")
  find_package(Python COMPONENTS Development NumPy)  
  #dump_variables("Python")
  if(NOT Python_Development_FOUND OR NOT Python_NumPy_FOUND)
    message(FATAL_ERROR "cannot find Python environment for USE_EMBEDDED_PYTHON, probably set Python_ROOT_DIR") 
  endif()
  message(STATUS "will link to embedded python ${Python_VERSION} using NumPy ${Python_NumPy_VERSION}")

  # we make use of these three variables
  assert_set(Python_INCLUDE_DIRS)
  assert_set(Python_NumPy_INCLUDE_DIRS)
  assert_set(Python_LIBRARIES) # e.g. /opt/homebrew/opt/python@3.11/Frameworks/Python.framework/Versions/3.11/lib/libpython3.11.dylib
  #cmake_print_variables(Python_INCLUDE_DIRS Python_NumPy_INCLUDE_DIRS Python_LIBRARIES)
  include_directories(${Python_INCLUDE_DIRS})
  include_directories(${Python_NumPy_INCLUDE_DIRS})
endif()

# these are optional external blas/lapack libs. 
# MKL and Apple's Accelerate is set in openmp_blas.cmake to have it before compile.cmake
if(USE_BLAS_LAPACK STREQUAL "NETLIB")
  set(USE_NETLIB 1)
  include("${CFSDEPS_DIR}/netlib/External_Netlib.cmake")
endif()

if(USE_BLAS_LAPACK STREQUAL "OPENBLAS")
  set(USE_OPENBLAS 1)
  include("${CFSDEPS_DIR}/openblas/External_OpenBLAS.cmake")
endif()

# preCICE is built from source via cfsdeps, UNLESS the developer points to an
# already-installed preCICE with -Dprecice_DIR=... (then the adapter falls back
# to find_package). Building preCICE pulls in eigen + libxml2 as build deps.
# Decided BEFORE zlib/boost/libxml2: those switch their build strategy
# (shared / -fPIC) when preCICE is built - see their External_*.cmake.
set(CFS_BUILD_PRECICE OFF)
if(USE_PRECICE AND NOT precice_DIR)
  set(CFS_BUILD_PRECICE ON)
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
if(USE_CGNS)
  include("${CFSDEPS_DIR}/cgns/External_CGNS.cmake")
endif(USE_CGNS)

if(USE_METIS)
  include("${CFSDEPS_DIR}/metis/External_METIS.cmake")
endif()

if(USE_GIDPOST)
  include("${CFSDEPS_DIR}/gidpost/External_GiDpost.cmake")
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

# preCICE optional MPI / PETSc features. openCFS does not build MPI itself and builds
# PETSc via cfsdeps only when USE_PETSC is set; preCICE's MPI/PETSc features are built
# against whatever openCFS/the system provides. Both flags only take EFFECT when preCICE
# is actually built (USE_PRECICE -> CFS_BUILD_PRECICE).
#
# USE_PRECICE_MPI is exposed UNCONDITIONALLY as a visible cache option (even while
# preCICE/OpenFOAM are off) so its state is always inspectable in ccmake - only the bool
# is always present, the forcing/defaulting logic below stays nested. It is tied to
# OpenFOAM, which is always built against a system MPI (WM_MPLIB=SYSTEMOPENMPI, whose
# presence is enforced in cfsdeps/openfoam/External_OpenFOAM.cmake): while USE_OPENFOAM
# is on the flag is forced ON and locked (and USE_OPENFOAM anyway requires USE_PRECICE,
# see the FATAL gate further below), otherwise it is a free choice defaulting to openCFS'
# own USE_MPI. _USE_PRECICE_MPI_OF_LOCK records OpenFOAM's last state so a user's own
# choice is preserved across plain reconfigures and only reset on the on->off transition.
#
# USE_PRECICE_PETSC stays nested (only shown/used when preCICE is built). PETSc mapping
# requires MPI, so enabling it forces USE_PRECICE_MPI on.

# --- MPI (always exposed) ---
option(USE_PRECICE_MPI "Build preCICE with MPI communication (against openCFS/system MPI); only effective when USE_PRECICE is on" "${USE_MPI}")
if(USE_OPENFOAM)
  set(USE_PRECICE_MPI ON CACHE BOOL "Build preCICE with MPI communication (against openCFS/system MPI); only effective when USE_PRECICE is on" FORCE)
  message(STATUS "preCICE: USE_PRECICE_MPI forced ON by USE_OPENFOAM (OpenFOAM always builds against a system MPI)")
elseif(_USE_PRECICE_MPI_OF_LOCK)
  # OpenFOAM was enabled (and forced MPI on) but is now off -> restore the free default
  set(USE_PRECICE_MPI "${USE_MPI}" CACHE BOOL "Build preCICE with MPI communication (against openCFS/system MPI); only effective when USE_PRECICE is on" FORCE)
endif()
set(_USE_PRECICE_MPI_OF_LOCK "${USE_OPENFOAM}" CACHE INTERNAL "tracks USE_OPENFOAM to flip USE_PRECICE_MPI back when OpenFOAM is disabled")

# --- PETSc (nested: only relevant when preCICE is built) ---
if(CFS_BUILD_PRECICE)
  set(_precice_petsc_default OFF)
  if(USE_PETSC)
    # openCFS builds PETSc in cfsdeps (see below) - reuse that one
    set(_precice_petsc_default ON)
  else()
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(_PRECICE_PETSC QUIET PETSc)
      if(_PRECICE_PETSC_FOUND)
        set(_precice_petsc_default ON)
      endif()
    endif()
  endif()
  option(USE_PRECICE_PETSC "Build preCICE with PETSc mapping (requires MPI; against openCFS/system PETSc)" ${_precice_petsc_default})

  # PETSc mapping needs MPI - force MPI on when PETSc is requested
  if(USE_PRECICE_PETSC AND NOT USE_PRECICE_MPI)
    message(STATUS "preCICE: USE_PRECICE_PETSC requires MPI - enabling USE_PRECICE_MPI")
    set(USE_PRECICE_MPI ON CACHE BOOL "Build preCICE with MPI communication (against openCFS/system MPI); only effective when USE_PRECICE is on" FORCE)
  endif()
endif()

# eigen is also a (header-only) build dependency of preCICE
if(USE_EIGEN OR CFS_BUILD_PRECICE)
  include("${CFSDEPS_DIR}/eigen/External_EIGEN.cmake")
endif()

# Find Library of Iterative Solvers
if(USE_LIS)
  include("${CFSDEPS_DIR}/lis/External_LIS.cmake")
endif()

if(USE_GINKGO)
  # we use json for the optional json configuration. Make it more universal once use want to use json independently
  include("${CFSDEPS_DIR}/nlohmann_json/External_nlohmann_json.cmake")
  include("${CFSDEPS_DIR}/ginkgo/External_Ginkgo.cmake")
endif()

if(USE_SUPERLU)
  include("${CFSDEPS_DIR}/superlu/External_SuperLU.cmake")
endif()

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
# libxml2 is an alternative for Xerces (and a required build dependency of preCICE)
#-------------------------------------------------------------------------------
if(USE_LIBXML2 OR CFS_BUILD_PRECICE)
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
if(USE_CGAL)
  # currently (CGAL 6.0.1) we still need gmp/mpfr on Unix, but it should work with boost alone
  if(NOT WIN32)
    include("${CFSDEPS_DIR}/gmp/External_GMP.cmake")
    include("${CFSDEPS_DIR}/mpfr/External_MPFR.cmake")
  endif()
  include("${CFSDEPS_DIR}/cgal/External_CGAL.cmake")
endif(USE_CGAL)

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

# nanoflann can efficiently do a k-nearest neighbor search
if(USE_NANOFLANN)
  include("${CFSDEPS_DIR}/nanoflann/External_NANOFLANN.cmake")
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
# C++ implementations for the topology optimiter MMA (dumas_mma) and a glonalized variant (dumas_gcmma)
#-----------------------------------------------------------------------------
if(USE_DUMAS)
  include("${CFSDEPS_DIR}/dumas/External_Dumas.cmake")
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

# preCICE coupling library (built against the cfs boost/eigen/libxml2 above).
# Exposes PRECICE_LIBRARY / PRECICE_INCLUDE_DIR consumed by source/Utils/preciceAdapter.
# Skipped when -Dprecice_DIR=... is given (developer-provided preCICE via find_package).
if(CFS_BUILD_PRECICE)
  include("${CFSDEPS_DIR}/precice/External_PRECICE.cmake")
endif()

# OpenFOAM + preCICE OpenFOAM adapter: the partner solver for real coupled
# openCFS<->OpenFOAM tests. Build order per precice.org: dependencies ->
# preCICE -> OpenFOAM -> OpenFOAM adapter (enforced via add_dependencies).
# Test-only (GPL): run as external programs, never linked or shipped with cfs.
if(USE_OPENFOAM)
  if(NOT CFS_BUILD_PRECICE)
    message(FATAL_ERROR "USE_OPENFOAM requires USE_PRECICE (without precice_DIR): "
      "the OpenFOAM adapter is built against the cfsdeps preCICE")
  endif()
  include("${CFSDEPS_DIR}/openfoam/External_OpenFOAM.cmake")
  include("${CFSDEPS_DIR}/openfoam-adapter/External_OpenFOAMAdapter.cmake")
endif()

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

