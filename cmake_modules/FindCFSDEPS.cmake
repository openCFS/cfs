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
# DEPRECATED
#
# OPTIONS 
# 
# USAGE 
# 
# NOTES
#
# AUTHOR
# Simon Triebenbacher simon@ibtriebenbacher.de (07/2006)


IF(WIN32)
  SET(WIN32_STYLE_FIND 1)
ENDIF(WIN32)
IF(MINGW)
  SET(WIN32_STYLE_FIND 0)
  SET(UNIX_STYLE_FIND 1)
ENDIF(MINGW)
IF(UNIX)
  SET(UNIX_STYLE_FIND 1)
ENDIF(UNIX)

IF(NOT UNIX_STYLE_FIND) 
  MESSAGE(FATAL_ERROR "FindCFSDEPS.cmake:  Platform unknown/unsupported by FindCFSDEPS.cmake. It's not UNIX compatible.")
ENDIF(NOT UNIX_STYLE_FIND) 


#---------------------------------------------------------------------------
# Candidates for root/base directory of CFSDEPS
#---------------------------------------------------------------------------
SET (CFSDEPS_POSSIBLE_ROOT_PATHS
  "$ENV{CFSDEPS_ROOT}"
  "/home/data/libraries/CFSDEPS/${CFS_ARCH_STR}"
  "/opt/lse"
  )

#    MESSAGE("DBG CFSDEPS_POSSIBLE_ROOT_PATHS: ${CFSDEPS_POSSIBLE_ROOT_PATHS}")

#-------------------------------------------------------------------------------
# Find root/base directory of CFSDEPS among possible directories.
# Candidates for root/base directory of CFSDEPS
# should have subdirs include and lib containing include/hdf5.h
#-------------------------------------------------------------------------------
FIND_PATH(CFSDEPS_ROOT_DIR  include/hdf5.h 
  ${CFSDEPS_POSSIBLE_ROOT_PATHS} 
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH)  

#-------------------------------------------------------------------------------
# Set paths for include and lib directories.
#-------------------------------------------------------------------------------
IF(${CFS_ARCH} STREQUAL "X86_64")
  SET(CFSDEPS_LIBRARY_DIR "${CFSDEPS_ROOT_DIR}/lib64") 
ELSE(${CFS_ARCH} STREQUAL "X86_64")
  SET(CFSDEPS_LIBRARY_DIR "${CFSDEPS_ROOT_DIR}/lib") 
ENDIF(${CFS_ARCH} STREQUAL "X86_64")

SET(CFSDEPS_INCLUDE_DIR "${CFSDEPS_ROOT_DIR}/include")


MARK_AS_ADVANCED(
  CFSDEPS_ROOT_DIR
  CFSDEPS_INCLUDE_DIR
  CFSDEPS_LIBRARY_DIR
  )


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

#-------------------------------------------------------------------------------
# Find ACML library
#-------------------------------------------------------------------------------
IF(CFS_SUBARCH STREQUAL "OPTERON")
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindACML.cmake")
ELSE(CFS_SUBARCH STREQUAL "OPTERON")
  #-----------------------------------------------------------------------------
  # Find mkl library
  #-----------------------------------------------------------------------------
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindIntelMKL.cmake")
ENDIF(CFS_SUBARCH STREQUAL "OPTERON")

IF(NOT ACML_FOUND AND NOT MKL_FOUND)
  #-----------------------------------------------------------------------------
  # If USE_BLAS option is defined find BLAS library
  #-----------------------------------------------------------------------------
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindBLAS.cmake")
  
  #-----------------------------------------------------------------------------
  # If USE_LAPACK option is defined find LAPACK library
  #-----------------------------------------------------------------------------
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindLAPACK.cmake")
ENDIF(NOT ACML_FOUND AND NOT MKL_FOUND)

#-------------------------------------------------------------------------------
# Search for Pardiso Library
#-------------------------------------------------------------------------------
IF(NOT MKL_FOUND)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindPardiso.cmake")
ENDIF(NOT MKL_FOUND)

#-------------------------------------------------------------------------------
# If USE_ARPACK option is defined find ARPACK library
#-------------------------------------------------------------------------------
IF(USE_ARPACK)
  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindARPACK.cmake")
ENDIF(USE_ARPACK)

#-------------------------------------------------------------------------------
# Find Boost
#-------------------------------------------------------------------------------
INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindBoostForCFS.cmake")

#-------------------------------------------------------------------------------
# If USE_PYTHON option is defined find Python library
#-------------------------------------------------------------------------------
#IF(USE_PYTHON)
#  INCLUDE("${CFS_SOURCE_DIR}/cmake_modules/FindPythonForNACS.cmake")
#ENDIF(USE_PYTHON)

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


#-------------------------------------------------------------------------------
# Set BLAS_LIBRARY and LAPACK_LIBRARY according to architecture
#-------------------------------------------------------------------------------
IF(CFS_ARCH STREQUAL "I386")
  IF(MKL_FOUND)
    SET(BLAS_MKL_SERIAL_LIB "${MKL_ROOT_DIR}/lib_serial/32/libmkl_ia32.a" CACHE FILEPATH "MKL Blas.")
    SET(BLAS_MKL_OPENMP_LIB "${MKL_ROOT_DIR}/lib/32/libmkl_ia32.a" CACHE FILEPATH "MKL Blas.")

    SET(LAPACK_MKL_SERIAL_LIB "${MKL_ROOT_DIR}/lib_serial/32/libmkl_lapack.a" CACHE FILEPATH "MKL Lapack.")
    SET(LAPACK_MKL_OPENMP_LIB "${MKL_ROOT_DIR}/lib/32/libmkl_lapack.a" CACHE FILEPATH "MKL Lapack.")

    SET(PARDISO_MKL_SERIAL_LIB "${MKL_ROOT_DIR}/lib_serial/32/libmkl_solver.a" CACHE FILEPATH "MKL Pardiso.")
    SET(PARDISO_MKL_OPENMP_LIB "${MKL_ROOT_DIR}/lib/32/libmkl_solver.a" CACHE FILEPATH "MKL Pardiso.")

    SET(MKL_GUIDE_OPENMP_LIB "${MKL_ROOT_DIR}/lib/32/libguide.a")

    MARK_AS_ADVANCED(BLAS_MKL_SERIAL_LIB)
    MARK_AS_ADVANCED(BLAS_MKL_OPENMP_LIB)
    MARK_AS_ADVANCED(LAPACK_MKL_SERIAL_LIB)
    MARK_AS_ADVANCED(LAPACK_MKL_OPENMP_LIB)
    MARK_AS_ADVANCED(PARDISO_MKL_SERIAL_LIB)
    MARK_AS_ADVANCED(PARDISO_MKL_OPENMP_LIB)
    
    IF(USE_OPENMP)
      SET(BLAS_LIBRARY "${BLAS_MKL_OPENMP_LIB};${MKL_GUIDE_OPENMP_LIB}")
      SET(LAPACK_LIBRARY "${LAPACK_MKL_OPENMP_LIB}")
      SET(PARDISO_LIBRARY "${PARDISO_MKL_OPENMP_LIB}")
    ELSE(USE_OPENMP)
      SET(BLAS_LIBRARY "${BLAS_MKL_SERIAL_LIB}")
      SET(LAPACK_LIBRARY "${LAPACK_MKL_SERIAL_LIB}")
      SET(PARDISO_LIBRARY "${PARDISO_MKL_SERIAL_LIB}")
    ENDIF(USE_OPENMP)

    IF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
      SET(BLAS_LIBRARY "${BLAS_LIBRARY};-lpthread")
    ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")

  ELSE(MKL_FOUND)
    SET(BLAS_LIBRARY "${BLAS_STANDARD_LIB}")
    SET(LAPACK_LIBRARY "${LAPACK_STANDARD_LIB}")
    SET(PARDISO_LIBRARY "${PARDISO_STANDARD_LIB}")
  ENDIF(MKL_FOUND)
ENDIF(CFS_ARCH STREQUAL "I386")

IF(CFS_ARCH STREQUAL "X86_64")
  IF(CFS_SUBARCH STREQUAL "OPTERON")
    IF(ACML_FOUND)
      IF(USE_OPENMP)
	SET(BLAS_LIBRARY "${BLAS_ACML_OPENMP_LIB}")
	SET(LAPACK_LIBRARY "${LAPACK_ACML_OPENMP_LIB}")
      ELSE(USE_OPENMP)
	SET(BLAS_LIBRARY "${BLAS_ACML_SERIAL_LIB}")
	SET(LAPACK_LIBRARY "${LAPACK_ACML_SERIAL_LIB}")
      ENDIF(USE_OPENMP)
    ELSE(ACML_FOUND)
      SET(BLAS_LIBRARY ${BLAS_STANDARD_LIB})
      SET(LAPACK_LIBRARY ${LAPACK_STANDARD_LIB})
    ENDIF(ACML_FOUND)

    SET(PARDISO_LIBRARY "${PARDISO_STANDARD_LIB}")

  ELSE(CFS_SUBARCH STREQUAL "OPTERON")
    IF(MKL_FOUND)
      SET(BLAS_MKL_SERIAL_LIB "${MKL_ROOT_DIR}/lib_serial/em64t/libmkl_em64t.a" CACHE FILEPATH "MKL Blas.")
      SET(BLAS_MKL_OPENMP_LIB "${MKL_ROOT_DIR}/lib/em64t/libmkl_em64t.a" CACHE FILEPATH "MKL Blas.")
      
      SET(LAPACK_MKL_SERIAL_LIB "${MKL_ROOT_DIR}/lib_serial/em64t/libmkl_lapack.a" CACHE FILEPATH "MKL Lapack.")
      SET(LAPACK_MKL_OPENMP_LIB "${MKL_ROOT_DIR}/lib/em64t/libmkl_lapack.a" CACHE FILEPATH "MKL Lapack.")
      
      SET(PARDISO_MKL_SERIAL_LIB "${MKL_ROOT_DIR}/lib_serial/em64t/libmkl_solver.a" CACHE FILEPATH "MKL Pardiso.")
      SET(PARDISO_MKL_OPENMP_LIB "${MKL_ROOT_DIR}/lib/em64t/libmkl_solver.a" CACHE FILEPATH "MKL Pardiso.")

      SET(MKL_GUIDE_OPENMP_LIB "${MKL_ROOT_DIR}/lib/em64t/libguide.a")

      MARK_AS_ADVANCED(BLAS_MKL_SERIAL_LIB)
      MARK_AS_ADVANCED(BLAS_MKL_OPENMP_LIB)
      MARK_AS_ADVANCED(LAPACK_MKL_SERIAL_LIB)
      MARK_AS_ADVANCED(LAPACK_MKL_OPENMP_LIB)
      MARK_AS_ADVANCED(PARDISO_MKL_SERIAL_LIB)
      MARK_AS_ADVANCED(PARDISO_MKL_OPENMP_LIB)
      
      IF(USE_OPENMP)
	SET(BLAS_LIBRARY "${BLAS_MKL_OPENMP_LIB};${MKL_GUIDE_OPENMP_LIB}")
	SET(LAPACK_LIBRARY "${LAPACK_MKL_OPENMP_LIB}")
	SET(PARDISO_LIBRARY "${PARDISO_MKL_OPENMP_LIB}")
      ELSE(USE_OPENMP)
	SET(BLAS_LIBRARY "${BLAS_MKL_SERIAL_LIB}")
	SET(LAPACK_LIBRARY "${LAPACK_MKL_SERIAL_LIB}")
	SET(PARDISO_LIBRARY "${PARDISO_MKL_SERIAL_LIB}")
      ENDIF(USE_OPENMP)
      
      IF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
	SET(BLAS_LIBRARY "${BLAS_LIBRARY};-lpthread")
      ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
    ELSE(MKL_FOUND)
      SET(BLAS_LIBRARY "${BLAS_STANDARD_LIB}")
      SET(LAPACK_LIBRARY "${LAPACK_STANDARD_LIB}")
      SET(PARDISO_LIBRARY "${PARDISO_STANDARD_LIB}")
    ENDIF(MKL_FOUND)
  ENDIF(CFS_SUBARCH STREQUAL "OPTERON")
ENDIF(CFS_ARCH STREQUAL "X86_64")

IF(CFSDEPS_INCLUDE_DIR OR CFSDEPS_CXX_FLAGS)
  ## found all we need.
  SET(CFSDEPS_FOUND 1)
  SET(CFSDEPS_INCLUDE_DIR "${CFSDEPS_INCLUDE_DIR}"
    CACHE PATH INTERNAL "CFSDEPS_INCLUDE_DIR")
  
  SET(CFSDEPS_CXX_FLAGS "${CFSDEPS_CXX_FLAGS}"
    CACHE INTERNAL "CFSDEPS_CXX_FLAGS")
  
  SET(CFSDEPS_LIBRARIES "${CFSDEPS_LIBRARY_DIR}"
    CACHE INTERNAL "CFSDEPS_LIBRARY_DIR")
  
  #    MESSAGE("CFSDEPS_INCLUDE_DIR: ${CFSDEPS_INCLUDE_DIR}")
  #    MESSAGE("CFSDEPS_LIBRARY_DIR: ${CFSDEPS_LIBRARY_DIR}")
  #    MESSAGE("DBG found CFSDEPS_FOUND: ${CFSDEPS_FOUND}")

ENDIF(CFSDEPS_INCLUDE_DIR OR CFSDEPS_CXX_FLAGS)
