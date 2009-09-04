SET (MKL_POSSIBLE_PATHS
  ${MKL_ROOT_DIR}
  $ENV{MKL_ROOT_DIR}
  ${MKL_ROOT_DIR_DEFAULT}
  # Paths at LSE and local paths
  /opt/intel/Compiler/11.0/059/Frameworks/mkl
  /opt/intel/Compiler/11.0/081/mkl
  /opt/intel/Compiler/11.0/074/mkl
  /opt/intel/mkl/10.0.5.025
  /home/data/programs/intel/mkl/10.0.5.025
  /opt/intel/mkl/9.1.023
  /opt/intel/mkl/9.1.021
  /home/data/programs/intel/mkl/9.1.021
  # Paths on Woody
  /apps/intel/Compiler/11.0/069/mkl
  /apps/intel/mkl/10.0.011
  /apps/intel/ict/3.0/cmkl/9.0
 )

FIND_FILE(MKL_ROOT_DIR
  NAMES include/mkl.h
  PATHS ${MKL_POSSIBLE_PATHS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  )

STRING(REPLACE "/include/mkl.h" "" MKL_ROOT_DIR "${MKL_ROOT_DIR}")

IF(NOT MKL_ROOT_DIR)
  MESSAGE(FATAL_ERROR "Could not find Intel MKL library. This script just "
    "checks explicitly for versions 9.1.02[1|3] and 10.0.5.025 of MKL. "
    "Please either install one of those versions into /opt/intel/mkl, set "
    "an environment variable MKL_ROOT_DIR, set '-DMKL_ROOT_DIR:PATH=...' on "
    "the CMake command line or adapt the "
    "MKL_POSSIBLE_PATHS at the beginning of cmake_modules/FindIntelMKL.cmake. "
    "You may also just copy MKL from LSE /home/data/programs/intel/mkl. "
    "Another option is to change CFS_BLAS_LAPACK and CFS_PARDISO to different "
    "implementations. GotoBLAS and ACML should provide similar performance like MKL.")
ELSE(NOT MKL_ROOT_DIR)
  SET(MKL_FOUND 1)
ENDIF(NOT MKL_ROOT_DIR)

SET(MKL_ROOT_DIR ${MKL_ROOT_DIR} CACHE PATH "Directory of MKL." FORCE)
SET(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include")

#-------------------------------------------------------------------------------
# Determine version of MKL 
# - relies (partially) on typical installation path namings for MKL
#-------------------------------------------------------------------------------
If(WIN32)
  STRING(REGEX REPLACE ".*MKL/"
    "" MKL_VERSION
    ${MKL_ROOT_DIR})
ELSE(WIN32)
  STRING(REGEX REPLACE ".*mkl/"
    "" MKL_VERSION10
    ${MKL_ROOT_DIR})
  STRING(REGEX REPLACE ".*Compiler/"
    "" MKL_VERSION11
    ${MKL_ROOT_DIR})
  IF(MKL_VERSION11 STREQUAL "${MKL_ROOT_DIR}")
    SET(MKL_VERSION ${MKL_VERSION10})
  ELSE(MKL_VERSION11 STREQUAL "${MKL_ROOT_DIR}")
    SET(MKL_VERSION ${MKL_VERSION11})
  ENDIF(MKL_VERSION11 STREQUAL "${MKL_ROOT_DIR}")
ENDIF(WIN32)
# MESSAGE("MKL_ROOT_DIR ${MKL_ROOT_DIR}")
# MESSAGE("MKL_VERSION10 ${MKL_VERSION10}")
# MESSAGE("MKL_VERSION11 ${MKL_VERSION11}")
# MESSAGE("MKL_VERSION ${MKL_VERSION}")

STRING(REGEX MATCH "[0-9]*"
  CFS_MKL_VERSION
  ${MKL_VERSION})
# MESSAGE("CFS_MKL_VERSION ${CFS_MKL_VERSION}")


#-------------------------------------------------------------------------------
# Set BLAS_LIBRARY and LAPACK_LIBRARY according to architecture
#-------------------------------------------------------------------------------
IF(CFS_ARCH STREQUAL "I386")
  IF(MKL_FOUND)
    SET(MKL_ARCH_DIR "lib/32")
    SET(LD "${MKL_ROOT_DIR}/${MKL_ARCH_DIR}")
    IF(CFS_DISTRO STREQUAL "MACOSX")
      IF(EXISTS "${MKL_ROOT_DIR}/../../lib/libiomp5.a")
	SET(MKL_LIB_IOMP "${MKL_ROOT_DIR}/../../lib/libiomp5.a")
      ENDIF(EXISTS "${MKL_ROOT_DIR}/../../lib/libiomp5.a")
      IF(EXISTS "${MKL_ROOT_DIR}/../../lib/libguide.a")
	SET(MKL_LIB_GUIDE "${MKL_ROOT_DIR}/../../lib/libguide.a")
      ENDIF(EXISTS "${MKL_ROOT_DIR}/../../lib/libguide.a")
    ENDIF(CFS_DISTRO STREQUAL "MACOSX")
    IF(EXISTS "${MKL_ROOT_DIR}/../lib/ia32/libiomp5.a")
      SET(MKL_LIB_IOMP "${MKL_ROOT_DIR}/../lib/ia32/libiomp5.a")
    ENDIF(EXISTS "${MKL_ROOT_DIR}/../lib/ia32/libiomp5.a")
    IF(EXISTS "${MKL_ROOT_DIR}/../lib/ia32/libguide.a")
      SET(MKL_LIB_GUIDE "${MKL_ROOT_DIR}/../lib/ia32/libguide.a")
    ENDIF(EXISTS "${MKL_ROOT_DIR}/../lib/ia32/libguide.a")
    IF(EXISTS "${LD}/libiomp5.a")
      SET(MKL_LIB_IOMP "${LD}/libiomp5.a")
    ENDIF(EXISTS "${LD}/libiomp5.a")
    IF(EXISTS "${LD}/libguide.a")
      SET(MKL_LIB_GUIDE "${LD}/libguide.a")
    ENDIF(EXISTS "${LD}/libguide.a")

    IF(CFS_MKL_VERSION LESS 10)
      SET(MKL_BLAS_LIB "${LD}/libmkl_ia32.a")
      SET(MKL_LAPACK_LIB "${LD}/libmkl_lapack.a")
      SET(MKL_PARDISO_LIB "${LD}/libmkl_solver.a")
    ELSE(CFS_MKL_VERSION LESS 10)
      IF(CFS_DISTRO STREQUAL "MACOSX")
        SET(MKL_BLAS_LIB "${LD}/libmkl_intel.a;${LD}/libmkl_intel_thread.a;${LD}/libmkl_core.a;${LD}/libmkl_intel_thread.a;${LD}/libmkl_core.a;${LD}/libmkl_intel_thread.a;${LD}/libmkl_core.a")
      ELSE(CFS_DISTRO STREQUAL "MACOSX")
	IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
          SET(MKL_BLAS_LIB "${LD}/libmkl_intel.a;${LD}/libmkl_core.a;${LD}/libmkl_intel_thread.a")
	ELSE(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
          SET(MKL_BLAS_LIB "${LD}/libmkl_gf.a;${LD}/libmkl_core.a;${LD}/libmkl_gnu_thread.a")
	ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC") 
      ENDIF(CFS_DISTRO STREQUAL "MACOSX")
    
      SET(MKL_LAPACK_LIB "")
      SET(MKL_PARDISO_LIB "${LD}/libmkl_solver.a")
    ENDIF(CFS_MKL_VERSION LESS 10) 

  ENDIF(MKL_FOUND)
ENDIF(CFS_ARCH STREQUAL "I386")

IF(CFS_ARCH STREQUAL "X86_64")
  IF(MKL_FOUND)
    SET(MKL_ARCH_DIR "lib/em64t")
    SET(LD "${MKL_ROOT_DIR}/${MKL_ARCH_DIR}")
    IF(EXISTS "${MKL_ROOT_DIR}/../lib/intel64/libiomp5.a")
      SET(MKL_LIB_IOMP "${MKL_ROOT_DIR}/../lib/intel64/libiomp5.a")
    ENDIF(EXISTS "${MKL_ROOT_DIR}/../lib/intel64/libiomp5.a")
    IF(EXISTS "${MKL_ROOT_DIR}/../lib/intel64/libguide.a")
      SET(MKL_LIB_GUIDE "${MKL_ROOT_DIR}/../lib/intel64/libguide.a")
    ENDIF(EXISTS "${MKL_ROOT_DIR}/../lib/intel64/libguide.a")
    IF(EXISTS "${LD}/libiomp5.a")
      SET(MKL_LIB_IOMP "${LD}/libiomp5.a")
    ENDIF(EXISTS "${LD}/libiomp5.a")
    IF(EXISTS "${LD}/libguide.a")
      SET(MKL_LIB_GUIDE "${LD}/libguide.a")
    ENDIF(EXISTS "${LD}/libguide.a")


    IF(CFS_MKL_VERSION LESS 10)
      SET(MKL_BLAS_LIB "${LD}/libmkl_em64t.a")
      SET(MKL_LAPACK_LIB "${LD}/libmkl_lapack.a")
      SET(MKL_PARDISO_LIB "${LD}/libmkl_solver.a")
    ELSE(CFS_MKL_VERSION LESS 10)
      IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
        SET(MKL_BLAS_LIB "${LD}/libmkl_intel_lp64.a;${LD}/libmkl_core.a;${LD}/libmkl_intel_thread.a")
      ELSE(CFS_CXX_COMPILER_NAME STREQUAL "ICC") 
        SET(MKL_BLAS_LIB "${LD}/libmkl_gf_lp64.a;${LD}/libmkl_core.a;${LD}/libmkl_gnu_thread.a")
      ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC") 
    
      SET(MKL_LAPACK_LIB "")
      SET(MKL_PARDISO_LIB "${LD}/libmkl_solver_lp64.a")
    ENDIF(CFS_MKL_VERSION LESS 10)
  ENDIF(MKL_FOUND)
ENDIF(CFS_ARCH STREQUAL "X86_64")

IF(CFS_ARCH STREQUAL "IA64")
  IF(MKL_FOUND)
    SET(MKL_ARCH_DIR "lib/64")
    SET(LD "${MKL_ROOT_DIR}/${MKL_ARCH_DIR}")
    IF(EXISTS "${MKL_ROOT_DIR}/../lib/ia64/libiomp5.a")
      SET(MKL_LIB_IOMP "${MKL_ROOT_DIR}/../lib/ia64/libiomp5.a")
    ENDIF(EXISTS "${MKL_ROOT_DIR}/../lib/ia64/libiomp5.a")
    IF(EXISTS "${MKL_ROOT_DIR}/../lib/ia64/libguide.a")
      SET(MKL_LIB_GUIDE "${MKL_ROOT_DIR}/../lib/ia64/libguide.a")
    ENDIF(EXISTS "${MKL_ROOT_DIR}/../lib/ia64/libguide.a")
    IF(EXISTS "${LD}/libiomp5.a")
      SET(MKL_LIB_IOMP "${LD}/libiomp5.a")
    ENDIF(EXISTS "${LD}/libiomp5.a")
    IF(EXISTS "${LD}/libguide.a")
      SET(MKL_LIB_GUIDE "${LD}/libguide.a")
    ENDIF(EXISTS "${LD}/libguide.a")


    IF(CFS_MKL_VERSION LESS 10)
      SET(MKL_BLAS_LIB "${LD}/libmkl_ipf.a")
      SET(MKL_LAPACK_LIB "${LD}/libmkl_lapack.a")
      SET(MKL_PARDISO_LIB "${LD}/libmkl_solver.a")
    ELSE(CFS_MKL_VERSION LESS 10)
      IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
        SET(MKL_BLAS_LIB "${LD}/libmkl_intel_lp64.a;${LD}/libmkl_core.a;${LD}/libmkl_intel_thread.a")
      ELSE(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
        SET(MKL_BLAS_LIB "${LD}/libmkl_gf_lp64.a;${LD}/libmkl_core.a;${LD}/libmkl_gnu_thread.a")
      ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")

      SET(MKL_LAPACK_LIB "")
      SET(MKL_PARDISO_LIB "${LD}/libmkl_solver_lp64.a")
    ENDIF(CFS_MKL_VERSION LESS 10)


  ENDIF(MKL_FOUND)
ENDIF(CFS_ARCH STREQUAL "IA64")

IF(CFS_MKL_VERSION LESS 10)
  IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
    # Let the Intel compiler decide which guide lib to use
    SET(MKL_BLAS_LIB "${MKL_BLAS_LIB};-lguide;-lpthread")
  ELSE(CFS_CXX_COMPILER_NAME STREQUAL "ICC") 
    # Use guide lib of MKL for GNU compiler
    SET(MKL_BLAS_LIB "${MKL_BLAS_LIB};${MKL_LIB_GUIDE};-lpthread;-lm")
  ENDIF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
ELSE(CFS_MKL_VERSION LESS 10)
  IF(USE_OPENMP AND
     CFS_CXX_COMPILER_NAME STREQUAL "GCC")
    # Using OpenMP and GCC: We should not need libiomp5.a in this case
    SET(MKL_BLAS_LIB "-Wl,--start-group;${MKL_BLAS_LIB};-Wl,--end-group;-lpthread")
  ELSE(USE_OPENMP AND
       CFS_CXX_COMPILER_NAME STREQUAL "GCC") 
    # In all other cases we need libiomp5.a
    IF(CFS_DISTRO STREQUAL "MACOSX")
      SET(MKL_BLAS_LIB "${MKL_BLAS_LIB};${MKL_LIB_IOMP};-lpthread;-lm")
    ELSE(CFS_DISTRO STREQUAL "MACOSX")
      SET(MKL_BLAS_LIB "-Wl,--start-group;${MKL_BLAS_LIB};-Wl,--end-group;${MKL_LIB_IOMP};-lpthread;-lm")
    ENDIF(CFS_DISTRO STREQUAL "MACOSX")
  ENDIF(USE_OPENMP AND
        CFS_CXX_COMPILER_NAME STREQUAL "GCC")
ENDIF(CFS_MKL_VERSION LESS 10)

SET(MKL_BLAS_LIB "${MKL_BLAS_LIB}" CACHE FILEPATH "MKL BLAS.")
IF(CFS_MKL_VERSION LESS 10)
  # LAPACK is a standalone lib in MKL 9
  SET(MKL_LAPACK_LIB "${MKL_LAPACK_LIB}" CACHE FILEPATH "MKL LAPACK.")
ELSE(CFS_MKL_VERSION LESS 10)
  # LAPACK and BLAS are contained in libmkl_core in MKL 10
  SET(MKL_LAPACK_LIB "${MKL_BLAS_LIB}" CACHE FILEPATH "MKL LAPACK.")
ENDIF(CFS_MKL_VERSION LESS 10)
SET(MKL_PARDISO_LIB "${MKL_PARDISO_LIB}" CACHE FILEPATH "MKL Pardiso.")

MARK_AS_ADVANCED(MKL_BLAS_LIB)
MARK_AS_ADVANCED(MKL_LAPACK_LIB)
MARK_AS_ADVANCED(MKL_PARDISO_LIB)
MARK_AS_ADVANCED(MKL_ROOT_DIR)

IF(CFS_BLAS_LAPACK STREQUAL "MKL")
  SET(BLAS_LIBRARY "${MKL_BLAS_LIB}")
  SET(LAPACK_LIBRARY "${MKL_LAPACK_LIB}")
ENDIF(CFS_BLAS_LAPACK STREQUAL "MKL")
SET(PARDISO_LIBRARY "${MKL_PARDISO_LIB}")
