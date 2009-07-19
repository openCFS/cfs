#-----------------------------------------------------------------------------
# These defaults are for Bastian's machine
#-----------------------------------------------------------------------------

# note compilers have to be set using CC=* and CXX=* (when using intel compiler)

SET(MKL_ROOT_DIR_DEFAULT "/opt/intel/Compiler/11.0/074/mkl")
SET(CFS_DEPS_ROOT_DEFAULT "/homelocal/bastian/workspace/cfs/cfsdeps")
SET(CFS_DEPS_CACHE_DIR_DEFAULT "/homelocal/bastian/workspace/cfs/cfsdeps.cache")

SET(CFS_BLAS_LAPACK_DEFAULT "MKL")
SET(CFS_PARDISO_DEFAULT "MKL")

#SET(USE_ILUPACK_DEFAULT ON)

SET(USE_SCPIP_DEFAULT ON)

IF(CFS_BINARY_DIR MATCHES /release$)
  SET(TESTSUITE_DIR "/homelocal/bastian/workspace/cfs/cfstest" CACHE STRING
    "Path where the CFS++ testsuite is located.")
ENDIF(CFS_BINARY_DIR MATCHES /release$)

IF(CFS_BINARY_DIR MATCHES /icc$)
  SET(TESTSUITE_DIR "/homelocal/bastian/workspace/cfs/cfstest" CACHE STRING
    "Path where the CFS++ testsuite is located.")
ENDIF(CFS_BINARY_DIR MATCHES /icc$)

IF(CFS_BINARY_DIR MATCHES /debug$)
  SET(DEBUG_DEFAULT ON)
ENDIF(CFS_BINARY_DIR MATCHES /debug$)

IF(CFS_BINARY_DIR MATCHES /openmp$)
  SET(TESTSUITE_DIR "/homelocal/bastian/workspace/cfs/cfstest" CACHE STRING
    "Path where the CFS++ testsuite is located.")
  SET(USE_OPENMP_DEFAULT ON)
ENDIF(CFS_BINARY_DIR MATCHES /openmp$)

IF(CFS_BINARY_DIR MATCHES /iccomp$)
  SET(TESTSUITE_DIR "/homelocal/bastian/workspace/cfs/cfstest" CACHE STRING
    "Path where the CFS++ testsuite is located.")
  SET(USE_OPENMP_DEFAULT ON)
ENDIF(CFS_BINARY_DIR MATCHES /iccomp$)

IF(CFS_BINARY_DIR MATCHES /nomkl$)
  SET(TESTSUITE_DIR "/homelocal/bastian/workspace/cfs/cfstest" CACHE STRING
    "Path where the CFS++ testsuite is located.")
  SET(CFS_PARDISO_DEFAULT "SCHENK")
  SET(USE_IPOPT_DEFAULT ON)
ENDIF(CFS_BINARY_DIR MATCHES /nomkl$)


