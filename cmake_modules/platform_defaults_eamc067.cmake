#-----------------------------------------------------------------------------
# These defaults are for Bastian's machine
#-----------------------------------------------------------------------------

# note compilers have to be set using CC=* and CXX=* (when using intel compiler)

SET(MKL_ROOT_DIR_DEFAULT "/opt/intel/Compiler/11.0/074/mkl")
SET(CFS_DEPS_ROOT_DEFAULT "/homelocal/bastian/workspace/cfs/deps")
SET(CFS_DEPS_CACHE_DIR_DEFAULT "/homelocal/bastian/workspace/cfs/deps.cache")

SET(CFS_BLAS_LAPACK_DEFAULT "MKL")
SET(CFS_PARDISO_DEFAULT "MKL")

SET(USE_ILUPACK_DEFAULT ON)

SET(USE_SCPIP_DEFAULT ON)

SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -mtune=native -march=native")

# release build is default

IF(CFS_BINARY_DIR MATCHES /debug$)
  SET(DEBUG_DEFAULT ON)
ELSE(CFS_BINARY_DIR_MATCHES /debug$)
  SET(TESTSUITE_DIR "/homelocal/bastian/workspace/cfs/test" CACHE STRING
    "Path where the CFS++ testsuite is located.")
ENDIF(CFS_BINARY_DIR MATCHES /debug$)

IF(CFS_BINARY_DIR MATCHES /nomkl$)
  SET(CFS_PARDISO_DEFAULT "SCHENK")
  SET(USE_IPOPT_DEFAULT ON)
ENDIF(CFS_BINARY_DIR MATCHES /nomkl$)


