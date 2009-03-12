#-----------------------------------------------------------------------------
# These defaults are for all machines at LSE.
#-----------------------------------------------------------------------------
SET(CFS_BLAS_LAPACK_DEFAULT "MKL")
SET(CFS_PARDISO_DEFAULT "MKL")

SET(USE_ILUPACK_DEFAULT ON)

#-----------------------------------------------------------------------------
# Just Fabian Wein needs SCPIP by default.
#-----------------------------------------------------------------------------
IF(CFS_BUILD_USER STREQUAL "fwein")
    SET(USE_SCPIP_DEFAULT ON)  
ENDIF(CFS_BUILD_USER STREQUAL "fwein")

SET(MKL_ROOT_DIR_DEFAULT "/home/data/programs/intel/mkl/10.0.5.025")

SET(CFS_DEPS_ROOT_DEFAULT "/home/data/libraries/CFSDEPS")
SET(CFS_DEPS_CACHE_DIR_DEFAULT "/home/data/libraries/CFSDEPSCACHE")
