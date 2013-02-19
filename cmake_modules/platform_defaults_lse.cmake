#-----------------------------------------------------------------------------
# These defaults are for all machines at LSE.
#-----------------------------------------------------------------------------
SET(CFS_BLAS_LAPACK_DEFAULT "MKL")
SET(USE_PARDISO_DEFAULT "ON")
SET(CFS_PARDISO_DEFAULT "MKL")

SET(USE_ILUPACK_DEFAULT ON)

#-----------------------------------------------------------------------------
# Just Fabian Wein needs SCPIP by default.
#-----------------------------------------------------------------------------
IF(CFS_BUILD_USER STREQUAL "fwein")
    SET(USE_SCPIP_DEFAULT ON)  
ENDIF(CFS_BUILD_USER STREQUAL "fwein")

SET(MKL_ROOT_DIR_DEFAULT "/home/shareAll/linux_bin/intel/Compiler/11.0/081/mkl")