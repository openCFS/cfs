#-----------------------------------------------------------------------------
# These defaults are for Simon Triebenbacher's MacBook Pro.
#-----------------------------------------------------------------------------
# SET(USE_ILUPACK_DEFAULT ON)

SET(CFS_DEPS_CACHE_DIR_DEFAULT "$ENV{HOME}/Documents/dev/CFSDEPSCACHE")

SET(CFS_BLAS_LAPACK_DEFAULT "MKL")
SET(CFS_PARDISO_DEFAULT "MKL")
SET(MKL_ROOT_DIR_DEFAULT "/opt/intel/Compiler/11.0/059/Frameworks/mkl")

SET(CFS_INTEL10_OPT_SWITCHES "-axT -parallel")
SET(CFS_INTEL11_OPT_SWITCHES "-xhost -parallel")
