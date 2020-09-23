#-----------------------------------------------------------------------------
# These defaults are for Simon Triebenbacher's MacBook Pro.
#-----------------------------------------------------------------------------
SET(USE_TCL_DEFAULT OFF)
SET(USE_PYTHON_DEFAULT ON)

SET(USE_GMSH_DEFAULT ON)

SET(USE_BLAS_LAPACK_DEFAULT "OPENBLAS")
# SET(USE_ILUPACK_DEFAULT ON)

SET(CFS_DEPS_ROOT_DEFAULT "$ENV{HOME}/Documents/dev/CFSDEPS")
SET(CFS_DEPS_CACHE_DIR_DEFAULT "$ENV{HOME}/Documents/dev/CFSDEPSCACHE")

SET(MKL_ROOT_DIR_DEFAULT "/opt/intel/Compiler/11.0/059/Frameworks/mkl")

SET(CFS_INTEL10_OPT_SWITCHES "-axT -parallel")
SET(CFS_INTEL11_OPT_SWITCHES "-xhost -parallel")
