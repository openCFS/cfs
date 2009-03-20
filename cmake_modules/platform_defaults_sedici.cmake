#-----------------------------------------------------------------------------
# These defaults are for sedici cluster at Applied Mechatronics Klagenfurt.
#-----------------------------------------------------------------------------
SET(CFS_BLAS_LAPACK_DEFAULT "MKL")
SET(CFS_PARDISO_DEFAULT "MKL")

SET(MKL_ROOT_DIR_DEFAULT "/opt/intel/Compiler/11.0/081/mkl")

SET(CFS_DEPS_CACHE_DIR_DEFAULT "/opt/CFSDEPSCACHE")

SET(CFS_INTEL11_OPT_SWITCHES "-xhost -parallel")

# Options to keep in mind for rom and sedici
# cf. http://dberkholz.wordpress.com/2007/10/12/new-gcc-hotness/
# 4.2 -mtune=native -march=native
# 4.3 x86_64 acml -mveclibabi=acml
