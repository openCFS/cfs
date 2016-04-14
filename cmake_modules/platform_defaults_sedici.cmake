#-----------------------------------------------------------------------------
# These defaults are for sedici cluster at Applied Mechatronics Klagenfurt.
#-----------------------------------------------------------------------------
SET(CFS_BLAS_LAPACK_DEFAULT "MKL")
SET(USE_PARDISO_DEFAULT "ON")

SET(MKL_ROOT_DIR_DEFAULT "/opt/pckg/intel/parallel_studio_xe_2011_sp1_update1/mkl")

# Options to keep in mind for rom and sedici
# cf. http://dberkholz.wordpress.com/2007/10/12/new-gcc-hotness/
# 4.2 -mtune=native -march=native

