#-----------------------------------------------------------------------------
# These defaults are for sedici cluster at Applied Mechatronics Klagenfurt.
#-----------------------------------------------------------------------------
SET(CFS_BLAS_LAPACK_DEFAULT "MKL")
SET(CFS_PARDISO_DEFAULT "MKL")

SET(USE_TCL_DEFAULT OFF)
SET(USE_PYTHON_DEFAULT ON)

SET(MKL_ROOT_DIR_DEFAULT "/opt/pckg/intel/Compiler/11.0/081/mkl")

#SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -mtune=native")
#SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -march=native")
#SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -combine")
#SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -fprefetch-loop-arrays")
#SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -funroll-all-loops")
#SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -ftree-parallelize-loops=4")
#SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -fwhole-program")
# SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -mabm")
# If ACML is activated the following option may be used
# SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -mveclibabi=acml")
SET(CFS_GCC43_OPT_SWITCHES "${CFS_GCC43_OPT_SWITCHES} -Wno-unused")

SET(CFS_INTEL11_OPT_SWITCHES "${CFS_INTEL11_OPT_SWITCHES} -xhost")
SET(CFS_INTEL11_OPT_SWITCHES "${CFS_INTEL11_OPT_SWITCHES} -parallel")

# Options to keep in mind for rom and sedici
# cf. http://dberkholz.wordpress.com/2007/10/12/new-gcc-hotness/
# 4.2 -mtune=native -march=native
# 4.3 x86_64 acml -mveclibabi=acml
