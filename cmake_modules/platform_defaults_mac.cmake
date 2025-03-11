# default Apple setting overwriting platform_defaults.cmake 
# will be overrided e.g. by a .cfs_platform_defaults.cmake in the user's home dire
# we assume an arm64 platform as Intel is more and more outdated 
# to use arm, have -DCMAKE_OSX_ARCHITECTURES=arm64 -DCFS_ARCH=ARM64 in cour cmake call

# Apple's Accelerate framework is very fast but cannot solve complex EV problems. 
# If you need the, use OPENBLAS 
set(USE_BLAS_LAPACK_DEFAULT "ACCELERATE")
set(USE_PARDISO_DEFAULT OFF)

 # To enable IPOPT one needs a HLS license as we assume no MKL. But only a few need ipopt anyway
set(USE_IPOPT_DEFAULT OFF) 

