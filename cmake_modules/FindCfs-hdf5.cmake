#-------------------------------------------------------------------------------
# Look for cfs-hdf5 header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("cfs-hdf5"
  "${CFS_BINARY_DIR}/include/hdf5Common.h"
  "${CFS_DEPS_ROOT}/cfs-hdf5/build_cfs-hdf5.pl"
  "build_cfs-hdf5.log")

SET(CFS_HDF5_LIB 
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libcfs-hdf5.a"
)
