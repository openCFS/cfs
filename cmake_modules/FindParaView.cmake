#-------------------------------------------------------------------------------
# Look for ParaView directory.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("ParaView 3.6.1 with HDF5 support"
  "${CFS_BINARY_DIR}/paraview"
  "${CFS_DEPS_ROOT}/paraview/build_paraview.pl"
  "build_paraview.log")
