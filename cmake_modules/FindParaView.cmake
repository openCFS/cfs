#-------------------------------------------------------------------------------
# Look for ParaView directory.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("ParaView with HDF5 support"
  "${CFS_BINARY_DIR}/paraview"
  "${CFS_DEPS_ROOT}/paraview/build_paraview.pl"
  "build_paraview.log")
