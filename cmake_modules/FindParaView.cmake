#-------------------------------------------------------------------------------
# Look for Qt4 directory. Latest ParaView version requires latest Qt4!
#-------------------------------------------------------------------------------
BUILD_EXTLIB("Qt4"
  "${CFS_BINARY_DIR}/qt4"
  "${CFS_DEPS_ROOT}/qt4/build_qt4.pl"
  "build_qt4.log")

#-------------------------------------------------------------------------------
# Look for ParaView directory.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("ParaView with HDF5 support"
  "${CFS_BINARY_DIR}/paraview"
  "${CFS_DEPS_ROOT}/paraview/build_paraview.pl"
  "build_paraview.log")
