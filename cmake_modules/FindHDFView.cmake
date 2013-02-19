#-------------------------------------------------------------------------------
# Look for HDFView directory.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("HDFViewer"
  "${CFS_BINARY_DIR}/hdfview"
  "${CFS_DEPS_ROOT}/hdfview/build_hdfview.pl"
  "build_hdfview.log")
