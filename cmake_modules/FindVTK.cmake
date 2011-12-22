SET(VTK_FOUND 0)

SET(VTK_VERSION "5.9")

#-------------------------------------------------------------------------------
# Look for VTK Customization header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("VTK ${VTK_VERSION}"
  "${CFS_BINARY_DIR}/vtk"
  "${CFS_DEPS_ROOT}/paraview/build_vtk.pl"
  "build_vtk.log")

#-------------------------------------------------------------------------------
# Determine paths of VTK libraries.
#-------------------------------------------------------------------------------
SET(VTK_ROOT_DIR
  "${CFS_BINARY_DIR}/vtk"
  CACHE
  PATH
  "Path to VTK ${VTK_VERSION} directory.")
  
MARK_AS_ADVANCED(VTK_ROOT_DIR)

SET(VTK_FOUND 1)
