SET(VTK_ROOT_DIR
  "/home/lse24/lse/strieben/packages/vtk_freiberg_static_nightly_230508_${CFS_ARCH}"
  CACHE
  PATH
  "Path to VTK 5.3 directory.")
  
MARK_AS_ADVANCED(VTK_ROOT_DIR)

INCLUDE_DIRECTORIES(${VTK_ROOT_DIR}/include/vtk-5.3)

SET(CPLREADER_OPENFOAM_SRCS
  FileReader_OPENFOAM.cc
  )

ADD_LIBRARY(cplreader_openfoam STATIC ${CPLREADER_OPENFOAM_SRCS})

TARGET_LINK_LIBRARIES(cplreader_openfoam
  ${VTK_ROOT_DIR}/lib/vtk-5.3/libvtkIO.a
  ${VTK_ROOT_DIR}/lib/vtk-5.3/libvtkFiltering.a
  ${VTK_ROOT_DIR}/lib/vtk-5.3/libvtkGraphics.a
  ${VTK_ROOT_DIR}/lib/vtk-5.3/libvtkCommon.a
  ${VTK_ROOT_DIR}/lib/vtk-5.3/libvtksys.a
  ${VTK_ROOT_DIR}/lib/vtk-5.3/libvtkzlib.a
  dl
)
