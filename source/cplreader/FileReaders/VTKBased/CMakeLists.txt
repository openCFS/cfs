INCLUDE_DIRECTORIES(${VTK_ROOT_DIR}/include/vtk-5.8)

SET(TARGET_LL
  ${VTK_ROOT_DIR}/lib/vtk-5.8/libvtkIO.a
  ${VTK_ROOT_DIR}/lib/vtk-5.8/libvtkFiltering.a
  ${VTK_ROOT_DIR}/lib/vtk-5.8/libvtkGraphics.a
  ${VTK_ROOT_DIR}/lib/vtk-5.8/libvtkCommon.a
  ${VTK_ROOT_DIR}/lib/vtk-5.8/libvtksys.a
  ${VTK_ROOT_DIR}/lib/vtk-5.8/libvtkzlib.a
  dl
)

IF(CPLREADER_OPENFOAM)
  SUBDIRS(OPENFOAM)
ENDIF(CPLREADER_OPENFOAM)

IF(CPLREADER_ENSIGHT)
  SUBDIRS(EnSight)
ENDIF(CPLREADER_ENSIGHT)
