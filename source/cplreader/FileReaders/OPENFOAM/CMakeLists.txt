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

SET_SOURCE_FILES_PROPERTIES(${CPLREADER_OPENFOAM_SRCS}
  PROPERTIES COMPILE_FLAGS "-Wno-deprecated")

