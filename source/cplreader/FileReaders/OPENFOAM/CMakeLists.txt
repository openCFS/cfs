INCLUDE_DIRECTORIES(${VTK_ROOT_DIR}/include/vtk-5.5)

SET(CPLREADER_OPENFOAM_SRCS
  FileReader_OPENFOAM.cc
  )

ADD_LIBRARY(cplreader_openfoam STATIC ${CPLREADER_OPENFOAM_SRCS})

TARGET_LINK_LIBRARIES(cplreader_openfoam
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtkIO.a
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtkFiltering.a
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtkGraphics.a
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtkCommon.a
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtksys.a
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtkzlib.a
  dl
)

SET_SOURCE_FILES_PROPERTIES(${CPLREADER_OPENFOAM_SRCS}
  PROPERTIES COMPILE_FLAGS "-Wno-deprecated")

