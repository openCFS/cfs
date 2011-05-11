INCLUDE_DIRECTORIES(${VTK_ROOT_DIR}/include/vtk-5.5)

SET(CPLREADER_CGNS_SRCS
  FileReader_CGNS.cc
  )

ADD_LIBRARY(cplreader_cgns STATIC ${CPLREADER_CGNS_SRCS})

TARGET_LINK_LIBRARIES(cplreader_cgns
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtkIO.a
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtkFiltering.a
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtkGraphics.a
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtkCommon.a
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtksys.a
  ${VTK_ROOT_DIR}/lib/vtk-5.5/libvtkzlib.a
  ${CGNS_LIBRARY}
  dl
)

SET_SOURCE_FILES_PROPERTIES(${CPLREADER_CGNS_SRCS}
  PROPERTIES COMPILE_FLAGS "-Wno-deprecated")

