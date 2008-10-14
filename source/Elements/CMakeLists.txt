# Make sure the compiler can find include files from our CFSMESHIO library.
#INCLUDE_DIRECTORIES(${CFSMESHIO_SOURCE_DIR}/include ${HDF5_INCLUDE_DIR})

# Make sure the linker can find the CFSMESHIO library once it is built.
#LINK_DIRECTORIES()

SET(ELEMENTS_SRCS jacobian.cc basefe.cc
  1D/line1fe.cc 1D/line2fe.cc 1D/linefe.cc
  2D/quad1fe.cc 2D/quad2fe.cc 2D/quad9fe.cc 2D/rectanglefe.cc
  2D/triangle1fe.cc 2D/triangle2fe.cc 2D/trianglefe.cc
  3D/hexa1FE.cc 3D/hexa2FE.cc 3D/hexa27FE.cc 3D/hexaFE.cc
  3D/pyra1FE.cc 3D/pyra2FE.cc 3D/pyraFE.cc
  3D/tetra1FE.cc 3D/tetra2FE.cc 3D/tetraFE.cc
  3D/wedge1FE.cc 3D/wedge2FE.cc 3D/wedgeFE.cc)

ADD_LIBRARY(elements STATIC ${ELEMENTS_SRCS})

TARGET_LINK_LIBRARIES(elements
  paramh
  logging
)
