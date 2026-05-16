INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

SET(SIMINOUTHDF5_SRCS
  CommonHDF5.cc
  SimInputHDF5.cc
  SimOutputHDF5.cc)

ADD_LIBRARY(siminouthdf5 STATIC ${SIMINOUTHDF5_SRCS})

add_dependencies(siminouthdf5 hdf5 zlib)

TARGET_LINK_LIBRARIES(siminouthdf5
  mesh
  ${HDF5_LIBRARY} 
  ${CMAKE_DL_LIBS}
  ${BOOST_LIBRARY}  
  ${ZLIB_LIBRARY} )


