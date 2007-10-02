SET(SIMINOUTHDF5_SRCS
  hdf5io.cc
  simInputHDF5.cc
  simOutputHDF5.cc
)

ADD_LIBRARY(siminouthdf5 STATIC ${SIMINOUTHDF5_SRCS})

TARGET_LINK_LIBRARIES(siminouthdf5 hdf5_cpp hdf5 boost_date_time )

IF(DEBUG)
  SET(CMAKE_C_FLAGS "-g")
ELSE(DEBUG)
  SET(CMAKE_C_FLAGS "")
ENDIF(DEBUG)
