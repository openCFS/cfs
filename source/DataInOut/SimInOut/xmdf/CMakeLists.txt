SET(SIMINOUTXMDF_SRCS
  hdf5io.cc
  simInputXMDF.cc
  simOutputXMDF.cc
  Xmdf.c
  xmdf_private.c
  xmdf_timestep.c
)

ADD_LIBRARY(siminoutxmdf STATIC ${SIMINOUTXMDF_SRCS})

TARGET_LINK_LIBRARIES(siminoutxmdf hdf5_cpp hdf5)

IF(DEBUG)
  SET(CMAKE_C_FLAGS "-g")
ELSE(DEBUG)
  SET(CMAKE_C_FLAGS "")
ENDIF(DEBUG)
