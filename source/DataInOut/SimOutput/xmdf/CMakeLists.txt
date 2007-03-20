SET(SIMOUTXMDF_SRCS
  simOutputXMDF.cc
  Xmdf.c
  xmdf_private.c
  xmdf_timestep.c
)

ADD_LIBRARY(simoutxmdf STATIC ${SIMOUTXMDF_SRCS})

IF(DEBUG)
  SET(CMAKE_C_FLAGS "-g")
ELSE(DEBUG)
  SET(CMAKE_C_FLAGS "")
ENDIF(DEBUG)

TARGET_LINK_LIBRARIES(simoutxmdf hdf5_cpp hdf5)