SET(SIMOUTGID_SRCS SimOutGiD.cc)

ADD_LIBRARY(simoutgid STATIC ${SIMOUTGID_SRCS})

# make sure cfsdeps are available (could be more specific, but this needs IFs)
ADD_DEPENDENCIES(simoutgid cfsdeps)

TARGET_LINK_LIBRARIES(simoutgid
  ${GIDPOST_LIBRARY}
  ${HDF5_LIBRARY}
  ${ZLIB_LIBRARY}
  )
