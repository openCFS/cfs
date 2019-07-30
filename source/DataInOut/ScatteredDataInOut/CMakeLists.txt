SET(SCATTEREDDATAINOUT_SRCS
  ScatteredDataReader.cc
  ScatteredDataReaderCSV.cc
  ScatteredDataReaderCSVT.cc
)

SET(TARGET_LL
  ${BOOST_SYSTEM_LIB}
  )

ADD_LIBRARY(scattereddatainout STATIC ${SCATTEREDDATAINOUT_SRCS})

ADD_DEPENDENCIES(scattereddatainout boost)

TARGET_LINK_LIBRARIES(scattereddatainout
  ${TARGET_LL}
  )
