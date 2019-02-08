SET(SCATTEREDDATAINOUT_SRCS
  ScatteredDataReader.cc
  ScatteredDataReaderCSV.cc
  ScatteredDataReaderCSVT.cc
)

SET(TARGET_LL
  ${BOOST_SYSTEM_LIB}
  )

IF(USE_CCMIO)
  LIST(APPEND
    SCATTEREDDATAINOUT_SRCS  
    ScatteredDataReaderCCM.cc
    )

  LIST(APPEND TARGET_LL
    ${CCMIO_LIBRARY}
    ${CGNS_LIBRARY}
    ${HDF5_LIBRARY}
    ${ZLIB_LIBRARY}
    )
ENDIF()

ADD_LIBRARY(scattereddatainout STATIC ${SCATTEREDDATAINOUT_SRCS})

ADD_DEPENDENCIES(scattereddatainout boost)

TARGET_LINK_LIBRARIES(scattereddatainout
  ${TARGET_LL}
  )
