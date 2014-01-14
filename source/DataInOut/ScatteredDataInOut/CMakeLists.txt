SET(SCATTEREDDATAINOUT_SRCS
  ScatteredDataReaderCSV.cc
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

TARGET_LINK_LIBRARIES(scattereddatainout
  ${TARGET_LL}
  )

IF(USE_CCMIO)
  ADD_EXECUTABLE(ScatteredDataReaderCCMDemo ScatteredDataReaderCCMDemo.cc)
  TARGET_LINK_LIBRARIES(ScatteredDataReaderCCMDemo 
    scattereddatainout
    ${TARGET_LL}
    cfsgeneral)
ENDIF()
