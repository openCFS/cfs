SET(RESULTS_SRCS 
  BaseResults.cc
  ResultFunctor.cc
  ResultInfo.cc
  MHTimeFreqResult.cc
  )


ADD_LIBRARY(results STATIC ${RESULTS_SRCS})

IF(MKL_INCLUDE_DIR)
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
ENDIF(MKL_INCLUDE_DIR)

SET(TARGET_LL
  paramh
  elemmapping
  matvec
  cfsgeneral
  utils
  febasis
#perhaps the following 
#will become necessay in the future
#  bilinforms
#  linforms
  )


TARGET_LINK_LIBRARIES(results
  ${TARGET_LL}
  )
