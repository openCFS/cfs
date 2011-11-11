SET(RESULTS_SRCS 
  BaseResults.cc
  ResultFunctor.cc
  ResultInfo.cc
  )


ADD_LIBRARY(results STATIC ${RESULTS_SRCS})

SET(TARGET_LL
  paramh
  elemmapping
  matvec
  cfsgeneral
  utils
#perhaps the following 
#will become necessay in the future
#  bilinforms
#  linforms
  )


TARGET_LINK_LIBRARIES(results
  ${TARGET_LL}
  )
