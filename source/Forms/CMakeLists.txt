SET(INTEGRATION_SRCS
  IntScheme.cc 
  )

ADD_LIBRARY(integration STATIC ${INTEGRATION_SRCS})


TARGET_LINK_LIBRARIES(integration
  elemmapping
)
