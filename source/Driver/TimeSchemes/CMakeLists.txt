SET(TIMESCHEME_SRCS
  TimeSchemeGLM.cc
  GLMSchemeLib.cc
  AdaptiveTimesteppingData.cc
  )

ADD_LIBRARY(timeschemes STATIC ${TIMESCHEME_SRCS})

SET(TARGET_LL
  datainout
  logging
  matvec
  cfsgeneral
)

TARGET_LINK_LIBRARIES(timeschemes ${TARGET_LL})
