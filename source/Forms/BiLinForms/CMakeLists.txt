SET(BILINFORMS_SRCS
  #nothing right now, each cc file
  #is included directly in the header 
  )

ADD_LIBRARY(bilinforms STATIC ${BILINFORMS_SRCS})


TARGET_LINK_LIBRARIES(bilinforms
  domain
  mathparser
  febasis
)
