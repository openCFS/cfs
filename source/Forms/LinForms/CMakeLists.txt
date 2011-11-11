SET(LINFORMS_SRCS
  #nothing right now, each cc file
  #is included directly in the header 
  )

ADD_LIBRARY(linforms STATIC ${FORMS_SRCS})


TARGET_LINK_LIBRARIES(linforms
  domain
  mathparser
  febasis
)
