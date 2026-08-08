SET(LINFORMS_SRCS
 SingleEntryInt.cc
 ContactLinInt.cc
)

ADD_LIBRARY(linforms STATIC ${LINFORMS_SRCS})


TARGET_LINK_LIBRARIES(linforms
  domain
  mathparser
  febasis
)
