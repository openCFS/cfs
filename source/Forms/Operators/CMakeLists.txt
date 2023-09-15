SET(OPERATORS_SRCS
  #No cc files yet but will be in the future
  DaeOperator.hh
  )

ADD_LIBRARY(operators STATIC ${OPERATORS_SRCS})


TARGET_LINK_LIBRARIES(operators
  matvec
  #THIS IS NOT GOOD, we should avoid this
  febasis
)