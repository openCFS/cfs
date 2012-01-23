SET(COEFFUNCTION_SRCS
  CoefFunction.cc
  CoefFunctionExpression.cc
  CoefFunctionTimeFreq.cc
)
#
ADD_LIBRARY(coeffunction STATIC ${COEFFUNCTION_SRCS})
#
#SET(TARGET_LL
#  elemmapping
#  matvec
#  mathparser
#  )
#
#TARGET_LINK_LIBRARIES(coeffunction
#  $(TARGET_LL)
#  )
#

