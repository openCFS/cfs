SET(COEFFUNCTION_SRCS
  CoefFunction.cc
  CoefFunctionExpression.cc
  CoefFunctionTimeFreq.cc
  CoefFunctionGrid.cc
)
#
ADD_LIBRARY(coeffunction STATIC ${COEFFUNCTION_SRCS})
#
#SET(TARGET_LL
#  datainout
#  )
#
#TARGET_LINK_LIBRARIES(coeffunction
#  $(TARGET_LL)
#  )


