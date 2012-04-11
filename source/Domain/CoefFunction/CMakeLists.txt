SET(COEFFUNCTION_SRCS
  CoefXpr.cc
  CoefFunction.cc
  CoefFunctionCompound.cc
  CoefFunctionConst.cc
  CoefFunctionExpression.cc
  CoefFunctionGrid.cc
  CoefFunctionTimeFreq.cc
  CoefFunctionSurf.cc
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


