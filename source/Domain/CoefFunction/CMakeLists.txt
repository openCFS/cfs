SET(COEFFUNCTION_SRCS
  CoefXpr.cc
  CoefFunction.cc
  CoefFunctionAccumulator.cc
  CoefFunctionApprox.cc
  CoefFunctionFormBased.cc
  CoefFunctionCompound.cc
  CoefFunctionConst.cc
  CoefFunctionExpression.cc
  CoefFunctionGrid.cc
  CoefFunctionGridNodal.cc
  CoefFunctionGridNodalDefault.cc
  CoefFunctionGridNodalInterp.cc
  CoefFunctionMeanFlowConvection.cc
  CoefFunctionImpedanceModel.cc
  CoefFunctionMulti.cc
  CoefFunctionOpt.cc
  CoefFunctionTimeFreq.cc
  CoefFunctionSurf.cc
  CoefFunctionPML.cc
  CoefFunctionScatteredData.cc
  CoefFunctionStabParams.cc
  
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

SET(TARGET_LL
  scattereddatainout
  ${MUPARSER_LIBRARY}
  )

TARGET_LINK_LIBRARIES(coeffunction
  ${TARGET_LL}
  )
