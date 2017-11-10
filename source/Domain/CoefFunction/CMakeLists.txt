SET(COEFFUNCTION_SRCS
  CoefXpr.cc
  CoefFunction.cc
  CoefFunctionAccumulator.cc
  CoefFunctionApprox.cc
  CoefFunctionFormBased.cc
  CoefFunctionCache.cc
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
  CoefFunctionMapping.cc
  CoefFunctionScatteredData.cc
  CoefFunctionStabParams.cc
  CoefFunctionHyst.cc
)

ADD_LIBRARY(coeffunction STATIC ${COEFFUNCTION_SRCS})

SET(TARGET_LL 
  scattereddatainout
  ${MUPARSER_LIBRARY} 
  ${LAPACK_LIBRARY} )

IF(USE_FLANN)
  # flann must be linked to because it uses LZ4 (directly included in the FLANN lib)
  LIST(APPEND TARGET_LL ${FLANN_LIBRARY})
ENDIF(USE_FLANN)

TARGET_LINK_LIBRARIES(coeffunction ${TARGET_LL})
