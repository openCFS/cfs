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
  CoefFunctionMulti.cc
  CoefFunctionTimeFreq.cc
  CoefFunctionSurf.cc
  CoefFunctionPML.cc
  CoefFunctionScatteredData.cc
  CoefFunctionStabParams.cc
  CoefFunctionMeanFlowConvection.cc
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

IF(USE_CCMIO)
  ADD_EXECUTABLE(readexample readexample.cc)
  TARGET_LINK_LIBRARIES(readexample
    ${CCMIO_LIBRARY}
    ${CGNS_LIBRARY}
    ${HDF5_LIBRARY}
    ${ZLIB_LIBRARY})
ENDIF()
