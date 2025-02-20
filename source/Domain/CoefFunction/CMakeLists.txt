SET(COEFFUNCTION_SRCS
  CoefXpr.cc
  CoefFunction.cc
  CoefFunctionAccumulator.cc
  CoefFunctionApprox.cc
  CoefFunctionFormBased.cc
  CoefFunctionCache.cc
  CoefFunctionCompound.cc
  CoefFunctionConst.cc
  CoefFunctionContactForceDensity.cc
  CoefFunctionDistance.cc
  CoefFunctionDummy.cc
  CoefFunctionExpression.cc
  CoefFunctionGrid.cc
  CoefFunctionGridNodal.cc
  CoefFunctionGridNodalSource.cc
  CoefFunctionGridNodalDefault.cc
  CoefFunctionGridNodalInterp.cc
  CoefFunctionGridElem.cc
  CoefFunctionGridElemDefault.cc
  CoefFunctionMeanFlowConvection.cc
  CoefFunctionImpedanceModel.cc
  CoefFunctionMulti.cc
  CoefFunctionOpt.cc
  CoefFunctionFileData.cc
  CoefFunctionTimeFreq.cc
  CoefFunctionSurf.cc
  CoefFunctionPML.cc
  CoefFunctionCurvilinearPML.cc
  CoefFunctionMapping.cc
  CoefFunctionScatteredData.cc
  CoefFunctionStabParams.cc
  CoefFunctionHyst.cc
  CoefFunctionHarmBalance.cc
  CoefFunctionComplexToReal.cc
  CoefFunctionDiagTensorFromScalar.cc
  CoefFunctionMaterialModel.cc
  CoefFunctionConversion.hh
  CoefFunctionSUPG.cc
  CoefFunctionSUPG.hh
)

if(USE_EMBEDDED_PYTHON)
  # see FindPrograms.cmake on how the Python_* stuff is set
  list(APPEND TARGET_LL ${Python_LIBRARIES})
  set(COEFFUNCTION_SRCS ${COEFFUNCTION_SRCS} CoefFunctionPython.cc)
else()
  # empty definitions to please linker
  set(COEFFUNCTION_SRCS ${COEFFUNCTION_SRCS} NoPythonCoefFunctionPython.cc)
endif()

ADD_LIBRARY(coeffunction STATIC ${COEFFUNCTION_SRCS})

IF(TARGET cgal)
  ADD_DEPENDENCIES(coeffunction cgal)
ENDIF()

SET(TARGET_LL 
  scattereddatainout
  materials
  preciceadapter
  ${MUPARSER_LIBRARY} 
  ${LAPACK_LIBRARY} )

IF(USE_FLANN)
  # flann must be linked to because it uses LZ4 (directly included in the FLANN lib)
  LIST(APPEND TARGET_LL ${FLANN_LIBRARY})
ENDIF(USE_FLANN)

TARGET_LINK_LIBRARIES(coeffunction ${TARGET_LL})
