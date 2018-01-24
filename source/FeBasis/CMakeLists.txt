SET(FEBASIS_SRC 
  BaseFE.cc
  FeSpace.cc
  FeSpaceHi.cc
  FeSpaceNodal.cc
  FeFunctions.cc
  FeHi.cc
  FeNodal.cc
  H1/FeSpaceH1Hi.cc
  H1/FeSpaceH1Nodal.cc
  HCurl/FeSpaceHCurlHi.cc
  H1/H1Elems.cc
  H1/H1ElemsLagExpl.cc
  H1/H1ElemsLagVar.cc
  H1/H1ElemsHi.cc
  HCurl/HCurlElemsHi.cc
  HCurl/HCurlElems.cc
  L2/FeSpaceL2.cc
  L2/FeSpaceL2Nodal.cc
  FeSpaceConst.cc
)

ADD_LIBRARY(febasis STATIC ${FEBASIS_SRC})

TARGET_LINK_LIBRARIES(febasis
  elemmapping
)

IF(USE_ADOLC)
  INCLUDE_DIRECTORIES(${ADOLC_INCLUDE_DIR})
  #LINK_DIRECTORIES(${ADOLC_LIBRARY})
ENDIF(USE_ADOLC)


